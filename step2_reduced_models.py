"""
Step 2: reduce the physics-based cell models to fast equivalent-circuit models (ECMs).

For each chemistry (Li-ion: Chen2020 NMC811/graphite-SiOx in PyBaMM SPMe; Na-ion:
Chayambuka2022 hard carbon/NVPF in PyBaMM DFN), a virtual pulse test is run on the
physics-based model:

  - a C/20 discharge gives a pseudo open-circuit voltage curve and the usable capacity;
  - 2C pulses of 10 s followed by 60 s rest, every 10% SOC, give R0, R1 and tau per SOC.

The Li-ion ECM (OCV(SOC) - R0*I - V1, with one RC branch) is then checked against the
physics-based model on a dynamic current profile derived from the WLTC drive cycle. The Na-ion
pulse parameters are used only to compare the two chemistries as power devices.
Resistances are stored normalised by capacity (ohm*Ah), so packs of any size can be built.

Usage:  python step2_reduced_models.py        (about 1 min)
"""
import json
import os

import numpy as np
import pandas as pd
import pybamm
from scipy.optimize import least_squares

import plotstyle as ps
from ecm import ECMCell

HERE = os.path.dirname(os.path.abspath(__file__))
CELLS = {
    "Li-ion": dict(model=pybamm.lithium_ion.SPMe, params="Chen2020", label="Li-ion (Chen2020, SPMe)"),
    "Na-ion": dict(model=pybamm.sodium_ion.BasicDFN, params="Chayambuka2022", label="Na-ion (Chayambuka2022, DFN)"),
}
SOC_POINTS = np.arange(0.9, 0.05, -0.1)


def solve(cell, experiment=None, current=None, t_end=None):
    p = pybamm.ParameterValues(cell["params"])
    if current is not None:
        p["Current function [A]"] = pybamm.Interpolant(current[0], current[1], pybamm.t)
        sim = pybamm.Simulation(cell["model"](), parameter_values=p)
        return sim.solve(np.arange(0, t_end, 1.0)), p
    sim = pybamm.Simulation(cell["model"](), parameter_values=p, experiment=pybamm.Experiment(experiment, period="1 second"))
    return sim.solve(), p


def extract(name, cell):
    # pseudo-OCV and capacity from a C/20 discharge
    sol, p = solve(cell, ["Discharge at C/20 until {} V".format(pybamm.ParameterValues(cell["params"])["Lower voltage cut-off [V]"])])
    q = sol["Discharge capacity [A.h]"].entries
    Q = float(q[-1])
    ocv_soc = 1 - q / Q
    ocv_v = sol["Voltage [V]"].entries
    # pulse test
    Qn = p["Nominal cell capacity [A.h]"]
    steps = []
    for k in range(len(SOC_POINTS)):
        steps += [f"Discharge at {0.2 * Q / Qn:.4f}C for 30 minutes",      # 10% of the usable capacity
                  "Rest for 30 minutes", "Discharge at 2C for 10 seconds", "Rest for 60 seconds"]
    sol, _ = solve(cell, steps)
    rows = []
    n_ok = len(sol.cycles) // 4                # stops early if a pulse hits the voltage cut-off
    for k in range(n_ok):
        c = sol.cycles[4 * k + 2]           # pulse
        r = sol.cycles[4 * k + 3]           # relaxation
        pre = sol.cycles[4 * k + 1]["Voltage [V]"].entries[-1]
        I = float(np.mean(c["Current [A]"].entries[1:]))
        v_p = c["Voltage [V]"].entries
        R0 = (pre - v_p[1]) / I if len(v_p) > 1 else np.nan  # instant drop at pulse start
        t_r = r["Time [s]"].entries - r["Time [s]"].entries[0]
        v_r = r["Voltage [V]"].entries

        def resid(x):
            R1, tau = np.exp(x)
            v1_0 = I * R1 * (1 - np.exp(-10 / tau))  # RC voltage at the end of the 10 s pulse
            return pre - v1_0 * np.exp(-t_r / tau) - v_r
        fit = least_squares(resid, x0=np.log([R0, 20.0]), bounds=(np.log([R0 * 1e-3, 0.5]), np.log([R0 * 100, 2000])))
        R1, tau = np.exp(fit.x)
        soc = 1 - float(sol.cycles[4 * k + 1]["Discharge capacity [A.h]"].entries[-1]) / Q
        rows.append(dict(soc=soc, R0_ohm_Ah=R0 * Q, R1_ohm_Ah=R1 * Q, tau_s=tau))
    tab = pd.DataFrame(rows).sort_values("soc")
    # remove the C/20 polarisation from the pseudo-OCV with the fitted resistances
    i_c20 = Q / 20
    r_tot = (np.interp(ocv_soc, tab.soc, tab.R0_ohm_Ah) + np.interp(ocv_soc, tab.soc, tab.R1_ohm_Ah)) / Q
    ocv_v = ocv_v + i_c20 * r_tot
    model = dict(name=name, label=cell["label"], capacity_Ah=Q,
                 ocv_soc=ocv_soc[::-1][::20].tolist(), ocv_v=ocv_v[::-1][::20].tolist(),  # ascending SOC, thinned
                 soc=tab.soc.tolist(), R0_ohm_Ah=tab.R0_ohm_Ah.tolist(), R1_ohm_Ah=tab.R1_ohm_Ah.tolist(),
                 tau_s=tab.tau_s.tolist(), v_min=float(p["Lower voltage cut-off [V]"]),
                 v_max=float(p["Upper voltage cut-off [V]"]))
    return model


def drive_current(Q, seconds=1800):
    """Cell current from the WLTC power profile, scaled so that the peak discharge is 2C."""
    w = pd.read_csv(pybamm.DataLoader().get_data("WLTC.csv"), comment="#", header=None).values
    t = np.arange(0, seconds + 1, 1.0)
    p = np.interp(t, w[:, 0], w[:, 1])
    return t, p / p.max() * 2 * Q


def main():
    os.makedirs(os.path.join(HERE, "models"), exist_ok=True)
    rows, val = [], {}
    for name, cell in CELLS.items():
        m = extract(name, cell)
        with open(os.path.join(HERE, "models", f"ecm_{name.replace('-', '').lower()}.json"), "w") as fh:
            json.dump(m, fh, indent=1)
        if name != "Li-ion":
            # only the Li-ion ECM is used downstream; the Na-ion parameters serve the technology screening
            rows.append(dict(cell=cell["label"], capacity_Ah=round(m["capacity_Ah"], 4),
                             R0_mid_ohm_Ah=round(float(np.interp(0.5, m["soc"], m["R0_ohm_Ah"])), 4),
                             R1_mid_ohm_Ah=round(float(np.interp(0.5, m["soc"], m["R1_ohm_Ah"])), 4),
                             tau_mid_s=round(float(np.interp(0.5, m["soc"], m["tau_s"])), 1),
                             ecm_vs_physics_rmse_mV=np.nan))
            continue
        # validation on a dynamic profile, starting from 80% SOC
        t, i = drive_current(m["capacity_Ah"])
        p = pybamm.ParameterValues(cell["params"])
        p["Current function [A]"] = pybamm.Interpolant(t, i, pybamm.t)
        sim = pybamm.Simulation(cell["model"](), parameter_values=p)
        sol = sim.solve(t, initial_soc=0.8)
        t = t[t <= sol["Time [s]"].entries[-1]]
        v_phys = sol["Voltage [V]"](t)
        # start the ECM at the state of charge whose OCV matches the model's initial rest voltage
        soc0 = float(np.interp(v_phys[0], m["ocv_v"], m["ocv_soc"]))
        e = ECMCell(m, soc0=soc0)
        v_ecm = np.array([e.step_current(x, 1.0) for x in i])[: len(v_phys)]
        rmse = float(np.sqrt(np.mean((v_ecm - v_phys) ** 2)) * 1000)
        rows.append(dict(cell=cell["label"], capacity_Ah=round(m["capacity_Ah"], 4),
                         R0_mid_ohm_Ah=round(float(np.interp(0.5, m["soc"], m["R0_ohm_Ah"])), 4),
                         R1_mid_ohm_Ah=round(float(np.interp(0.5, m["soc"], m["R1_ohm_Ah"])), 4),
                         tau_mid_s=round(float(np.interp(0.5, m["soc"], m["tau_s"])), 1),
                         ecm_vs_physics_rmse_mV=round(rmse, 1)))
        val[name] = (t[: len(v_phys)], v_phys, v_ecm, m)
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(HERE, "results", "step2_ecm_parameters.csv"), index=False)
    print(df.to_string(index=False))

    t, vp, ve, m = val["Li-ion"]
    fig, ax = ps.subplots(figsize=(5.5, 2.6))
    ax.plot(t / 60, vp, color=ps.C["black"], lw=1, label="Physics-based SPMe (PyBaMM)")
    ax.plot(t / 60, ve, color=ps.C["orange"], lw=1, ls="--", label="Reduced ECM")
    ax.set(xlabel="Time (min)", ylabel="Cell voltage (V)")
    ax.legend(frameon=False, loc="lower left", bbox_to_anchor=(0.0, 1.0), ncol=2)
    ps.save(fig, os.path.join(HERE, "figures", "fig2_ecm_validation.png"))


if __name__ == "__main__":
    main()
