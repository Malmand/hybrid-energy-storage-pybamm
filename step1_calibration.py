"""
Step 1: calibrate and check a physics-based Li-ion model against experimental data.

Data: constant-current discharges at 1C and 5C of the Kokam pouch cell measured by Ecker et al.
(J. Electrochem. Soc. 162 (2015) A1836), distributed with PyBaMM (pybamm.DataLoader).
Model: PyBaMM SPMe with the Ecker2015 parameter set.

Two-stage calibration, each stage using the data that carries the information:
  1. capacity balance from the 1C curve (close to equilibrium): negative electrode active
     material volume fraction;
  2. transport and kinetics from the 5C curve: negative and positive particle diffusivities and
     the exchange-current density (log-scaled multipliers).
Each stage is then checked on the other curve, so a gain on one rate does not hide a loss on the other.

Usage:  python step1_calibration.py        (about 1 min)
"""
import os

import numpy as np
import pandas as pd
import pybamm
from scipy.optimize import least_squares, minimize_scalar

import plotstyle as ps

HERE = os.path.dirname(os.path.abspath(__file__))
NAMES = ["Negative particle diffusivity [m2.s-1]", "Positive particle diffusivity [m2.s-1]",
         "Negative electrode exchange-current density [A.m-2]", "Positive electrode exchange-current density [A.m-2]"]
AM = "Negative electrode active material volume fraction"


def load(name):
    return pd.read_csv(pybamm.DataLoader().get_data(name), header=None).values


def parameters(am_scale=1.0, x=(0, 0, 0)):
    base = pybamm.ParameterValues("Ecker2015")
    p = base.copy()
    a = 10 ** np.asarray(x, float)
    scale = lambda f, k: (lambda *args: k * f(*args))
    p.update({AM: base[AM] * am_scale,
              NAMES[0]: scale(base[NAMES[0]], a[0]), NAMES[1]: scale(base[NAMES[1]], a[1]),
              NAMES[2]: scale(base[NAMES[2]], a[2]), NAMES[3]: scale(base[NAMES[3]], a[2])})
    return p


def simulate(p, c_rate):
    sim = pybamm.Simulation(pybamm.lithium_ion.SPMe(), parameter_values=p,
                            experiment=pybamm.Experiment([f"Discharge at {c_rate}C until 2.5 V"]))
    sol = sim.solve()
    return sol["Time [s]"].entries, sol["Voltage [V]"].entries


def residual(p, data, c_rate):
    try:
        t, v = simulate(p, c_rate)
    except Exception:
        return np.ones(len(data))
    return np.interp(data[:, 0], t, v, right=2.5) - data[:, 1]


def rmse_mv(p, data, c_rate):
    return float(np.sqrt(np.mean(residual(p, data, c_rate) ** 2)) * 1000)


def main():
    d1, d5 = load("Ecker_1C.csv"), load("Ecker_5C.csv")
    stage1 = minimize_scalar(lambda s: rmse_mv(parameters(s), d1, 1), bounds=(0.9, 1.1),
                             method="bounded", options={"xatol": 1e-3})
    am = stage1.x
    stage2 = least_squares(lambda x: residual(parameters(am, x), d5, 5), x0=[0, 0, 0],
                           bounds=([-1.5] * 3, [1.5] * 3), diff_step=0.05, max_nfev=60)
    cases = {"Literature parameters": parameters(),
             "After stage 1 (capacity, 1C)": parameters(am),
             "After stage 2 (transport, 5C)": parameters(am, stage2.x)}
    rows = [dict(parameters=k, rmse_1C_mV=round(rmse_mv(p, d1, 1), 1), rmse_5C_mV=round(rmse_mv(p, d5, 5), 1))
            for k, p in cases.items()]
    table = pd.DataFrame(rows)
    fitted = pd.DataFrame([dict(parameter=AM, multiplier=round(am, 4))] +
                          [dict(parameter=n, multiplier=round(10 ** v, 3)) for n, v in
                           zip(["Negative particle diffusivity", "Positive particle diffusivity",
                                "Exchange-current density (both electrodes)"], stage2.x)])
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    table.to_csv(os.path.join(HERE, "results", "step1_calibration_rmse.csv"), index=False)
    fitted.to_csv(os.path.join(HERE, "results", "step1_fitted_parameters.csv"), index=False)
    print(table.to_string(index=False)); print(fitted.to_string(index=False))

    fig, ax = ps.subplots(1, 2, figsize=(7.2, 2.8))
    for k, (d, c) in enumerate([(d1, 1), (d5, 5)]):
        ax[k].plot(d[:, 0] / 60, d[:, 1], "o", ms=3, mfc="none", color=ps.C["black"], label="Experiment (Ecker 2015)")
        for (name, p), col, ls in zip(cases.items(), [ps.C["blue"], ps.C["sky"], ps.C["orange"]], [":", "-.", "-"]):
            if name.startswith("After stage 1"):
                continue
            t, v = simulate(p, c)
            ax[k].plot(t / 60, v, color=col, ls=ls, lw=1.2, label="SPMe, " + name[0].lower() + name[1:])
        ax[k].set(xlabel="Time (min)", ylabel="Voltage (V)", ylim=(2.5, 4.25))
        ax[k].set_title(f"{c}C discharge", fontsize=8)
        ps.panel(ax[k], "AB"[k])
    h, l = ax[0].get_legend_handles_labels()
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    fig.legend(h, l, frameon=False, loc="upper center", ncol=3)
    fig.savefig(os.path.join(HERE, "figures", "fig1_calibration.png"))


if __name__ == "__main__":
    main()
