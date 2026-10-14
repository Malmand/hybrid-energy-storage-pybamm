"""
Step 3: energy management and sizing of the battery-supercapacitor HESS, and FMU check.

Application: a stationary 1 MW / 1 MWh Li-ion battery delivering frequency containment reserve
(FCR). FCR makes the battery follow every fluctuation of grid frequency, so it sees many small,
fast current reversals. A supercapacitor bank can absorb the fast part.

  - Frequency: synthetic 2 h signal at 1 s, a slow and a fast Ornstein-Uhlenbeck component
    (standard deviations 20 and 6 mHz, time constants 60 s and 3 s).
  - FCR request: 1 MW at +/-200 mHz, linear in between (no deadband).
  - Energy management: the battery follows a low-pass filtered request (time constant tau);
    the supercapacitor supplies the rest, with a feedback term that pulls its state of charge
    back to 50%.
  - Sizing: supercapacitor strings (150 cells of 3000 F in series) x filter time constant.
  - Metrics: battery RMS C-rate, battery current reversals per hour (micro-cycles), FCR
    tracking error, losses.

The best design is then re-run through the FMU (fmu/HessFmu.fmu) with a co-simulation loop in
fmpy, as Simulink would run it, and compared with the direct Python run.

Usage:  python step3_ems_sizing.py        (about 1 min)
"""
import os

import numpy as np
import pandas as pd

import plotstyle as ps
from hess_plant import HessPlant

HERE = os.path.dirname(os.path.abspath(__file__))
DT, HOURS, P_FCR = 1.0, 2, 1.0e6


def frequency(seed=3):
    rng = np.random.default_rng(seed)
    n = int(HOURS * 3600 / DT)
    f = np.zeros(n)
    slow = fast = 0.0
    for k in range(1, n):
        slow += -slow / 60 * DT + 0.020 * np.sqrt(2 / 60 * DT) * rng.standard_normal()
        fast += -fast / 3 * DT + 0.006 * np.sqrt(2 / 3 * DT) * rng.standard_normal()
        f[k] = 50 + slow + fast
    return f


def fcr_request(f):
    return -P_FCR * np.clip((f - 50) / 0.2, -1, 1)          # positive = discharge


class FilterEMS:
    def __init__(self, tau, k_soc=2.0e5, has_sc=True):
        self.tau, self.k, self.has_sc, self.p_lp = tau, k_soc, has_sc, 0.0

    def split(self, p_req, soc_sc):
        if not self.has_sc:
            return p_req, 0.0
        self.p_lp += (p_req - self.p_lp) * DT / self.tau
        p_bat = self.p_lp - self.k * (soc_sc - 0.5)       # low SC charge -> battery supplies more
        return p_bat, p_req - p_bat


def run(p_req, sc_parallel, tau):
    plant = HessPlant(bat_kWh=1000, bat_series=200, sc_series=150, sc_parallel=sc_parallel)
    ems = FilterEMS(tau, has_sc=sc_parallel > 0)
    n = len(p_req)
    out = {k: np.zeros(n) for k in ("p_bat", "p_sc", "i_bat", "soc_sc", "loss")}
    soc_sc = plant.sc.soc if sc_parallel else 0.5
    for k in range(n):
        pb_set, ps_set = ems.split(p_req[k], soc_sc)
        r = plant.step(pb_set, ps_set, DT)
        # anything the supercapacitor cannot deliver falls back to the battery in the next step
        out["p_bat"][k], out["p_sc"][k], out["i_bat"][k] = r["p_bat"], r["p_sc"], r["i_bat"]
        out["soc_sc"][k], out["loss"][k] = r["soc_sc"], r["loss_bat"] + r["loss_sc"]
        soc_sc = r["soc_sc"] if sc_parallel else 0.5
    q = plant.bat.Q
    i = out["i_bat"]
    sign = np.sign(np.where(np.abs(i) > 0.01 * q, i, 0))
    nz = sign[sign != 0]
    reversals = int(np.sum(nz[1:] != nz[:-1])) / HOURS
    delivered = out["p_bat"] + out["p_sc"]
    return out, dict(sc_strings=sc_parallel, sc_usable_kWh=round(plant.sc.energy_usable_kWh, 2), tau_s=tau,
                     bat_rms_c_rate=float(np.sqrt(np.mean(i ** 2)) / q),
                     bat_reversals_per_h=reversals,
                     fcr_error_rms_pct=float(np.sqrt(np.mean((delivered - p_req) ** 2)) / P_FCR * 100),
                     losses_kWh=float(out["loss"].sum() * DT / 3.6e6),
                     sc_soc_min=float(np.nanmin(out["soc_sc"])) if sc_parallel else np.nan)


def fmu_check(p_req, sc_parallel, tau):
    from fmpy import extract, read_model_description
    from fmpy.fmi2 import FMU2Slave
    path = os.path.join(HERE, "fmu", "HessFmu.fmu")
    md = read_model_description(path)
    vr = {v.name: v.valueReference for v in md.modelVariables}
    unzip = extract(path)
    fmu = FMU2Slave(guid=md.guid, unzipDirectory=unzip, modelIdentifier=md.coSimulation.modelIdentifier,
                    instanceName="hess")
    fmu.instantiate()
    fmu.setupExperiment(startTime=0.0)
    fmu.enterInitializationMode()
    fmu.setReal([vr["sc_parallel"], vr["bat_kWh"], vr["bat_series"], vr["sc_series"], vr["soc0"]],
                [float(sc_parallel), 1000.0, 200.0, 150.0, 0.5])
    fmu.exitInitializationMode()
    ems, losses = FilterEMS(tau), 0.0
    soc_sc = (0.8 ** 2 - 0.5 ** 2) / (1 - 0.5 ** 2)          # initial state of charge of the bank
    for k, p in enumerate(p_req):
        pb, psc = ems.split(p, soc_sc)
        fmu.setReal([vr["P_bat_set"], vr["P_sc_set"]], [pb, psc])
        fmu.doStep(currentCommunicationPoint=k * DT, communicationStepSize=DT)
        soc_bat, soc_sc, loss = fmu.getReal([vr["SOC_bat"], vr["SOC_sc"], vr["P_loss"]])
        losses += loss * DT / 3.6e6
    fmu.terminate()
    fmu.freeInstance()
    return soc_bat, soc_sc, losses


def main():
    f = frequency()
    p_req = fcr_request(f)
    rows, keep = [], {}
    _, base = run(p_req, 0, 1.0)
    rows.append(base)
    for strings in (5, 10, 20, 40):
        for tau in (10, 30, 120):
            out, m = run(p_req, strings, tau)
            rows.append(m)
            keep[(strings, tau)] = out
    df = pd.DataFrame(rows)
    b = df.iloc[0]
    df["rms_reduction_pct"] = (1 - df.bat_rms_c_rate / b.bat_rms_c_rate) * 100
    df["reversal_reduction_pct"] = (1 - df.bat_reversals_per_h / b.bat_reversals_per_h) * 100
    df = df.round(4)
    df.to_csv(os.path.join(HERE, "results", "step3_sizing.csv"), index=False)
    pd.set_option("display.width", 200)
    print(df.to_string(index=False))

    # pick the smallest bank that keeps FCR error below 1% and cuts reversals by at least 80%
    ok = df[(df.sc_strings > 0) & (df.fcr_error_rms_pct < 1) & (df.reversal_reduction_pct >= 80)]
    best = ok.sort_values(["sc_usable_kWh", "bat_rms_c_rate"]).iloc[0] if len(ok) else df.iloc[-1]
    s, tau = int(best.sc_strings), float(best.tau_s)
    out = keep[(s, tau)]
    soc_b, soc_s, loss_fmu = fmu_check(p_req, s, tau)
    direct_loss = out["loss"].sum() * DT / 3.6e6
    check = pd.DataFrame([dict(design=f"{s} strings, tau {tau:.0f} s", direct_losses_kWh=direct_loss,
                               fmu_losses_kWh=loss_fmu, direct_soc_sc_end=out["soc_sc"][-1], fmu_soc_sc_end=soc_s)])
    check.to_csv(os.path.join(HERE, "results", "step3_fmu_check.csv"), index=False)
    print(check.to_string(index=False))

    # figure
    fig, ax = ps.subplots(1, 2, figsize=(7.2, 2.9))
    w = slice(3600, 3600 + 600)
    t = np.arange(len(p_req))[w] / 60
    ax[0].plot(t, p_req[w] / 1e3, color=ps.C["black"], lw=0.8, label="FCR request")
    ax[0].plot(t, out["p_bat"][w] / 1e3, color=ps.C["blue"], lw=1.2, label="Battery")
    ax[0].plot(t, out["p_sc"][w] / 1e3, color=ps.C["orange"], lw=0.8, ls="--", label="Supercapacitor")
    ax[0].set(xlabel="Time (min)", ylabel="Power (kW)")
    ax[0].legend(frameon=False, loc="lower left", bbox_to_anchor=(0.0, 1.02), ncol=3)
    ps.panel(ax[0], "A", y=1.1)
    marks = {10: ("o", ps.C["blue"], "-"), 30: ("s", ps.C["orange"], "--"), 120: ("^", ps.C["green"], "-.")}
    for tau_, (mk, c, ls) in marks.items():
        g = df[(df.tau_s == tau_) & (df.sc_strings > 0)]
        ax[1].plot(g.sc_usable_kWh, g.reversal_reduction_pct, marker=mk, color=c, ls=ls, lw=1, ms=4,
                   label=f"tau = {tau_} s")
    ax[1].set(xlabel="Usable supercapacitor energy (kWh)", ylabel="Fewer battery current reversals (%)")
    ax[1].legend(frameon=False, loc="lower right")
    ps.panel(ax[1], "B", y=1.1)
    ps.save(fig, os.path.join(HERE, "figures", "fig3_ems_sizing.png"))


if __name__ == "__main__":
    main()
