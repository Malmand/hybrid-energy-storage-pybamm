"""
Step 4: repeat the sizing of step 3 with measured grid frequency instead of the synthetic signal.

Data: Continental European frequency at 1 s for 2019, recorded by TransnetBW and pre-processed by
Kruse, Schäfer and Witthaut (Zenodo, doi:10.5281/zenodo.5105820). The TransnetBW data are not openly
licensed, so the repository holds only what tools/prepare_frequency.py derives from them: yearly
statistics, hourly statistics and three 2 h excerpts (the 2 h periods with median, 90th and 99th
percentile variability), which are used here.

Usage:  python step4_measured_frequency.py     (about 2 min)
"""
import json
import os

import numpy as np
import pandas as pd

import plotstyle as ps
import step3_ems_sizing as s3

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data", "frequency")
GRID = [(0, 1.0)] + [(s, t) for s in (5, 10, 20, 40) for t in (10, 30, 120)]


def sweep(p_req):
    rows = [s3.run(p_req, s, t)[1] for s, t in GRID]
    df = pd.DataFrame(rows)
    b = df.iloc[0]
    df["rms_reduction_pct"] = (1 - df.bat_rms_c_rate / b.bat_rms_c_rate) * 100
    df["reversal_reduction_pct"] = (1 - df.bat_reversals_per_h / b.bat_reversals_per_h) * 100
    return df


def signal_stats(f):
    d = f - 50
    x = d - d.mean()
    ac = lambda l: float(np.dot(x[:-l], x[l:]) / np.dot(x, x))
    return dict(std_mHz=d.std() * 1e3, step_std_mHz=np.diff(d).std() * 1e3, autocorr_60s=ac(60), autocorr_300s=ac(300))


def main():
    year = json.load(open(os.path.join(DATA, "freq_stats.json")))
    win = np.load(os.path.join(DATA, "freq_windows.npz"))
    signals = {"synthetic (step 3)": s3.frequency()}
    for k in sorted(win.files):
        q, day, hm = k.split("_")
        signals[f"measured, {q[1:]}th pct variability ({day[:4]}-{day[4:6]}-{day[6:]} {hm[:2]}:{hm[2:]})"] = \
            win[k].astype(float)
    rows, sweeps = [], {}
    for name, f in signals.items():
        df = sweep(s3.fcr_request(f))
        df.insert(0, "signal", name)
        sweeps[name] = df
        st = signal_stats(f)
        base = df.iloc[0]
        ok = df[(df.sc_strings > 0) & (df.fcr_error_rms_pct < 1) & (df.reversal_reduction_pct >= 80)]
        best = ok.sort_values(["sc_usable_kWh", "bat_rms_c_rate"]).iloc[0] if len(ok) else None
        rows.append(dict(signal=name, **{k: round(v, 3) for k, v in st.items()},
                         battery_only_reversals_per_h=base.bat_reversals_per_h,
                         battery_only_rms_c_rate=round(base.bat_rms_c_rate, 4),
                         smallest_bank_for_80pct=(f"{int(best.sc_strings)} strings, tau {best.tau_s:.0f} s"
                                                  if best is not None else "none")))
    pd.concat(sweeps.values()).round(4).to_csv(os.path.join(HERE, "results", "step4_sizing_measured.csv"), index=False)
    summ = pd.DataFrame(rows)
    summ.to_csv(os.path.join(HERE, "results", "step4_signal_comparison.csv"), index=False)
    pd.set_option("display.width", 250)
    print(summ.to_string(index=False))
    print(f"\nwhole year 2019: std {year['std_mHz']:.1f} mHz, 1 s step std {year['step_std_mHz']:.2f} mHz, "
          f"autocorrelation 60 s {year['autocorr']['60']:.2f}, 300 s {year['autocorr']['300']:.2f}, "
          f"FCR request reversals {year['fcr_reversals_per_h']:.0f} per hour, "
          f"{year['fcr_throughput_full_cycles_per_day']:.2f} full-power equivalent cycles per day")

    # figure: A, 10 min of synthetic vs measured; B, reversals vs usable SC energy for tau = 30 s
    fig, ax = ps.subplots(1, 2, figsize=(7.2, 2.9))
    names = list(signals)
    styles = [(ps.C["black"], "-", "o"), (ps.C["blue"], "-", "s"), (ps.C["orange"], "--", "^"), (ps.C["green"], "-.", "D")]
    t = np.arange(600) / 60
    for (c, ls, mk), n in zip(styles[:2], [names[0], names[2]]):
        ax[0].plot(t, (signals[n][3600:4200] - 50) * 1e3, color=c, ls=ls, lw=0.8,
                   label="Synthetic" if n == names[0] else "Measured (90th pct)")
    ax[0].set(xlabel="Time (min)", ylabel="Frequency deviation (mHz)")
    ax[0].legend(frameon=False, loc="lower left", bbox_to_anchor=(0.0, 1.02), ncol=2)
    ps.panel(ax[0], "A", y=1.1)
    for (c, ls, mk), n in zip(styles, names):
        g = sweeps[n]
        b = g.iloc[0].bat_reversals_per_h
        g = g[(g.tau_s == 30) & (g.sc_strings > 0)]
        lab = "Synthetic" if n == names[0] else "Measured, " + n.split(", ")[1].split(" (")[0].replace("th pct variability", "th pct")
        ax[1].plot(np.r_[0, g.sc_usable_kWh], np.r_[b, g.bat_reversals_per_h], marker=mk, color=c, ls=ls, lw=1, ms=4, label=lab)
    ax[1].set(xlabel="Usable supercapacitor energy (kWh)", ylabel="Battery current reversals per hour")
    ax[1].legend(frameon=False, loc="upper right", fontsize=6)
    ps.panel(ax[1], "B", y=1.1)
    ps.save(fig, os.path.join(HERE, "figures", "fig4_measured_frequency.png"))


if __name__ == "__main__":
    main()
