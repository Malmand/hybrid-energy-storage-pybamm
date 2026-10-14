"""Condense one year of measured Continental European grid frequency (1 s) for the FCR study.

Source: Kruse, Schäfer, Witthaut, Pre-Processed Power Grid Frequency Time Series, Zenodo,
doi:10.5281/zenodo.5105820 (2019_cleansed.zip, TransnetBW recordings). The TransnetBW data are not
openly licensed, so only statistics and short excerpts used for the analysis are produced here.

Usage: python prepare_frequency.py 2019_cleansed.zip
Writes: freq_stats.json, freq_hourly.csv (per-hour statistics), freq_windows.npz (2 h excerpts)
"""
import json
import sys
import zipfile
import io
import numpy as np
import pandas as pd

ZIP = sys.argv[1]
P_FCR, FULL = 1.0, 0.2            # MW, full activation at 200 mHz


def load(zf, name):
    with zf.open(name) as fh:
        raw = fh.read()
    if name.endswith(".zip"):
        inner = zipfile.ZipFile(io.BytesIO(raw))
        n = [m for m in inner.namelist() if m.lower().endswith(".csv")][0]
        raw = inner.read(n)
    df = pd.read_csv(io.BytesIO(raw), header=None, names=["t", "f"], skiprows=1)
    df["t"] = pd.to_datetime(df.t, errors="coerce")
    df["f"] = pd.to_numeric(df.f, errors="coerce")
    return df.dropna(subset=["t"]).set_index("t").sort_index()


def main():
    zf = zipfile.ZipFile(ZIP)
    names = zf.namelist()
    print("archive contents:", names[:50])
    tb = [n for n in names if "transnet" in n.lower() and (n.lower().endswith(".csv") or n.lower().endswith(".zip"))]
    print("TransnetBW files:", tb)
    s = pd.concat([load(zf, n) for n in tb]).sort_index()
    s = s[~s.index.duplicated()]
    f = s.f
    print(f"{len(f)} samples, {f.isna().mean():.3%} NaN, {f.index.min()} to {f.index.max()}")
    df = (f - 50.0)
    req = -np.clip(df / FULL, -1, 1)                                 # p.u. of FCR power
    h = pd.DataFrame({"std_mHz": df.groupby(df.index.floor("h")).std() * 1e3,
                      "mean_mHz": df.groupby(df.index.floor("h")).mean() * 1e3,
                      "absmax_mHz": df.abs().groupby(df.index.floor("h")).max() * 1e3,
                      "throughput_pu_h": req.abs().groupby(req.index.floor("h")).sum() / 3600,
                      "valid": f.notna().groupby(f.index.floor("h")).mean()})
    h.to_csv("freq_hourly.csv", float_format="%.5g")
    sign = np.sign(req.where(req.abs() > 0.01, 0)).replace(0, np.nan).dropna()
    rev = (sign.diff().abs() > 0).sum() / (len(f.dropna()) / 3600)
    d1 = df.diff().dropna()
    # autocorrelation of the deviation at a few lags (1 s data)
    x = df.dropna().values[:5_000_000]
    x = x - x.mean()
    ac = {str(l): float(np.dot(x[:-l], x[l:]) / np.dot(x, x)) for l in (1, 3, 10, 30, 60, 300, 900)}
    stats = dict(samples=int(len(f)), nan_frac=float(f.isna().mean()),
                 std_mHz=float(df.std() * 1e3), p99_abs_mHz=float(df.abs().quantile(0.99) * 1e3),
                 max_abs_mHz=float(df.abs().max() * 1e3), frac_beyond_200mHz=float((df.abs() > 0.2).mean()),
                 frac_beyond_100mHz=float((df.abs() > 0.1).mean()),
                 step_std_mHz=float(d1.std() * 1e3), fcr_reversals_per_h=float(rev),
                 fcr_throughput_full_cycles_per_day=float(req.abs().sum() / 3600 / (len(f.dropna()) / 86400) / 2),
                 autocorr=ac,
                 minute_of_hour_mean_abs_mHz={str(m): float(v) for m, v in
                                              (df.abs().groupby(df.index.minute).mean() * 1e3).items()})
    json.dump(stats, open("freq_stats.json", "w"), indent=1)
    print(json.dumps({k: v for k, v in stats.items() if k != "minute_of_hour_mean_abs_mHz"}, indent=1))
    # 2 h excerpts: median, 90th and 99th percentile variability, complete data only
    hh = h[h.valid > 0.999].copy()
    roll = hh.std_mHz.rolling(2).mean().dropna()
    pick = {}
    for q in (0.5, 0.9, 0.99):
        t_end = (roll - roll.quantile(q)).abs().idxmin()
        t0 = t_end - pd.Timedelta(hours=1)
        w = f.loc[t0:t0 + pd.Timedelta(seconds=7199)]
        if len(w) == 7200 and w.notna().all():
            pick[f"q{int(q * 100)}_{t0:%Y%m%d_%H%M}"] = w.values.astype(np.float32)
    np.savez_compressed("freq_windows.npz", **pick)
    print("windows:", list(pick))


if __name__ == "__main__":
    main()
