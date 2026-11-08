"""KIT BPE blocking-electrode EIS (Zenodo 7428858) -> R_b, sigma(T), Arrhenius/VTF, heating vs cooling.
File format (checked 2026-10-29): 'cycle number', 'Re(Z)/Ohm', '-Im(Z)/Ohm' ; 54 cycles x 114 points.
Each cycle = 57-point sweep out + 57-point sweep back (Nyquist traced twice) -> both halves analysed.
R_b: Re(Z) at the -Im minimum between bulk arc and blocking spike (intercept method);
     cross-check: Kasa circle fit on the bulk arc (points before the minimum) -> right real-axis crossing.
"""
import sys, glob, re as _re, numpy as np, pandas as pd
from scipy.optimize import curve_fit

L_UM = {1: 95, 2: 140, 3: 75, 4: 70, 5: 145, 6: 85, 7: 105, 8: 140, 9: 90}
A_CM2 = 0.785
KB = 8.617333e-5

def cyc_to_T(c):
    return (5 + 10*((c-1)//3), "heating", (c-1) % 3 + 1) if c <= 27 else (85 - 10*((c-28)//3), "cooling", (c-28) % 3 + 1)

def circle_right(x, y):
    A = np.c_[x, y, np.ones_like(x)]
    c = np.linalg.lstsq(A, x**2 + y**2, rcond=None)[0]
    xc, yc = c[0]/2, c[1]/2
    r2 = c[2] + xc**2 + yc**2
    return xc + np.sqrt(max(r2 - yc**2, 0)) if r2 > yc**2 else np.nan

def rb_half(re, mim):
    # valley = minimum phase angle atan(-Im/Re): bulk-arc end / blocking-spike start.
    # (global min of -Im fails at low T, where the high-frequency start of the huge bulk arc has small -Im)
    ph = np.arctan2(mim, re)
    i = int(np.argmin(ph))
    rb_int = re[i]
    rb_circ = circle_right(re[:i+1], mim[:i+1]) if i >= 5 else np.nan
    return rb_int, rb_circ, np.degrees(ph[i])

def analyse(fn):
    n = int(_re.search(r"EIS_(\d)_", fn).group(1))
    d = pd.read_csv(fn, sep="\t").iloc[:, :3]; d.columns = ["c", "re", "mim"]
    rows = []
    for c, g in d.groupby("c"):
        c = int(c); re_, mim = g.re.to_numpy(), g.mim.to_numpy()
        h = len(re_)//2
        a = rb_half(re_[:h], mim[:h]); b = rb_half(re_[h:][::-1], mim[h:][::-1])
        T, dirn, rep = cyc_to_T(c)
        rb = np.mean([a[0], b[0]])
        rows.append(dict(sample=n, file=fn.split("/")[-1].split("-", 1)[1], cycle=c, T_C=T, direction=dirn, rep=rep,
                         Rb_int_out=a[0], Rb_int_back=b[0], Rb_circle=np.nanmean([a[1], b[1]]),
                         valley_phase_deg=a[2], Rb=rb, L_um=L_UM[n],
                         sigma=L_UM[n]*1e-4/(rb*A_CM2)))
    return pd.DataFrame(rows)

def vtf(T, A, B, T0): return A*T**-0.5*np.exp(-B/(T-T0))

if __name__ == "__main__":
    import os
    out_dir = os.path.join(os.path.dirname(__file__), "..", "results")
    os.makedirs(out_dir, exist_ok=True)
    R = pd.concat([analyse(f) for f in sorted(glob.glob(sys.argv[1] + "/*EIS_[1-9]_BPE*.txt"))])
    R["cation"] = R.file.str.extract(r"BPE-(\w+)TFSI")[0]; R["EO"] = R.file.str.extract(r"_(\d+)-1")[0].astype(int)
    R.to_csv(os.path.join(out_dir, "eis_anchor_results.csv"), index=False)
    G = R.groupby(["sample", "cation", "EO", "T_C", "direction"]).agg(
        sigma=("sigma", "mean"), sigma_sd=("sigma", "std"), Rb=("Rb", "mean"),
        circ_vs_int=("Rb_circle", lambda s: np.nan)).reset_index()
    G["circ_vs_int"] = R.groupby(["sample", "T_C", "direction"]).apply(
        lambda g: np.nanmean(g.Rb_circle)/g.Rb.mean()).values
    G.to_csv(os.path.join(out_dir, "eis_anchor_sigma_by_T.csv"), index=False)
    fits = []
    for (s, cat, eo, dirn), g in G.groupby(["sample", "cation", "EO", "direction"]):
        T = g.T_C.to_numpy() + 273.15; sg = g.sigma.to_numpy()
        hi = T >= 328.15  # Arrhenius on 55-85 C (amorphous regime)
        p = np.polyfit(1/T[hi], np.log(sg[hi]*T[hi]), 1)
        try:
            pv, _ = curve_fit(vtf, T, sg, p0=[1, 1000, 200], bounds=([1e-6, 50, 100], [1e6, 5000, T.min()-5]), maxfev=20000)
            r2 = 1 - np.sum((np.log10(vtf(T, *pv)) - np.log10(sg))**2)/np.sum((np.log10(sg)-np.log10(sg).mean())**2)
        except Exception:
            pv, r2 = [np.nan]*3, np.nan
        fits.append(dict(sample=s, cation=cat, EO=eo, direction=dirn, Ea_55_85_eV=-p[0]*KB,
                         VTF_A=pv[0], VTF_B_K=pv[1], VTF_T0_K=pv[2], VTF_R2_log=r2))
    F = pd.DataFrame(fits); F.to_csv(os.path.join(out_dir, "eis_anchor_temperature_fits.csv"), index=False)
    pd.set_option("display.width", 250)
    P = G.pivot_table(index=["cation", "EO", "T_C"], columns="direction", values="sigma")
    P["cool/heat"] = P.cooling/P.heating
    print(P.map(lambda v: f"{v:.3g}").to_string())
    print(F.round(4).to_string())
    print("rep CV max %:", (G.sigma_sd/G.sigma*100).max().round(2), " circle/intercept range:", G.circ_vs_int.min().round(3), G.circ_vs_int.max().round(3))
