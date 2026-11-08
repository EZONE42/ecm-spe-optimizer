"""Bruce–Vincent t+ from KIT Zenodo 7428858 files (Re/-Im only, no frequency column).
R_b  = Re(Z) at first local minimum of -Im (end of bulk arc)
R_tot= Re(Z) where the interfacial arc meets the real axis: circle fit (Kasa) on arc points,
       cross-checked with the low-frequency valley of -Im
R_i  = R_tot - R_b
"""
import sys, glob, numpy as np, pandas as pd

def load_eis(fn):
    d = pd.read_csv(fn, sep="\t", skiprows=1).dropna()
    return d.iloc[:, 0].to_numpy(float), d.iloc[:, 1].to_numpy(float)   # Re, -Im

def circle_fit(x, y):
    A = np.c_[x, y, np.ones_like(x)]
    b = x**2 + y**2
    c = np.linalg.lstsq(A, b, rcond=None)[0]
    xc, yc = c[0]/2, c[1]/2
    r = np.sqrt(c[2] + xc**2 + yc**2)
    h = np.sqrt(max(r**2 - yc**2, 0))
    return xc - h, xc + h, xc, yc, r

def analyse(fn):
    re, mim = load_eis(fn)
    i1 = next(i for i in range(1, len(mim)-1) if mim[i] <= mim[i-1] and mim[i] <= mim[i+1])
    Rb = re[i1]
    ip = i1 + int(np.argmax(mim[i1:]))                    # interfacial arc apex
    after = mim[ip:]
    iv = None
    for j in range(1, len(after)-1):
        if after[j] <= after[j-1] and after[j] <= after[j+1] and after[j] < 0.2*mim[ip]:
            iv = ip + j; break
    end = iv if iv is not None else len(re)-1
    sel = slice(i1 + 3, end)                              # arc points
    lo, hi, xc, yc, r = circle_fit(re[sel], mim[sel])
    Rtot_valley = re[iv] if iv is not None else np.nan
    return dict(file=fn.split("/")[-1], Rb=Rb, Rtot_circle=hi, Rtot_valley=Rtot_valley,
                Rtot_last=re[-1], arc_left=lo, depression_deg=np.degrees(np.arcsin(-yc/r)) if r else np.nan)

def load_ca(fn):
    d = pd.read_csv(fn, sep="\t").dropna()
    t_h, I_mA = d.iloc[:, 0].to_numpy(float), d.iloc[:, 1].to_numpy(float)
    return t_h*3600, I_mA*1e-3

def tplus(Rb0, Ri0, Ri_ss, I0, Iss, dV=0.05):
    return Iss*(dV - I0*Ri0) / (I0*(dV - Iss*Ri_ss))

if __name__ == "__main__":
    root = sys.argv[1]
    rows = [analyse(f) for f in sorted(glob.glob(root + "/*EIS_[0-9]_*.txt"))]
    E = pd.DataFrame(rows); print(E.round(2).to_string())
    pairs = {"Na": ("EIS_3", "EIS_4", "CA_2", 0.5026), "K": ("EIS_5", "EIS_6", "CA_3", 2.0106)}
    out = []
    for cat, (b, a, ca, Acm2) in pairs.items():
        B = E[E.file.str.contains(b)].iloc[0]; Aa = E[E.file.str.contains(a)].iloc[0]
        t, I = load_ca(glob.glob(root + f"/*{ca}_*.txt")[0])
        Iss = I[t >= t[-1]-600].mean()
        I0_meas = I[1]                                   # first point after step (row 0 = t≈0)
        for lab, kb, ka in [("circle", "Rtot_circle", "Rtot_circle"), ("valley/last", "Rtot_valley", "Rtot_valley")]:
            Rt0 = B[kb] if not np.isnan(B[kb]) else B["Rtot_last"]
            Rts = Aa[ka] if not np.isnan(Aa[ka]) else Aa["Rtot_last"]
            Ri0, Ris = Rt0 - B.Rb, Rts - Aa.Rb
            I0_calc = 0.05/(B.Rb + Ri0)
            out.append(dict(cation=cat, R_est=lab, Rb0=B.Rb, Ri0=Ri0, Rb_ss=Aa.Rb, Ri_ss=Ris,
                            I0_calc_uA=I0_calc*1e6, I0_meas_uA=I0_meas*1e6, Iss_uA=Iss*1e6,
                            Iss_over_I0=Iss/I0_calc,
                            t_plus_calcI0=tplus(B.Rb, Ri0, Ris, I0_calc, Iss),
                            t_plus_measI0=tplus(B.Rb, Ri0, Ris, I0_meas, Iss),
                            frac_dV_on_bulk=B.Rb/(B.Rb+Ri0),
                            ASR_int_per_side=Ri0*Acm2/2, CA_hours=t[-1]/3600))
    T = pd.DataFrame(out); pd.set_option("display.width", 250)
    print(T.round(4).to_string())
    import os
    out_dir = os.path.join(os.path.dirname(__file__), "..", "results")
    os.makedirs(out_dir, exist_ok=True)
    T.to_csv(os.path.join(out_dir, "bv_results_NaK.csv"), index=False)
    E.to_csv(os.path.join(out_dir, "bv_eis_intercepts.csv"), index=False)
