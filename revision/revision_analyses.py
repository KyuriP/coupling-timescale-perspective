"""
Revision analyses for "Rethinking Group Differences in Psychopathology Networks:
A Slow-Fast Perspective on Context and Symptom Activation".

Python port of revision_analyses.R (same estimator, same sections), written to
run where R is unavailable. Variable names and coding follow
R/01_prepare_analysis_data.R in the project repository.

Sections
  C0  Original SRL, median split: must reproduce manuscript Table 1 (check)
  A   Revised SRL: adversity coding of work variables, no healthcare cultural
      distance, winsorized activity/pack-years, minimum-completion rule
  B   Non-overlapping SRL (drops indicators that may reflect depression)
  C   Exact moderated Ising model, continuous SRL on activation (h) and
      interactions (J); LRTs, AIC/BIC, equal-slope test, implied probabilities
  D   Domain-specific moderators (separately and jointly)
  E   Symptom presence defined as item >= 2
  F   Ethnicity: moderation of J, SRL x ethnicity on h, models within groups
  G   Parameters for a redrawn Fig. 1

Model:  h_i(n) = Wh[n] . Bh[:, i],  J_ij(n) = WJ[n] . BJ[:, ij],
exact likelihood over all 2^9 = 512 symptom states, Newton-Raphson with the
exact Hessian (observed = expected information for this exponential family).

Usage:  python3 revision_analyses.py [--selftest] [--data data/helius.sav] [--out DIR]
"""
import argparse, itertools, json, os, sys, time, warnings
import numpy as np
warnings.filterwarnings("ignore", category=RuntimeWarning)   # nanmean of all-NaN rows
import pandas as pd
from scipy import stats, optimize
from scipy.special import logsumexp

SYMPTOMS = ["anh", "dep", "slp", "ene", "app", "glt", "con", "mot", "sui"]
PHQ_VARS = dict(anh="H1_WlbvRecent1", dep="H1_WlbvRecent2", slp="H1_WlbvRecent3",
                ene="H1_WlbvRecent4", app="H1_WlbvRecent5", glt="H1_WlbvRecent6",
                con="H1_WlbvRecent7", mot="H1_WlbvRecent_89", sui="H1_WlbvRecent10")
PHQSUM = "H1_PHQ9_sumscore"
MISSING = [-9, -1]

# ---------------------------------------------------------------------------
# 1. Exact (moderated) Ising estimator
# ---------------------------------------------------------------------------
class Space:
    def __init__(self, p):
        self.p = p
        S = np.array(list(itertools.product([0, 1], repeat=p)), dtype=float)
        self.pairs = np.array(list(itertools.combinations(range(p), 2)))
        self.m = len(self.pairs)
        S2 = S[:, self.pairs[:, 0]] * S[:, self.pairs[:, 1]]
        self.Tst = np.hstack([S, S2])
        sets = [(i,) for i in range(p)] + [tuple(r) for r in self.pairs]
        K = len(sets); keys = {}; cols = []
        self.uidx = np.zeros((K, K), dtype=int)
        for k in range(K):
            for l in range(k, K):
                u = tuple(sorted(set(sets[k]) | set(sets[l])))
                if u not in keys:
                    keys[u] = len(cols); cols.append(S[:, list(u)].prod(axis=1))
                self.uidx[k, l] = self.uidx[l, k] = keys[u]
        self.Tsub = np.column_stack(cols)


def pair_products(X, pairs):
    return X[:, pairs[:, 0]] * X[:, pairs[:, 1]]


class Map:
    """theta -> phi = c(vec(Bh), vec(BJ)) in R's column-major order."""
    def __init__(self, h_names, J_names, p, m, h_common=()):
        self.qh, self.qJ, self.p, self.m = len(h_names), len(J_names), p, m
        nphi = self.qh * p + self.qJ * m
        tid = np.zeros(nphi, dtype=int); labs = []; common = {}
        for i in range(p):
            for a, nm in enumerate(h_names):
                idx = a + i * self.qh
                if nm in h_common:
                    if nm not in common:
                        labs.append(f"h*:{nm}"); common[nm] = len(labs) - 1
                    tid[idx] = common[nm]
                else:
                    labs.append(f"h{i+1}:{nm}"); tid[idx] = len(labs) - 1
        for k in range(m):
            for a, nm in enumerate(J_names):
                idx = self.qh * p + a + k * self.qJ
                labs.append(f"J{k+1}:{nm}"); tid[idx] = len(labs) - 1
        self.M = np.zeros((nphi, len(labs))); self.M[np.arange(nphi), tid] = 1
        self.labels = labs

    def split(self, phi):
        Bh = phi[: self.qh * self.p].reshape(self.p, self.qh).T
        BJ = phi[self.qh * self.p:].reshape(self.m, self.qJ).T
        return Bh, BJ


def ising_eval(theta, X2, Wh, WJ, mp, sp, w, hessian=False):
    p, m = sp.p, sp.m; K = p + m
    Bh, BJ = mp.split(mp.M @ theta)
    eta = np.hstack([Wh @ Bh, WJ @ BJ])
    E = eta @ sp.Tst.T
    lz = logsumexp(E, axis=1)
    ll = float(np.sum(w * ((eta * X2).sum(1) - lz)))
    P = np.exp(E - lz[:, None]); del E
    mu = P @ sp.Tst
    R = (X2 - mu) * w[:, None]
    g_phi = np.concatenate([(Wh.T @ R[:, :p]).ravel(order="F"),
                            (WJ.T @ R[:, p:]).ravel(order="F")])
    out = dict(ll=ll, grad=mp.M.T @ g_phi)
    if hessian:
        EU = P @ sp.Tsub; del P
        D = lambda k: Wh if k < p else WJ
        off = lambda k: k * mp.qh if k < p else mp.qh * p + (k - p) * mp.qJ
        q = lambda k: mp.qh if k < p else mp.qJ
        H = np.zeros((mp.M.shape[0],) * 2)
        for k in range(K):
            Dk = D(k) * w[:, None]
            C = EU[:, sp.uidx[k, k:]] - mu[:, [k]] * mu[:, k:]
            for j, l in enumerate(range(k, K)):
                blk = -(Dk * C[:, [j]]).T @ D(l)
                ik = slice(off(k), off(k) + q(k)); il = slice(off(l), off(l) + q(l))
                H[ik, il] = blk
                if l != k: H[il, ik] = blk.T
        out["H"] = mp.M.T @ H @ mp.M
    return out


class Fit:
    pass


def fit_ising(X, Wh, WJ, h_names, J_names, h_common=(), w=None, label="",
              tol=1e-8, maxit=100):
    X = np.asarray(X, float); n, p = X.shape
    if np.isnan(X).any() or np.isnan(Wh).any() or np.isnan(WJ).any():
        raise ValueError(f"[{label}] missing values in symptoms or design matrices")
    sp = SPACES.setdefault(p, Space(p))
    X2 = np.hstack([X, pair_products(X, sp.pairs)])
    w = np.ones(n) if w is None else w
    mp = Map(h_names, J_names, p, sp.m, h_common)
    theta = np.zeros(len(mp.labels))
    pm = np.clip(X.mean(0), 1e-3, 1 - 1e-3)
    for i in range(p):
        lab = f"h{i+1}:(Intercept)"
        if lab in mp.labels: theta[mp.labels.index(lab)] = np.log(pm[i] / (1 - pm[i]))
    ev = ising_eval(theta, X2, Wh, WJ, mp, sp, w, True)
    for it in range(1, maxit + 1):
        try: step = np.linalg.solve(-ev["H"], ev["grad"])
        except np.linalg.LinAlgError: step = np.linalg.pinv(-ev["H"]) @ ev["grad"]
        s = 1.0
        while True:
            new = ising_eval(theta + s * step, X2, Wh, WJ, mp, sp, w)
            if np.isfinite(new["ll"]) and new["ll"] >= ev["ll"] - 1e-10: break
            s /= 2
            if s < 1e-8: break
        theta = theta + s * step
        dec = s * float(step @ ev["grad"])
        ev = ising_eval(theta, X2, Wh, WJ, mp, sp, w, True)
        if abs(dec) < tol and np.max(np.abs(ev["grad"])) < 1e-4: break
    try: V = np.linalg.inv(-ev["H"])
    except np.linalg.LinAlgError: V = np.linalg.pinv(-ev["H"])
    f = Fit()
    f.theta = pd.Series(theta, index=mp.labels); f.se = pd.Series(np.sqrt(np.clip(np.diag(V), 0, None)), index=mp.labels)
    f.ll = ev["ll"]; f.k = len(theta); f.n = w.sum()
    f.AIC = -2 * f.ll + 2 * f.k; f.BIC = -2 * f.ll + np.log(f.n) * f.k
    f.max_grad = float(np.max(np.abs(ev["grad"]))); f.iter = it; f.label = label
    f.p, f.m, f.pairs, f.map = p, sp.m, sp.pairs, mp
    f.h_names, f.J_names = list(h_names), list(J_names)
    return f


SPACES = {}


def design(df, terms):
    cols = [np.ones(len(df))] + [df[t].to_numpy(float) for t in terms]
    return np.column_stack(cols), ["(Intercept)"] + list(terms)


def fit_terms(X, df, h_terms, J_terms, label, h_common=()):
    Wh, hn = design(df, h_terms); WJ, jn = design(df, J_terms)
    return fit_ising(X, Wh, WJ, hn, jn, h_common, label=label)


def lrt(restricted, full, test=""):
    st = 2 * (full.ll - restricted.ll); df = full.k - restricted.k
    return dict(test=test, restricted=restricted.label, full=full.label,
                chisq=st, df=df, p=stats.chi2.sf(st, df))


def fit_table(*fits):
    return pd.DataFrame([dict(model=f.label, k=f.k, LL=f.ll, AIC=f.AIC, BIC=f.BIC,
                              max_grad=f.max_grad, iter=f.iter) for f in fits])


def h_coef(f, term):
    labs = [f"h{i+1}:{term}" for i in range(f.p)]
    if labs[0] not in f.theta.index: labs = [f"h*:{term}"] * f.p
    est = f.theta[labs].to_numpy(); se = f.se[labs].to_numpy()
    pv = 2 * stats.norm.sf(np.abs(est / se))
    return pd.DataFrame(dict(symptom=SYMPTOMS[: f.p], term=term, est=est, se=se,
                             lo=est - 1.96 * se, hi=est + 1.96 * se, OR=np.exp(est),
                             OR_lo=np.exp(est - 1.96 * se), OR_hi=np.exp(est + 1.96 * se),
                             p=pv, p_holm=holm(pv)))


def J_coef(f, term):
    labs = [f"J{k+1}:{term}" for k in range(f.m)]
    est = f.theta[labs].to_numpy(); se = f.se[labs].to_numpy()
    pv = 2 * stats.norm.sf(np.abs(est / se))
    edges = [f"{SYMPTOMS[a]}--{SYMPTOMS[b]}" for a, b in f.pairs]
    return pd.DataFrame(dict(edge=edges, term=term, est=est, se=se, lo=est - 1.96 * se,
                             hi=est + 1.96 * se, p=pv, p_holm=holm(pv)))


def holm(p):
    p = np.asarray(p); o = np.argsort(p); n = len(p)
    adj = np.maximum.accumulate((n - np.arange(n)) * p[o]).clip(max=1)
    out = np.empty(n); out[o] = adj; return out


def four_models(X, df, mod, cov_h=(), cov_J=(), tag=""):
    mod = [mod] if isinstance(mod, str) else list(mod)
    cov_h, cov_J = list(cov_h), list(cov_J)
    full = fit_terms(X, df, cov_h + mod, cov_J + mod, tag + "both moderated")
    hOnly = fit_terms(X, df, cov_h + mod, cov_J, tag + "activation moderated only")
    JOnly = fit_terms(X, df, cov_h, cov_J + mod, tag + "interactions moderated only")
    none = fit_terms(X, df, cov_h, cov_J, tag + "neither moderated")
    return dict(fits=dict(full=full, hOnly=hOnly, JOnly=JOnly, none=none),
                table=fit_table(full, hOnly, JOnly, none),
                lrt=pd.DataFrame([lrt(hOnly, full, "interaction moderation"),
                                  lrt(JOnly, full, "activation moderation")]))


def implied_marginals(f, wh, wj):
    sp = SPACES[f.p]; Bh, BJ = f.map.split(f.map.M @ f.theta.to_numpy())
    eta = np.concatenate([np.asarray(wh) @ Bh, np.asarray(wj) @ BJ])
    e = sp.Tst @ eta; pr = np.exp(e - e.max()); pr /= pr.sum()
    return pr @ sp.Tst[:, : f.p]


# ---------------------------------------------------------------------------
# 2. Self-test (same checks as the R version)
# ---------------------------------------------------------------------------
def selftest():
    rng = np.random.default_rng(20261020)
    p = 9; sp = SPACES.setdefault(p, Space(p)); m = sp.m; n = 6000
    z = rng.normal(size=n); g = rng.binomial(1, .5, n).astype(float)
    Bh = np.vstack([rng.normal(-1.5, .5, p), rng.uniform(.1, .5, p), rng.normal(0, .2, p)])
    BJ = np.vstack([rng.normal(.4, .4, m), np.r_[np.zeros(m - 3), .5, -.4, .3]])
    Wh = np.column_stack([np.ones(n), z, g]); WJ = np.column_stack([np.ones(n), z])
    eta = np.hstack([Wh @ Bh, WJ @ BJ]); E = eta @ sp.Tst.T
    P = np.exp(E - logsumexp(E, axis=1)[:, None])
    st = np.array([rng.choice(len(sp.Tst), p=r) for r in P])
    X = sp.Tst[st, :p]; df = pd.DataFrame(dict(z=z, g=g))
    hn, jn = ["(Intercept)", "z", "g"], ["(Intercept)", "z"]
    mp = Map(hn, jn, p, m); X2 = np.hstack([X, pair_products(X, sp.pairs)])
    th = rng.normal(0, .1, len(mp.labels)); w = np.ones(n)
    ev = ising_eval(th, X2, Wh, WJ, mp, sp, w, True)
    idx = rng.choice(len(th), 5, replace=False); eps = 1e-5
    fd = [(ising_eval(th + eps * np.eye(len(th))[j], X2, Wh, WJ, mp, sp, w)["ll"] -
           ising_eval(th - eps * np.eye(len(th))[j], X2, Wh, WJ, mp, sp, w)["ll"]) / (2 * eps) for j in idx]
    print("Self-test 1: max |grad - finite diff| =", np.max(np.abs(ev["grad"][idx] - fd)))
    e0 = np.eye(len(th))[idx[0]] * eps
    fdH = (ising_eval(th + e0, X2, Wh, WJ, mp, sp, w)["grad"] - ising_eval(th - e0, X2, Wh, WJ, mp, sp, w)["grad"]) / (2 * eps)
    print("              max |Hessian col - finite diff| =", np.max(np.abs(ev["H"][:, idx[0]] - fdH)))
    f = fit_terms(X, df, ["z", "g"], ["z"], "sim")
    hz = h_coef(f, "z")
    print(f"Self-test 2: h slope corr = {np.corrcoef(Bh[1], hz.est)[0,1]:.3f}, "
          f"coverage = {np.mean((Bh[1] >= hz.lo) & (Bh[1] <= hz.hi)):.2f}")
    fg = fit_terms(X, df, ["g"], ["g"], "dummy")
    f0 = fit_terms(X[g == 0], df[g == 0], [], [], "g0"); f1 = fit_terms(X[g == 1], df[g == 1], [], [], "g1")
    print(f"Self-test 3: LL dummy {fg.ll:.4f} vs g0+g1 {f0.ll + f1.ll:.4f}")
    fm = four_models(X, df, "z", cov_h=["g"], tag="sim: ")
    print(fm["table"].to_string()); print(fm["lrt"].to_string())


# ---------------------------------------------------------------------------
# 3. Data preparation
# ---------------------------------------------------------------------------
def zscore(x):
    return (x - np.nanmean(x)) / np.nanstd(x, ddof=1)


def winsor(x, lo=.01, hi=.99):
    a, b = np.nanquantile(x, [lo, hi]); return np.clip(x, a, b)


def row_mean_min(Z, min_obs):
    nobs = np.sum(~np.isnan(Z), axis=1)
    mu = np.nanmean(Z, axis=1)
    mu[nobs < min_obs] = np.nan
    return mu


def load(path):
    import pyreadstat
    cols = ["H1_geslacht", "H1_lft", "H1_InkHhMoeite", "H1_WerkSit", "H1_Opleid",
            "H1_Arbeidsparticipatie", "H1_BeroepsNiveau", "H1_PCS12", "H1_Squash_totmwk",
            "H1_PackYears", "H1_AlcoholConsumption", "H1_SBSQ_meanscore", "H1_SSQT",
            "H1_Discr_sumscore", "H1_StressWerk", "H1_StressHuish", "H1_MoeiteNLtaal",
            "H1_CultDistMeanScore6_corrected", "H1_etniciteit", PHQSUM] + list(PHQ_VARS.values())
    df, meta = pyreadstat.read_sav(path, usecols=cols)
    df = df.apply(pd.to_numeric, errors="coerce")
    df = df.mask(df.isin(MISSING))
    return df


def prepare(raw):
    r = raw.copy()
    ethn = r["H1_etniciteit"].where(r["H1_etniciteit"] != 6)
    dutch = (ethn == 1)
    sh = r["H1_StressHuish"].where(r["H1_StressHuish"].isin([1, 2, 3, 4]))
    sw = r["H1_StressWerk"].where(r["H1_StressWerk"].isin([1, 2, 3, 4]))
    cult = r["H1_CultDistMeanScore6_corrected"].where(lambda x: x >= 0)
    lang = r["H1_MoeiteNLtaal"]
    out = pd.DataFrame(index=r.index)
    out["age"] = r["H1_lft"]; out["female"] = (r["H1_geslacht"] == 2).astype(float)
    out.loc[r["H1_geslacht"].isna(), "female"] = np.nan
    out["ethn"] = ethn; out["nonDutch"] = np.where(ethn.isna(), np.nan, (ethn != 1).astype(float))
    out["PHQsum"] = r[PHQSUM]
    for s, v in PHQ_VARS.items(): out[s] = r[v]

    # ---- original SRL, exactly as R/01_prepare_analysis_data.R ----
    o = dict(inc=r["H1_InkHhMoeite"], work=r["H1_WerkSit"], edu=r["H1_Opleid"], arb=r["H1_Arbeidsparticipatie"],
             occ=r["H1_BeroepsNiveau"], pcs=r["H1_PCS12"], act=r["H1_Squash_totmwk"], smk=r["H1_PackYears"],
             alc=r["H1_AlcoholConsumption"], lit=r["H1_SBSQ_meanscore"], sup=r["H1_SSQT"],
             dis=r["H1_Discr_sumscore"], sw=sw, sh=sh,
             lang=lang.where(~(dutch & lang.isna()), 0), cult=cult.where(~(dutch & cult.isna()), 0))
    z = {k: zscore(v.to_numpy(float)) for k, v in o.items()}
    for k in ["edu", "occ", "pcs", "act", "lit", "sup"]: z[k] = -z[k]
    ses = np.nanmean(np.column_stack([z[k] for k in ["inc", "work", "arb", "edu", "occ", "lit", "sup"]]), 1)
    stv = np.nanmean(np.column_stack([z[k] for k in ["dis", "lang", "cult", "sw", "sh"]]), 1)
    hlt = np.nanmean(np.column_stack([z[k] for k in ["pcs", "act", "smk", "alc"]]), 1)
    out["srl_orig"] = np.nanmean(np.column_stack([ses, stv, hlt]), 1)

    # ---- revised SRL ----
    ws = r["H1_WerkSit"]; ap = r["H1_Arbeidsparticipatie"]
    rv = dict(
        ses=dict(inc=(r["H1_InkHhMoeite"], 1), work=(ws.isin([6, 7, 8]).astype(float).where(ws.notna()), 1),
                 edu=(r["H1_Opleid"], -1), arb=(ap.isin([3, 4]).astype(float).where(ap.notna()), 1),
                 occ=(r["H1_BeroepsNiveau"], -1), lit=(r["H1_SBSQ_meanscore"], -1), sup=(r["H1_SSQT"], -1)),
        str=dict(dis=(r["H1_Discr_sumscore"], 1), lang=(lang.where(~(dutch & lang.isna()), 0), 1),
                 sw=(sw, 1), sh=(sh, 1)),
        hlt=dict(pcs=(r["H1_PCS12"], -1), act=(pd.Series(winsor(r["H1_Squash_totmwk"].to_numpy(float))), -1),
                 smk=(pd.Series(winsor(r["H1_PackYears"].to_numpy(float))), 1), alc=(r["H1_AlcoholConsumption"], 1)))
    OVERLAP = {"work", "arb", "sup", "sw", "sh", "pcs", "act"}
    miss_rows = []
    def build(drop=()):
        D = {}; nind = 0
        for dom, inds in rv.items():
            Zd = np.column_stack([d * zscore(np.asarray(x, float)) for k, (x, d) in inds.items() if k not in drop])
            D[dom] = row_mean_min(Zd, int(np.ceil(0.5 * Zd.shape[1])))
            nind = nind + np.sum(~np.isnan(Zd), 1)
        Dm = np.column_stack(list(D.values()))
        return row_mean_min(Dm, 2), D, nind
    out["srl"], D, out["n_ind"] = build()
    for dom in D: out[dom] = D[dom]
    out["srl_nov"], _, _ = build(OVERLAP)
    for dom, inds in rv.items():
        for k, (x, d) in inds.items():
            miss_rows.append(dict(domain=dom, indicator=k, overlap=k in OVERLAP))
    return out, pd.DataFrame(miss_rows)


# ---------------------------------------------------------------------------
# 4. Main
# ---------------------------------------------------------------------------
def main(data, outdir, sections):
    want = lambda x: (sections == "all") or (x in sections.split(","))
    os.makedirs(outdir, exist_ok=True)
    log = open(os.path.join(outdir, f"results_{sections.replace(',', '-')}.txt"), "w")
    def say(*a):
        s = " ".join(str(x) for x in a); print(s); log.write(s + "\n"); log.flush()
    pd.set_option("display.width", 200); pd.set_option("display.max_columns", 30)
    t0 = time.time()
    say("Revision analyses", time.strftime("%Y-%m-%d %H:%M"))
    raw = load(data); say("Records in source dataset:", len(raw))
    d, inds = prepare(raw)
    items = d[SYMPTOMS].to_numpy(float)
    complete = ~np.isnan(items).any(1)
    res = {}

    # ---- C0: original SRL, replicating Table 1 ----
    s_desc = d["srl_orig"].notna() & d["PHQsum"].notna()
    s_net = s_desc & complete
    med_o = d.loc[s_desc, "srl_orig"].median()
    if want("C0"):
        say(f"\n== C0: original SRL. descriptive n = {s_desc.sum()}, network n = {s_net.sum()} "
            f"(manuscript: 23,689 / 23,283)")
        g_o = (d.loc[s_net, "srl_orig"] > med_o).astype(float).to_numpy()
        say(f"   Low/High n = {(g_o == 0).sum()}/{(g_o == 1).sum()} (manuscript: 11,717/11,566)")
        X1o = (items[s_net.to_numpy()] >= 1).astype(float)
        mg0 = four_models(X1o, pd.DataFrame(dict(g=g_o)), "g", tag="orig split: ")
        say(mg0["table"].to_string()); say(mg0["lrt"].to_string())
        say("   Table 1 LL: -93555.20 / -93614.78 / -93611.51 / -95229.03")
        mg0["table"].to_csv(os.path.join(outdir, "C0_original_table1_check.csv"), index=False)

    # ---- A: revised SRL and analytic sample ----
    flow = [("Source records", len(d))]
    k1 = d["srl"].notna() & d["PHQsum"].notna(); flow.append(("Revised SRL computable + PHQ total", int(k1.sum())))
    k2 = k1 & d[["age", "female", "nonDutch"]].notna().all(1); flow.append(("+ age, sex, ethnicity observed", int(k2.sum())))
    k3 = k2 & complete; flow.append(("+ all nine PHQ items (network sample)", int(k3.sum())))
    flow = pd.DataFrame(flow, columns=["step", "n"])
    net = d[k3].reset_index(drop=True); Xo = items[k3.to_numpy()]
    X1 = (Xo >= 1).astype(float); X2c = (Xo >= 2).astype(float)
    net["age_z"] = zscore(net["age"].to_numpy(float))
    for v in ["srl", "srl_nov", "ses", "str", "hlt", "srl_orig"]: net[v + "_z"] = zscore(net[v].to_numpy(float))
    if want("A"):
        say("\n== A: exclusion flow (revised SRL) =="); say(flow.to_string(index=False))
        flow.to_csv(os.path.join(outdir, "A_exclusion_flow.csv"), index=False)
        say("   original-SRL sample dropped by the completion rule:", int((s_desc & ~k1).sum()))
        say(f"   indicators per person: min {net.n_ind.min():.0f}, median {net.n_ind.median():.0f}, max {net.n_ind.max():.0f}")
        say("   corr(revised, original) = %.3f; corr(revised, non-overlap) = %.3f" %
            (net.srl.corr(net.srl_orig), net.srl.corr(net.srl_nov)))
        for v in ["srl", "srl_nov"]:
            ct = pd.crosstab(net[v] > d.loc[k1, v].median(), net.srl_orig > med_o)
            say(f"   median-split agreement {v} vs original: {np.trace(ct.values) / ct.values.sum():.3f}")
    COV_H = ["age_z", "female", "nonDutch"]
    if os.environ.get("DEBUG_EXPORT"):
        pd.concat([pd.DataFrame(X1, columns=SYMPTOMS), net[COV_H + ["srl_z"]]], axis=1).to_csv(os.environ["DEBUG_EXPORT"], index=False)

    def run_mod(X, mod, tag, file):
        ok = net[mod].notna().to_numpy()
        if not ok.all(): say(f"   {mod}: using {ok.sum()} of {len(ok)} participants with this score observed")
        dfm, X = net[ok].reset_index(drop=True), X[ok]
        fm = four_models(X, dfm, mod, cov_h=COV_H, cov_J=["nonDutch"], tag=tag)
        cs = fit_terms(X, dfm, COV_H + [mod], ["nonDutch", mod], tag + "common slope", h_common=[mod])
        fm["lrt"] = pd.concat([fm["lrt"], pd.DataFrame([lrt(cs, fm["fits"]["full"], "equal SRL slopes on h")])], ignore_index=True)
        fm["h"] = h_coef(fm["fits"]["full"], mod); fm["J"] = J_coef(fm["fits"]["full"], mod)
        say(fm["table"].to_string()); say(fm["lrt"][["test", "chisq", "df", "p"]].to_string())
        say("SRL effect on activation (log-odds per SD; full model):"); say(fm["h"].round(4).to_string())
        say("Edges with Holm p < .05 for moderation:"); say(fm["J"][fm["J"].p_holm < .05].round(4).to_string())
        for k in ["table", "lrt", "h", "J"]: fm[k].to_csv(os.path.join(outdir, f"{file}_{k}.csv"), index=False)
        say(f"   [{time.time() - t0:.0f}s]")
        return fm

    gr = (net["srl"] > d.loc[k1, "srl"].median()).astype(float).to_numpy()
    if want("C0b"):
        say("\n== C0b: revised SRL, median split ==")
        mgr = four_models(X1, pd.DataFrame(dict(g=gr)), "g", tag="revised split: ")
        say(mgr["table"].to_string()); say(mgr["lrt"].to_string())
        mgr["table"].to_csv(os.path.join(outdir, "C0b_revised_split.csv"), index=False)

    full_cache = {}
    def main_full():
        if "f" not in full_cache:
            full_cache["f"] = fit_terms(X1, net, COV_H + ["srl_z"], ["nonDutch", "srl_z"], "C: both moderated")
        return full_cache["f"]

    if want("C"):
        say("\n== C: moderated Ising, revised SRL (continuous), item >= 1 ==")
        modC = run_mod(X1, "srl_z", "C: ", "C_srl"); full_cache["f"] = modC["fits"]["full"]
        prof = lambda zv: ([1, 0, 0, 1, zv], [1, 1, zv])   # mean age, male, non-Dutch
        imp = pd.DataFrame({lab: implied_marginals(main_full(), *prof(zv)) for lab, zv in [("SRL -1 SD", -1), ("SRL +1 SD", 1)]}, index=SYMPTOMS)
        say("Model-implied symptom probabilities (mean age, male, non-Dutch):"); say(imp.round(3).to_string())
        imp.to_csv(os.path.join(outdir, "C_implied_marginals.csv"))

    if want("Corig"):
        say("\n== C-orig: moderated Ising with the ORIGINAL SRL (continuous) ==")
        run_mod(X1, "srl_orig_z", "C-orig: ", "C_srl_original")

    if want("B"):
        say("\n== B: non-overlapping SRL ==")
        run_mod(X1, "srl_nov_z", "B: ", "B_srl_nonoverlap")

    for v in ["ses_z", "str_z", "hlt_z"]:
        if want("D" + v[:3]):
            say(f"\n== D: domain {v} alone =="); run_mod(X1, v, f"D[{v}]: ", f"D_{v}")
    DOMS = ["ses_z", "str_z", "hlt_z"]
    if want("Djoint"):
        say("\n== D: domains jointly ==")
        ok = net[DOMS].notna().all(1).to_numpy(); dn, Xd = net[ok].reset_index(drop=True), X1[ok]
        say(f"   participants with all three domain scores: {ok.sum()}")
        dj = fit_terms(Xd, dn, COV_H + DOMS, ["nonDutch"] + DOMS, "D joint")
        rows = [lrt(fit_terms(Xd, dn, COV_H + [x for x in DOMS if x != v], ["nonDutch"] + DOMS, f"drop {v} from h"), dj,
                    f"activation moderation by {v} given others") for v in DOMS]
        rows += [lrt(fit_terms(Xd, dn, COV_H + DOMS, ["nonDutch"] + [x for x in DOMS if x != v], f"drop {v} from J"), dj,
                     f"interaction moderation by {v} given others") for v in DOMS]
        djl = pd.DataFrame(rows); say(djl[["test", "chisq", "df", "p"]].to_string())
        djh = pd.concat([h_coef(dj, v) for v in DOMS]); say(djh.round(4).to_string())
        djl.to_csv(os.path.join(outdir, "D_joint_lrt.csv"), index=False); djh.to_csv(os.path.join(outdir, "D_joint_h.csv"), index=False)

    if want("E"):
        say("\n== E: symptom presence = item >= 2 ==")
        say("prevalence:", dict(zip(SYMPTOMS, X2c.mean(0).round(3))))
        run_mod(X2c, "srl_z", "E: ", "E_cut2")
    if want("Esplit"):
        say("\n== E: item >= 2, revised-SRL median split ==")
        mge = four_models(X2c, pd.DataFrame(dict(g=gr)), "g", tag="E split: ")
        say(mge["table"].to_string()); say(mge["lrt"][["test", "chisq", "df", "p"]].to_string())

    if want("F"):
        say("\n== F: ethnicity ==")
        full = main_full()
        f1 = fit_terms(X1, net, COV_H + ["srl_z"], ["srl_z"], "no ethnicity on J")
        net["srl_x_nonDutch"] = net["srl_z"] * net["nonDutch"]
        f2 = fit_terms(X1, net, COV_H + ["srl_z", "srl_x_nonDutch"], ["nonDutch", "srl_z"], "SRL x ethnicity on h")
        fl = pd.DataFrame([lrt(f1, full, "ethnicity moderation of interactions"), lrt(full, f2, "SRL x ethnicity on activation")])
        say(fl[["test", "chisq", "df", "p"]].to_string()); say(h_coef(f2, "srl_x_nonDutch").round(4).to_string())
        fl.to_csv(os.path.join(outdir, "F_lrt.csv"), index=False)
    for e, nm in [(0, "Dutch"), (1, "nonDutch")]:
        if want("F" + nm):
            s = (net["nonDutch"] == e).to_numpy()
            say(f"\n== F: within {nm} (n = {s.sum()}) ==")
            sub = net[s].reset_index(drop=True)
            fm = four_models(X1[s], sub, "srl_z", cov_h=["age_z", "female"], tag=f"F[{nm}]: ")
            say(fm["table"].to_string()); say(fm["lrt"][["test", "chisq", "df", "p"]].to_string())
            hc = h_coef(fm["fits"]["full"], "srl_z"); say(hc.round(4).to_string())
            hc.to_csv(os.path.join(outdir, f"F_within_{nm}_h.csv"), index=False)
            mg = four_models(X1[s], pd.DataFrame(dict(g=gr[s])), "g", tag=f"F[{nm}] split: ")
            say(mg["table"].to_string()); say(mg["lrt"][["test", "chisq", "df", "p"]].to_string())

    if want("G"):
        say("\n== G: Fig. 1 parameters ==")
        full = main_full()
        sp = SPACES[9]; Bh, BJ = full.map.split(full.map.M @ full.theta.to_numpy())
        h0 = Bh[0]; Jv = BJ[0].copy()
        def dist(h, Jvec):
            e = sp.Tst @ np.concatenate([h, Jvec]); pr = np.exp(e - e.max()); return pr / pr.sum()
        def phi(pr):
            S = sp.Tst[:, :9]; mu = pr @ S; C = (S * pr[:, None]).T @ S - np.outer(mu, mu)
            return C / np.sqrt(np.outer(np.diag(C), np.diag(C)))
        gc = [k for k, (a, b) in enumerate(sp.pairs) if (SYMPTOMS[a], SYMPTOMS[b]) == ("glt", "con")][0]
        JB = Jv.copy(); JB[gc] += 1.0
        target = phi(dist(h0, JB)); iu = np.triu_indices(9, 1)
        obj = lambda h: np.sum((phi(dist(h, Jv)) - target)[iu] ** 2)
        hC = optimize.minimize(obj, h0, method="BFGS").x
        say("RMS phi difference scenario B vs C (group 2): %.4f" % np.sqrt(np.mean((phi(dist(hC, Jv)) - target)[iu] ** 2)))
        pd.DataFrame(dict(symptom=SYMPTOMS, h_group1=h0, h_scenC_group2=hC)).to_csv(os.path.join(outdir, "G_fig1_activation.csv"), index=False)
        for nm, M_ in [("group1", phi(dist(h0, Jv))), ("scenB_group2", target), ("scenC_group2", phi(dist(hC, Jv)))]:
            pd.DataFrame(M_, index=SYMPTOMS, columns=SYMPTOMS).to_csv(os.path.join(outdir, f"G_phi_{nm}.csv"))
    say(f"\nDone in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--data", default="data/helius.sav")
    ap.add_argument("--out", default="results/revision_analyses")
    ap.add_argument("--sections", default="all", help="comma list: C0,A,C0b,C,Corig,B,Dses,Dstr,Dhlt,Djoint,E,Esplit,F,FDutch,FnonDutch,G")
    a = ap.parse_args()
    if a.selftest: selftest()
    else: main(a.data, a.out, a.sections)
