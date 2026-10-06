"""Independent Bayesian fit of Model 1: same data, same structural and residual model, estimated by
Hamiltonian Monte Carlo (NumPyro, NUTS) with the random effects integrated by sampling rather than by the
Laplace approximation.

Structures: "primary" (clearance correlation rho and volume correlation r_V, as Model 1) and "full" (all 6
correlations, LKJ prior). Weakly informative priors: log typical clearance ~ N(log 3 L/h, 1), log typical volume
~ N(log 20 L, 1), between-patient SD ~ half-normal(1), residual SD ~ half-normal(0.5), correlations uniform
on (-1, 1) (LKJ concentration 1 for the full matrix). 4 chains x 1,500 warm-up + 2,000 draws, seed 20261005.
Outputs: outputs/model1_bayes_summary.csv
"""
import csv
import os
import sys
import time

import numpy as np

os.environ.setdefault("XLA_FLAGS", "--xla_force_host_platform_device_count=4")
import jax  # noqa: E402
import jax.numpy as jnp  # noqa: E402
import numpyro  # noqa: E402
import numpyro.distributions as dist  # noqa: E402
from numpyro.diagnostics import summary as np_summary  # noqa: E402
from numpyro.infer import MCMC, NUTS  # noqa: E402

jax.config.update("jax_enable_x64", True)
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import joint_popk_nlme as M  # noqa: E402

OUT = M.OUT
SEED = 20261005
TAU = M.TAU
DOSE = (M.DOSE["caz"], M.DOSE["avi"])


def arrays():
    subs = M.load()
    nmax = max(len(s.times["caz"]) for s in subs)
    n = len(subs)
    t = np.zeros((n, nmax)); mask = np.zeros((n, nmax), bool)
    y = np.zeros((2, n, nmax)); tinf = np.zeros(n)
    for i, s in enumerate(subs):
        k = len(s.times["caz"])
        assert np.array_equal(s.times["caz"], s.times["avi"])
        t[i, :k] = s.times["caz"]; mask[i, :k] = True; tinf[i] = s.t_inf
        y[0, i, :k] = s.logconc["caz"]; y[1, i, :k] = s.logconc["avi"]
    return jnp.array(t), jnp.array(mask), jnp.array(y), jnp.array(tinf)


def css(t, cl, v, dose, tinf):
    k = cl / v
    rate = dose / tinf
    c0 = (rate / cl) * (1 - jnp.exp(-k * tinf)) * jnp.exp(-k * (TAU - tinf)) / (1 - jnp.exp(-k * TAU))
    t_in = jnp.minimum(t, tinf)
    during = (rate / cl) * (1 - jnp.exp(-k * t_in)) + c0 * jnp.exp(-k * t_in)
    c_end = (rate / cl) * (1 - jnp.exp(-k * tinf)) + c0 * jnp.exp(-k * tinf)
    return jnp.where(t <= tinf, during, c_end * jnp.exp(-k * jnp.maximum(t - tinf, 0.0)))


def model(t, mask, y, tinf, structure="primary"):
    n = t.shape[0]
    log_cl = numpyro.sample("log_cl", dist.Normal(jnp.log(3.0), 1.0).expand([2]).to_event(1))
    log_v = numpyro.sample("log_v", dist.Normal(jnp.log(20.0), 1.0).expand([2]).to_event(1))
    om = numpyro.sample("omega", dist.HalfNormal(1.0).expand([4]).to_event(1))
    sig = numpyro.sample("sigma", dist.HalfNormal(0.5).expand([2]).to_event(1))
    c = numpyro.sample("c", dist.Uniform(-1.0, 1.0))
    u = numpyro.sample("u", dist.Normal(0.0, 1.0).expand([n, 4]).to_event(2))
    if structure == "primary":
        rho = numpyro.sample("rho", dist.Uniform(-1.0, 1.0))
        r_v = numpyro.sample("r_v", dist.Uniform(-1.0, 1.0))
        z0 = u[:, 0]; z1 = rho * u[:, 0] + jnp.sqrt(1 - rho ** 2) * u[:, 1]
        z2 = u[:, 2]; z3 = r_v * u[:, 2] + jnp.sqrt(1 - r_v ** 2) * u[:, 3]
        z = jnp.stack([z0, z1, z2, z3], axis=1)
    else:
        L = numpyro.sample("L_omega", dist.LKJCholesky(4, concentration=1.0))
        z = u @ L.T
        R = L @ L.T
        numpyro.deterministic("rho", R[0, 1])
        numpyro.deterministic("r_v", R[2, 3])
        numpyro.deterministic("r_cl_v_caz", R[0, 2])
        numpyro.deterministic("r_cl_v_avi", R[1, 3])
    cl = jnp.exp(log_cl)[None, :] * jnp.exp(om[None, :2] * z[:, :2])
    v = jnp.exp(log_v)[None, :] * jnp.exp(om[None, 2:] * z[:, 2:])
    p0 = jnp.log(jnp.clip(css(t, cl[:, :1], v[:, :1], DOSE[0], tinf[:, None]), 1e-10))
    p1 = jnp.log(jnp.clip(css(t, cl[:, 1:], v[:, 1:], DOSE[1], tinf[:, None]), 1e-10))
    e0 = (y[0] - p0) / sig[0]
    e1 = (y[1] - p1) / sig[1]
    q = (e0 ** 2 - 2 * c * e0 * e1 + e1 ** 2) / (1 - c ** 2)
    ll = -0.5 * q - jnp.log(2 * jnp.pi * sig[0] * sig[1] * jnp.sqrt(1 - c ** 2))
    numpyro.factor("loglik", jnp.sum(jnp.where(mask, ll, 0.0)))
    numpyro.deterministic("cl_typ", jnp.exp(log_cl))


def run(structure, key):
    t, mask, y, tinf = arrays()
    mcmc = MCMC(NUTS(model, target_accept_prob=0.95, max_tree_depth=10), num_warmup=1500, num_samples=2000,
                num_chains=4, chain_method="parallel", progress_bar=False)
    t0 = time.time()
    mcmc.run(key, t, mask, y, tinf, structure=structure, extra_fields=("diverging",))
    s = mcmc.get_samples(group_by_chain=True)
    div = int(np.sum(mcmc.get_extra_fields()["diverging"]))
    names = ["rho", "r_v", "c", "omega", "sigma", "cl_typ"] + (["r_cl_v_caz", "r_cl_v_avi"] if structure == "full" else [])
    summ = np_summary({k: np.asarray(s[k]) for k in names}, prob=0.95)
    rows = []
    for k in names:
        arr = np.asarray(s[k]).reshape(8000, -1)
        for j in range(arr.shape[1]):
            lab = k if arr.shape[1] == 1 else f"{k}[{j}]"
            d = arr[:, j]
            st = summ[k]
            rhat = float(np.ravel(st["r_hat"])[j]); ess = float(np.ravel(st["n_eff"])[j])
            rows.append(dict(structure=structure, parameter=lab, median=float(np.median(d)), mean=float(d.mean()),
                             q2_5=float(np.percentile(d, 2.5)), q97_5=float(np.percentile(d, 97.5)), r_hat=rhat,
                             ess_bulk=ess, divergences=div, seconds=time.time() - t0))
    rho = np.asarray(s["rho"]).ravel()
    rows.append(dict(structure=structure, parameter="P(rho >= 0.94)", median=float(np.mean(rho >= 0.94)), mean="",
                     q2_5="", q97_5="", r_hat="", ess_bulk="", divergences=div, seconds=time.time() - t0))
    rows.append(dict(structure=structure, parameter="P(rho > 0)", median=float(np.mean(rho > 0)), mean="",
                     q2_5="", q97_5="", r_hat="", ess_bulk="", divergences=div, seconds=time.time() - t0))
    return rows


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    print("jax", jax.__version__, "numpyro", numpyro.__version__, "devices", jax.local_device_count(), flush=True)
    rows = []
    for i, structure in enumerate(("primary", "full")):
        r = run(structure, jax.random.PRNGKey(SEED + i))
        rows += r
        for x in r:
            if x["parameter"] in ("rho", "r_v", "c", "P(rho >= 0.94)", "r_cl_v_caz", "r_cl_v_avi") or x["parameter"].startswith(("omega", "sigma")):
                print(f"  {structure:7} {x['parameter']:16} median {x['median']:.4f}  95% {x['q2_5'] if x['q2_5'] == '' else round(x['q2_5'], 4)} "
                      f"to {x['q97_5'] if x['q97_5'] == '' else round(x['q97_5'], 4)}  r_hat {x['r_hat']}  ess {x['ess_bulk']}  div {x['divergences']}", flush=True)
    path = os.path.join(OUT, "model1_bayes_summary.csv")
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.6f}" if isinstance(v, float) else v) for k, v in r.items()})
    print("wrote", os.path.relpath(path, os.path.dirname(OUT)))


if __name__ == "__main__":
    main()
