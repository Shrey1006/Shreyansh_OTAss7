#!/usr/bin/env python3
"""OPTIONAL plotting helper (not part of the algorithms). Reads results/fronts/*.csv and results/true_*.csv
written by ./moo log and saves PNG figures to results/plots/.  Needs: python3, matplotlib, numpy."""
import glob, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

os.makedirs("results/plots", exist_ok=True)
def load(p):
    try:
        return np.loadtxt(p, delimiter=",", skiprows=1 if "fronts" in p and "P8" not in p else 0, ndmin=2)
    except Exception:
        return None
for f in sorted(glob.glob("results/fronts/alg*_*.csv")):
    name = os.path.basename(f)[:-4]; prob = name.split("_", 1)[1]
    d = load(f)
    if d is None or d.shape[1] < 2: continue
    plt.figure(figsize=(4.5, 3.6))
    t = "results/true_%s.csv" % prob
    if os.path.exists(t) and os.path.getsize(t) > 0 and np.loadtxt(t, delimiter=",", ndmin=2).shape[1] >= 2:
        T = np.loadtxt(t, delimiter=",", ndmin=2); plt.plot(T[:, 0], T[:, 1], "k-", lw=1, label="true front")
    plt.scatter(d[:, 0], d[:, 1], s=8, label=name)
    plt.xlabel("f1"); plt.ylabel("f2"); plt.legend(fontsize=7); plt.tight_layout()
    plt.savefig("results/plots/%s.png" % name, dpi=130); plt.close()
if os.path.exists("results/P13_trajectories.csv"):
    D = np.genfromtxt("results/P13_trajectories.csv", delimiter=",", names=True)
    plt.figure(figsize=(5, 3.6))
    for k in ["fixed_L1", "fixed_L2", "mgda_L1", "mgda_L2", "gradnorm_L1", "gradnorm_L2"]: plt.plot(D["iter"], D[k], label=k)
    plt.yscale("log"); plt.legend(fontsize=7); plt.xlabel("iteration"); plt.ylabel("task loss"); plt.tight_layout()
    plt.savefig("results/plots/P13_losses.png", dpi=130)
print("plots saved in results/plots/")
