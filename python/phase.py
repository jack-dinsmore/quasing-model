import matplotlib.pyplot as plt
import numpy as np
import os
from util import load, get_stats_steep

WINDOW = (0.25, 0.55)
plt.style.use("root")

def post_process(load_result, eta):
    betas, mag, sus = load_result
    t2 = np.tan((eta + 1.) / 4. * np.pi)
    t1 = 1 / t2
    return (betas / t1, mag, sus)

def get_phase(start_name):
    etas = []
    temp_criticals = []
    for f in os.listdir("../data/output/"):
        if f.startswith(start_name):
            eta = float(f[len(start_name):-4])
            betas, mag, sus = post_process(load(f), eta)
            temp = get_stats_steep(betas, mag)
            etas.append(eta)
            temp_criticals.append(temp)
    etas = np.array(etas)
    temp_criticals = np.array(temp_criticals)
    argsort = np.argsort(etas)
    return etas[argsort], temp_criticals[argsort]

def plot(ax, name, label, color, marker):
    etas, transitions = get_phase(f"{name}-")

    ax.plot(etas, transitions, color=color, marker=marker, markerfacecolor=color, markeredgecolor='k', zorder=-1, label=label)
    ax.set_xlabel("Anisotropy $\eta$")
    ax.set_ylabel("Critical Temperature $T_c$")
    ax.set_xlim(-1,1)
    ax.set_ylim(0,None)

betas, mag, sus = load("ising-square.dat")
print("Square Ising", get_stats_steep(betas, mag))
betas, mag, sus = load("xy-square.dat")
print("Square XY", get_stats_steep(betas, mag))
betas, mag, sus = load("heisenberg-square.dat")
print("Square Heisenberg", get_stats_steep(betas, mag))
betas, mag, sus = load("ising-penrose.dat")
print("Penrose Ising", get_stats_steep(betas, mag))
betas, mag, sus = load("xy-penrose.dat")
print("Penrose XY", get_stats_steep(betas, mag))
betas, mag, sus = load("heisenberg-penrose.dat")
print("Penrose Heisenberg", get_stats_steep(betas, mag))


fig, ax = plt.subplots(figsize=(6, 4))
plot(ax, "rect-ising", "Ising", "gold", "o")
plot(ax, "rect-xy", "XY", "orange", "^")
plot(ax, "rect-heisenberg", "Heisenberg", "darkgoldenrod", "d")
fig.legend(ncol=3, loc="upper center", bbox_to_anchor=(0.5, 1.06))
fig.tight_layout()
fig.savefig(f"../figs/rect-phase.png")
fig.savefig(f"../figs/rect-phase.pdf")

fig, ax = plt.subplots(figsize=(6, 4))
plot(ax, "einstein-ising", "Ising", "gold", "o")
plot(ax, "einstein-xy", "XY", "orange", "^")
plot(ax, "einstein-heisenberg", "Heisenberg", "darkgoldenrod", "d")
fig.legend(ncol=3, loc="upper center", bbox_to_anchor=(0.5, 1.06))
fig.tight_layout()
fig.savefig(f"../figs/einstein-phase.png")
fig.savefig(f"../figs/einstein-phase.pdf")