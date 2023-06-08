import matplotlib.pyplot as plt
import numpy as np
import os
from util import load, get_stats_steep

WINDOW = (0.25, 0.55)


def get_phase(start_name):
    t2s = []
    beta_criticals = []
    for f in os.listdir("../data/output/"):
        if f.startswith(start_name):
            t2 = float(f[len(start_name):-4])
            betas, mag, sus = load(f)
            beta = get_stats_steep(betas, mag)
            t2s.append(t2)
            beta_criticals.append(beta)
    t2s = np.array(t2s)
    beta_criticals = np.array(beta_criticals)
    argsort = np.argsort(t2s)
    return t2s[argsort], beta_criticals[argsort]

def plot(name):
    t2s, transitions = get_phase(f"{name}-")

    fig, ax_transitions = plt.subplots(figsize=(6, 4))

    ax_transitions.scatter(t2s, 1/transitions, color="C0")
    ax_transitions.plot(t2s, 1/transitions, color="C0")
    ax_transitions.set_xlabel("t2")
    ax_transitions.set_ylabel("Tc")

    # ax_transitions.set_xlim(0.1,2)
    # ax_transitions.set_ylim(0,1)
    fig.tight_layout()
    fig.savefig(f"../figs/{name}-phase.png")

betas, mag, sus = load("ising-square.dat")
print("Square Ising", get_stats_steep(betas, mag))
betas, mag, sus = load("xy-square.dat")
print("Square XY", get_stats_steep(betas, mag))
betas, mag, sus = load("ising-penrose.dat")
print("Penrose Ising", get_stats_steep(betas, mag))
betas, mag, sus = load("xy-penrose.dat")
print("Penrose XY", get_stats_steep(betas, mag))

plot("rect-ising")
plot("einstein-ising")
plot("rect-xy")
plot("einstein-xy")