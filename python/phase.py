import matplotlib.pyplot as plt
import numpy as np
import os
from util import load, get_stats_steep

WINDOW = (0.25, 0.55)
plt.style.use("root")

def middle(a,b,c):
    # Returns the center of the quadratic interpolation between three points
    num = a[0]**2 * (b[1] - c[1]) + b[0]**2 * (c[1] - a[1]) + c[0]**2 * (a[1] - b[1])
    denom = a[1] * (b[0] - c[0]) + b[1] * (c[0] - a[0]) + c[1] * (a[0] - b[0])
    return -0.5 * num / denom

def max_middle(xs, ys):
    s = np.array([xs,ys]).transpose()
    argmax = np.argmax(ys)
    return middle(s[argmax-1], s[argmax], s[argmax+1])
    
def max_left_right(xs, ys):
    res = 0
    for mask in [xs > 0, xs < 0]:
        argmax = np.argmax(ys[mask])
        s = np.array([xs,ys]).transpose()[mask]
        res += np.abs(middle(s[argmax-1], s[argmax], s[argmax+1]))
    return res / 2

def get_phase(start_name):
    etas = []
    temp_criticals = []
    for f in os.listdir("../data/output/"):
        if f.startswith(start_name):
            eta = float(f[len(start_name):-4])
            betas, mag, sus = load(f)
            temp = max(0,get_stats_steep(betas, mag))
            etas.append(eta)
            temp_criticals.append(temp)
    etas = np.array(etas)
    temp_criticals = np.array(temp_criticals)
    argsort = np.argsort(etas)
    return etas[argsort], temp_criticals[argsort]

def plot(ax, name, label, color, marker):
    etas, transitions = get_phase(f"{name}-")

    if name in ["rect-ising", "rect-xy"]:
        print(name, max_left_right(etas, transitions))

    if name in ["einstein-ising", "einstein-xy", "einstein-heisenberg"]:
        print(name, max_middle(etas, transitions))

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
betas, mag, sus = load("einstein-ising--0.32999998.dat")
print("Einstein Ising", get_stats_steep(betas, mag))
betas, mag, sus = load("einstein-xy--0.32999998.dat")
print("Einstein XY", get_stats_steep(betas, mag))
betas, mag, sus = load("einstein-heisenberg--0.32999998.dat")
print("Einstein Heisenberg", get_stats_steep(betas, mag))
print()



fig, ax = plt.subplots()
plot(ax, "rect-ising", "Ising", "gold", "o")
plot(ax, "rect-xy", "XY", "orange", "^")
plot(ax, "rect-heisenberg", "Heisenberg", "darkgoldenrod", "d")
fig.legend(ncol=3, loc="upper center", bbox_to_anchor=(0.5, 1.06))
fig.tight_layout()
fig.savefig(f"../figs/rect-phase.png")
fig.savefig(f"../figs/rect-phase.pdf")

fig, ax = plt.subplots()
plot(ax, "einstein-ising", "Ising", "gold", "o")
plot(ax, "einstein-xy", "XY", "orange", "^")
plot(ax, "einstein-heisenberg", "Heisenberg", "darkgoldenrod", "d")
fig.legend(ncol=3, loc="upper center", bbox_to_anchor=(0.5, 1.06))
fig.tight_layout()
fig.savefig(f"../figs/einstein-phase.png")
fig.savefig(f"../figs/einstein-phase.pdf")