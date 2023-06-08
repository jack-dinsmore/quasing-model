import matplotlib.pyplot as plt
import numpy as np
from util import load

# plt.rcParams["text.usetex"] = True
plt.rcParams["font.size"] = 12
plt.rcParams["font.family"] = "Courier"
    
fig, (ax_mag, ax_sus) = plt.subplots(figsize=(6, 7), nrows=2, sharex=True)
lines = [
    ("ising-square.dat", "crimson", "dodgerblue", "o"),
    ("xy-square.dat", "lightsalmon", "cyan", "^"),
]
for filename, red, blue, marker in lines:
    betas, mag, sus = load(filename)
    ax_mag.plot(1/betas, mag, color=blue, zorder=-10)
    ax_mag.scatter(1/betas, mag, facecolor=blue, edgecolor='k', marker=marker)
    ax_sus.plot(1/betas, sus, color=red, zorder=-10)
    ax_sus.scatter(1/betas, sus, facecolor=red, edgecolor='k', marker=marker)
ax_sus.set_xlabel("Critical temperature $T_c$")
ax_mag.set_ylabel("Magnetization $M$")
ax_sus.set_ylabel("Susceptibility $\chi$")
ax_sus.set_ylim(1e-1,None)
ax_mag.set_ylim(0,1)
ax_sus.set_yscale('log')
ax_sus.set_xlim(1/betas[-1], 1/betas[0])
ax_sus.plot([],[],color='k', marker='o', markerfacecolor='lightgray', markeredgecolor='k',label="Ising model")
ax_sus.plot([],[],color='k', marker='^', markerfacecolor='lightgray', markeredgecolor='k',label="XY model")
fig.legend(ncol=2, loc="upper center")
# fig.tight_layout()
fig.savefig("../figs/square.png", bbox_inches="tight")




fig, (ax_mag, ax_sus) = plt.subplots(figsize=(6, 7), nrows=2, sharex=True)
lines = [
    ("ising-penrose.dat", "crimson", "dodgerblue", "o"),
    ("xy-penrose.dat", "lightsalmon", "cyan", "^"),
]
for filename, red, blue, marker in lines:
    betas, mag, sus = load(filename)
    ax_mag.plot(1/betas, mag, color=blue, zorder=-10)
    ax_mag.scatter(1/betas, mag, facecolor=blue, edgecolor='k', marker=marker)
    ax_sus.plot(1/betas, sus, color=red, zorder=-10)
    ax_sus.scatter(1/betas, sus, facecolor=red, edgecolor='k', marker=marker)
ax_sus.set_xlabel("Critical temperature $T_c$")
ax_mag.set_ylabel("Magnetization $M$")
ax_sus.set_ylabel("Susceptibility $\chi$")
ax_mag.set_ylim(0,1)
ax_sus.set_ylim(1e-1,None)
ax_sus.set_yscale('log')
ax_sus.set_xlim(1/betas[-1], 1/betas[0])
ax_sus.plot([],[],color='k', marker='o', markerfacecolor='lightgray', markeredgecolor='k',label="Ising model")
ax_sus.plot([],[],color='k', marker='^', markerfacecolor='lightgray', markeredgecolor='k',label="XY model")
fig.legend(ncol=2, loc="upper center")
fig.savefig("../figs/penrose.png", bbox_inches="tight")