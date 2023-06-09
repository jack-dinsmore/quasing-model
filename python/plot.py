import matplotlib.pyplot as plt
import numpy as np
from util import load, get_stats_steep
plt.style.use("root")

def display_crystal(ising_name, xy_name, save_name):
    fig, (ax_mag, ax_sus) = plt.subplots(figsize=(6, 7), nrows=2, sharex=True)
    lines = [
        (f"{ising_name}.dat", "crimson", "dodgerblue", "o"),
        (f"{xy_name}.dat", "lightsalmon", "cyan", "^"),
    ]
    t_crits = []
    for filename, red, blue, marker in lines:
        if filename[:4] == "None": continue
        betas, mag, sus = load(filename)
        t_crit = get_stats_steep(betas, mag)
        t_crits.append(t_crit)
        ax_sus.axvline(t_crit, linestyle="dotted", linewidth=1, color='k', zorder=-2)
        ax_mag.axvline(t_crit, linestyle="dotted", linewidth=1, color='k', zorder=-2)
        ax_mag.plot(1/betas, mag, color=blue, zorder=-1)
        ax_mag.scatter(1/betas, mag, facecolor=blue, marker=marker)
        ax_sus.plot(1/betas, sus, color=red, zorder=-1)
        ax_sus.scatter(1/betas, sus, facecolor=red, marker=marker)
    ax_sus.set_xlabel("Critical temperature $T_c$")
    ax_mag.set_ylabel("Magnetization $M$")
    ax_sus.set_ylabel("Susceptibility $\chi$")
    ax_sus.set_ylim(1e-1,None)
    # ax_mag.set_ylim(0,1)
    ax_sus.set_yscale('log')
    ax_sus.set_xlim(1/betas[-1], 1/betas[0])
    for t_crit, (_,_,_,marker) in zip(t_crits, lines):
        ax_mag.scatter([t_crit], ax_mag.get_ylim()[1], facecolor="lightgray", clip_on=False, zorder=10, marker=marker)
        ax_sus.scatter([t_crit], ax_sus.get_ylim()[0], facecolor="lightgray", clip_on=False, zorder=10, marker=marker)
    ax_sus.plot([],[],color='k', marker='o', markerfacecolor='lightgray',label="Ising model")
    ax_sus.plot([],[],color='k', marker='^', markerfacecolor='lightgray',label="XY model")
    fig.legend(ncol=2, loc="upper center", bbox_to_anchor=(0.5, 0.95))
    plt.subplots_adjust(hspace=0.05)
    fig.savefig(f"../figs/{save_name}.png", bbox_inches="tight")
    fig.savefig(f"../figs/{save_name}.pdf", bbox_inches="tight")

display_crystal("rect-heisenberg-1.00000000", None, "heisenberg")
# display_crystal("tim-square", None, "test")
display_crystal("ising-square", "xy-square", "square")
display_crystal("ising-penrose", "xy-penrose", "penrose")