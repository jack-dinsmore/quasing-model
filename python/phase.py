import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.optimize import minimize

WINDOW = (0.25, 0.55)

def load(filename):
    with open(f"../data/output/{filename}", 'r') as f:
        betas = f.readline().split(',')[:-1]
        betas = np.array([float(b) for b in betas])
        sus = f.readline().split(',')[:-1]
        sus = np.array([float(s) for s in sus])
        argsort = np.argsort(betas)
        betas = betas[argsort]
        sus = sus[argsort]
        mid_betas = []
        sus_deriv = []
        for i in range(len(sus) - 1):
            dbeta = betas[i+1] - betas[i]
            if dbeta == 0: continue
            dsus = sus[i+1] - sus[i]
            mid_betas.append((betas[i] + betas[i+1]) / 2)
            sus_deriv.append(dsus/dbeta)

        mask = (WINDOW[0] < betas) & (betas < WINDOW[1])
        mask = sus < 0.25
        sus = sus[mask]
        betas = betas[mask]
        return betas, sus, mid_betas, sus_deriv
    
def model_fn(param, exponent, scale, xs):
    tau = param / xs - 1
    output = (np.abs(tau)**exponent) * scale
    output[tau > 0] = 0
    return output

def chisq(params, xs, ys):
    param, exponent, scale = params
    return np.sum((ys - model_fn(param, exponent, scale, xs))**2 / 2)

def get_stats_fit(betas, sus):
    x0 = 0.5, 0.5, 0.5
    result = minimize(chisq, x0, args=(betas, sus), bounds=(
        (0.25,0.75),
        (0,3),
        (0,1)
    ))

    fig, ax = plt.subplots()
    ax.plot(betas, sus)
    ax.plot(betas, model_fn(result.x[0], result.x[1], result.x[2],betas))
    ax.plot(betas, model_fn(0.28, 0.5, 1,betas))
    fig.savefig("test.png")
    return result.x[0], result.x[1]

def get_stats_steep(betas, sus):
    max_slope = None
    max_slope_beta = None
    max_slope_sus = None
    if len(sus) <= 1: return np.nan, 1
    for i in range(len(sus) - 1):
        dbeta = betas[i+1] - betas[i]
        if dbeta == 0: continue
        dsus = sus[i+1] - sus[i]
        if max_slope is None or dsus/dbeta > max_slope:
            max_slope = dsus/dbeta
            max_slope_sus = (sus[i] + sus[i+1]) / 2
            max_slope_beta = (betas[i] + betas[i+1]) / 2
    intercept = max_slope_beta - max_slope_sus / max_slope
    return intercept, 1


def get_phase_attributes(start_name):
    t2s = []
    transitions = []
    exponents = []
    for f in os.listdir("../data/output/"):
        if f.startswith(start_name):
            t2 = float(f[len(start_name):-4])
            betas, sus, mid_betas, sus_deriv = load(f)
            if t2 < 1:
                fig,ax = plt.subplots()
                ax.plot(betas, sus)
                fig.savefig("test.png")
            transition, exponent = get_stats_steep(betas, sus)
            t2s.append(t2)
            transitions.append(transition)
            exponents.append(exponent)
    t2s = np.array(t2s)
    exponents = np.array(exponents)
    transitions = np.array(transitions)
    argsort = np.argsort(t2s)
    return t2s[argsort], transitions[argsort], exponents[argsort]

def plot(name):
    t2s, transitions, exponents = get_phase_attributes(f"{name}-")

    fig, (ax_transitions, ax_exponents) = plt.subplots(figsize=(6, 8), nrows=2, sharex=True)

    ax_transitions.scatter(t2s, 1/transitions, color="C0")
    ax_transitions.plot(t2s, 1/transitions, color="C0")
    ax_exponents.scatter(t2s, exponents, color="C1")
    ax_exponents.set_xlabel("t2")
    ax_exponents.set_ylabel("Critical Exponent")
    ax_transitions.set_ylabel("Tc")

    ax_transitions.set_xlim(0.1,2)
    ax_transitions.set_ylim(0,1)
    fig.tight_layout()
    fig.savefig(f"../figs/{name}-phase.png")

# plot("rect-ising")
# plot("einstein-ising")
# plot("rect-xy")
plot("einstein-xy")