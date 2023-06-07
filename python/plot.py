import matplotlib.pyplot as plt
import numpy as np

def load(filename):
    with open(f"../data/output/{filename}", 'r') as f:
        betas = f.readline().split(',')[:-1]
        betas = np.array([float(b) for b in betas])
        mag = f.readline().split(',')[:-1]
        mag = np.array([float(s) for s in mag])
        sus = f.readline().split(',')[:-1]
        sus = np.array([float(s) for s in sus])
        argsort = np.argsort(betas)
        betas = betas[argsort]
        mag = mag[argsort]
        sus = sus[argsort]
        return betas, mag, sus
    
betas, mag, sus = load("ising-square.dat")

fig, (ax_mag, ax_sus) = plt.subplots(figsize=(6, 10), nrows=2, sharex=True)

ax_mag.scatter(betas, mag, color="C0")
ax_mag.plot(betas, mag, color="C0")
ax_sus.scatter(betas, sus, color="C1")
ax_sus.plot(betas, sus, color="C1")
ax_sus.set_xlabel("beta")
ax_mag.set_ylabel("Magnetization")
ax_sus.set_ylabel("Susceptibility")

# ax_mag.set_xlim(0.35, 0.55)

fig.tight_layout()
fig.savefig("../figs/params.png")