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
    
def get_stats_steep(betas, mag):
    max_slope = None
    max_slope_temp = None
    max_slope_mag= None
    temps = 1/betas
    if len(mag) <= 1: return np.nan, 1
    for i in range(len(mag) - 1):
        dtemp = temps[i+1] - temps[i]
        if dtemp == 0: continue
        dmag = mag[i+1] - mag[i]
        if max_slope is None or dmag/dtemp < max_slope:
            max_slope = dmag/dtemp
            max_slope_mag = (mag[i] + mag[i+1]) / 2
            max_slope_temp = (temps[i] + temps[i+1]) / 2
    intercept = max_slope_temp + (0.5 - max_slope_mag) / max_slope
    return intercept
