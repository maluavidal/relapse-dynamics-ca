import numpy as np
import time
from rules import evolve, S, E, R, C, T, P

# --- CONFIGURATION ---
# --- OPTIMIZED CONFIGURATION ---
SHARED = {
    "N": 100,
    "steps": 500,
    "runs": 100,
    "implement_at": 150
}

# EXPERIMENT A: Ideal (Success)
# - High support (5%), Good recovery (alpha=0.05), Low burnout (eta=0.1)
EXP_A = {
    "params": {
        'beta': 0.25, 'gamma': 0.25, 'sigma': 0.60,
        'alpha': 0.05, 'eta': 0.10, 'delta': 1e-4,
        'tau': 0.10, 'theta': 0.25, 'rho': 0.05
    },
    "p_supp": 0.05, "p_init_E": 0.02
}

# EXPERIMENT B: Chronic (Failure)
# - Low support (1%), Poor recovery (alpha=0.015), High burnout (eta=0.5)
# - CHANGED: Lowered delta (0.005) to stop flickering, raised beta (0.35) to maintain pressure
EXP_B = {
    "params": {
        'beta': 0.35, 'gamma': 0.25, 'sigma': 0.10,
        'alpha': 0.015, 'eta': 0.50, 'delta': 0.005,
        'tau': 0.05, 'theta': 0.15, 'rho': 0.10
    },
    "p_supp": 0.01, "p_init_E": 0.02
}

# EXPERIMENT C: Intervention
# Starts with B params. At t=150, we simulate a "Public Health Initiative":
# 1. More Mentors (p_supp 1% -> 6%)
# 2. Better Prevention (sigma 0.10 -> 0.60)
# 3. Better Treatment (alpha 0.015 -> 0.05) <--- CRITICAL CHANGE
EXP_C_AFTER_PARAMS = {
    'beta': 0.35,
    'gamma': 0.25,
    'sigma': 0.60,
    'alpha': 0.05,
    'eta': 0.10,
    'delta': 0.005,
    'tau': 0.15,
    'theta': 0.30,
    'rho': 0.05
}

EXP_C_CONFIG = {
    "p_supp_before": 0.01,
    "p_supp_after": 0.06,
}

# Fraction of initially S who has minimal susceptibility (susceptibility 0); rest have 1
P_LOW_SUSCEPTIBILITY = 0.25

def init_grid(N, p_supp, p_init_E, p_low_susceptibility=P_LOW_SUSCEPTIBILITY):
    grid = np.zeros((N, N), dtype=int)
    total = N * N
    grid.flat[np.random.choice(total, int(total * p_supp), replace=False)] = P
    avail = np.where(grid.flat == S)[0]
    grid.flat[np.random.choice(avail, int(total * p_init_E), replace=False)] = E
    # Per-cell susceptibility: S with prob p_low_susceptibility get 0 (non-drinker), else 1
    susceptibility = np.ones((N, N), dtype=float)
    mask_S = (grid == S)
    r = np.random.random((N, N))
    susceptibility[mask_S] = (r[mask_S] >= p_low_susceptibility).astype(float)
    return grid, np.zeros((N, N), dtype=int), susceptibility

def add_support(grid, target_p):
    N = grid.shape[0]
    curr_p = np.sum(grid == P) / (N*N)
    if curr_p >= target_p: return grid

    needed = int((target_p * N*N) - np.sum(grid == P))
    avail = np.where(grid.flat == S)[0]
    if len(avail) > 0:
        idx = np.random.choice(avail, min(len(avail), needed), replace=False)
        grid.flat[idx] = P
    return grid

def run_set(name, params, p_supp, intervention=False):
    print(f"Running {name}...")
    history = np.zeros((SHARED['runs'], SHARED['steps'], 4))  # R, E, C, T

    for r in range(SHARED['runs']):
        grid, k, susceptibility = init_grid(SHARED['N'], p_supp, EXP_A['p_init_E'])
        current_params = params.copy()

        for t in range(SHARED['steps']):
            if intervention and t == SHARED['implement_at']:
                grid = add_support(grid, EXP_C_CONFIG['p_supp_after'])
                current_params = EXP_C_AFTER_PARAMS.copy()

            history[r, t, 0] = np.sum(grid == R)
            history[r, t, 1] = np.sum(grid == E)
            history[r, t, 2] = np.sum(grid == C)
            history[r, t, 3] = np.sum(grid == T)
            grid, k = evolve(grid, k, current_params, susceptibility)

    # Save with names expected by plot_figures.py / plot_results.py
    out = {"ExpA": "experiment_a_results.npy", "ExpB": "experiment_b_results.npy", "ExpC": "experiment_c_results.npy"}[name]
    np.save(out, history)
    print(f"Saved {out}")

if __name__ == "__main__":
    run_set("ExpA", EXP_A['params'], EXP_A['p_supp'])
    run_set("ExpB", EXP_B['params'], EXP_B['p_supp'])
    run_set("ExpC", EXP_B['params'], EXP_C_CONFIG['p_supp_before'], intervention=True)