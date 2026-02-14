"""
Run CA to a chosen time, then compute box-counting fractal dimension of the R-set.
Usage: from src/acdr run  python3 compute_fractal_dimension.py
"""

import numpy as np
from rules import evolve
from run_experiments import init_grid, EXP_A, EXP_B, EXP_C_CONFIG, EXP_C_AFTER_PARAMS, add_support
from spatial_clusters import fractal_dimension_boxcounting, R_STATE

N = 100
steps = 300
p_init_E = 0.02


def run_to_time(params, p_supp, steps, intervention=False, implement_at=150):
    grid, k, susceptibility = init_grid(N, p_supp, p_init_E)
    current_params = params.copy()
    for t in range(steps):
        if intervention and t == implement_at:
            grid = add_support(grid, EXP_C_CONFIG["p_supp_after"])
            current_params = EXP_C_AFTER_PARAMS.copy()
        grid, k = evolve(grid, k, current_params, susceptibility)
    return grid


def report(grid, label, t):
    n_r = int(np.sum(grid == R_STATE))
    d = fractal_dimension_boxcounting(grid, state=R_STATE)
    if np.isnan(d):
        print(f"  {label} (t={t}): D = — (R-count = {n_r}; too few for box-counting)")
    else:
        print(f"  {label} (t={t}): D = {d:.4f}  (R-count = {n_r})")


if __name__ == "__main__":
    print("Fractal dimension (box-counting) of R-set.")
    print("(Exp A and C often have zero or very few R at t=300; reporting at t=100 too.)\n")

    print("Exp B (chronic) — high relapse sustained:")
    grid_b = run_to_time(EXP_B["params"], EXP_B["p_supp"], steps)
    report(grid_b, "Exp B", steps)

    print("Exp A (ideal) — relapse controlled by t=300:")
    grid_a = run_to_time(EXP_A["params"], EXP_A["p_supp"], steps)
    report(grid_a, "Exp A", steps)

    print("Exp C (intervention) — relapse reduced after t=150:")
    grid_c = run_to_time(EXP_B["params"], EXP_C_CONFIG["p_supp_before"], steps, intervention=True, implement_at=150)
    report(grid_c, "Exp C", steps)

    print("\nAt earlier time t=100 (all experiments still have non-trivial R):")
    grid_b_100 = run_to_time(EXP_B["params"], EXP_B["p_supp"], 100)
    grid_a_100 = run_to_time(EXP_A["params"], EXP_A["p_supp"], 100)
    grid_c_100 = run_to_time(EXP_B["params"], EXP_C_CONFIG["p_supp_before"], 100, intervention=False)
    report(grid_b_100, "Exp B", 100)
    report(grid_a_100, "Exp A", 100)
    report(grid_c_100, "Exp C (pre-intervention)", 100)
