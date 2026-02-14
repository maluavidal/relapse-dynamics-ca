import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import matplotlib.patches as mpatches
from rules import evolve, S, E, R, C, T, P

# Experiment A (Ideal) — same as run_experiments.EXP_A
PARAMS = {
    "beta": 0.25,
    "gamma": 0.25,
    "sigma": 0.60,
    "alpha": 0.05,
    "eta": 0.10,
    "delta": 1e-4,
    "tau": 0.10,
    "theta": 0.25,
    "rho": 0.05,
}
CONFIG = {
    "N": 100,
    "p_supp": 0.05,
    "p_init_E": 0.02,
    "snapshots": [0, 20, 50, 100, 200, 400],
}
# Fraction of S with minimal susceptibility (non-drinkers); same as run_experiments
P_LOW_SUSCEPTIBILITY = 0.25

# Experiment C (Intervention) — same as run_experiments: B until t=150, then add support + after params
IMPLEMENT_AT = 150
EXP_C_BEFORE = {"beta": 0.35, "gamma": 0.25, "sigma": 0.10, "alpha": 0.015, "eta": 0.50, "delta": 0.005, "tau": 0.05, "theta": 0.15, "rho": 0.10}
EXP_C_AFTER = {"beta": 0.35, "gamma": 0.25, "sigma": 0.60, "alpha": 0.05, "eta": 0.10, "delta": 0.005, "tau": 0.15, "theta": 0.30, "rho": 0.05}
EXP_C_P_SUPP_BEFORE = 0.01
EXP_C_P_SUPP_AFTER = 0.06
SNAPSHOTS_C = [0, 50, 100, 150, 250, 400]  # bracket intervention at 150

def add_support(grid, target_p):
    N = grid.shape[0]
    total = N * N
    need = max(0, int(total * target_p) - np.sum(grid == P))
    if need <= 0:
        return grid
    avail = np.where(grid.flat == S)[0]
    need = min(need, len(avail))
    if need > 0:
        grid.flat[np.random.choice(avail, need, replace=False)] = P
    return grid

def init_sim_visual(N, p_supp, p_init_E, p_low_susceptibility=P_LOW_SUSCEPTIBILITY):
    """Initial condition matching run_experiments.init_grid: random E, plus susceptibility for S."""
    grid = np.zeros((N, N), dtype=int)
    total = N * N
    grid.flat[np.random.choice(total, int(total * p_supp), replace=False)] = P
    avail = np.where(grid.flat == S)[0]
    num_E = max(1, min(len(avail), int(total * p_init_E)))
    grid.flat[np.random.choice(avail, num_E, replace=False)] = E
    susceptibility = np.ones((N, N), dtype=float)
    mask_S = (grid == S)
    r = np.random.random((N, N))
    susceptibility[mask_S] = (r[mask_S] >= p_low_susceptibility).astype(float)
    return grid, susceptibility

def main():
    print("Generating Figure 3 (Spatial)...")
    fig, axes = plt.subplots(2, 3, figsize=(10, 6))
    axes = axes.flatten()

    # Colors: S=Dark, E=Yellow, R=Red, C=Green, P=Blue
    cmap = ListedColormap(["#1a1a2e", "#f4d03f", "#e74c3c", "#2ecc71", "#e377c2", "#3498db"])
    legend_patches = [
        mpatches.Patch(color="#1a1a2e", label="Susceptible"),
        mpatches.Patch(color="#f4d03f", label="Exposed"),
        mpatches.Patch(color="#e74c3c", label="Relapsed"),
        mpatches.Patch(color="#2ecc71", label="Recovered"),
        mpatches.Patch(color="#e377c2", label="Treatment"),
        mpatches.Patch(color="#3498db", label="Support"),
    ]

    grid, susceptibility = init_sim_visual(CONFIG['N'], CONFIG['p_supp'], CONFIG['p_init_E'])
    k_counts = np.zeros((CONFIG['N'], CONFIG['N']), dtype=int)

    plot_idx = 0
    for t in range(401):
        if t in CONFIG['snapshots']:
            ax = axes[plot_idx]
            ax.imshow(grid, cmap=cmap, vmin=0, vmax=5, interpolation="nearest")
            ax.set_title(f"t = {t}")
            ax.axis("off")
            plot_idx += 1
        grid, k_counts = evolve(grid, k_counts, PARAMS, susceptibility)

    fig.legend(handles=legend_patches, loc="lower center", ncol=6, frameon=False)
    plt.suptitle("Spatial Dynamics of Recovery (Experiment A)")
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.1)
    plt.savefig("Fig3_Spatial.png", dpi=300)
    print("Saved Fig3_Spatial.png")
    plt.show()

    # --- Figure 4: Experiment C (Intervention) — spatial before/after ---
    print("Generating Figure 4 (Spatial, Intervention)...")
    fig2, axes2 = plt.subplots(2, 3, figsize=(10, 6))
    axes2 = axes2.flatten()

    grid, susceptibility = init_sim_visual(CONFIG["N"], EXP_C_P_SUPP_BEFORE, CONFIG["p_init_E"])
    k_counts = np.zeros((CONFIG["N"], CONFIG["N"]), dtype=int)
    params = dict(EXP_C_BEFORE)
    snapshot_set = set(SNAPSHOTS_C)
    max_t = max(SNAPSHOTS_C)
    plot_idx = 0

    for t in range(max_t + 1):
        if t == IMPLEMENT_AT:
            grid = add_support(grid, EXP_C_P_SUPP_AFTER)
            params = dict(EXP_C_AFTER)
        if t in snapshot_set:
            ax = axes2[plot_idx]
            ax.imshow(grid, cmap=cmap, vmin=0, vmax=5, interpolation="nearest")
            ax.set_title(f"t = {t}" + (" (intervention)" if t == IMPLEMENT_AT else ""))
            ax.axis("off")
            plot_idx += 1
        grid, k_counts = evolve(grid, k_counts, params, susceptibility)

    fig2.legend(handles=legend_patches, loc="lower center", ncol=6, frameon=False)
    plt.suptitle("Spatial Dynamics: Experiment C (Intervention at t=150)")
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.1)
    plt.savefig("Fig4_Spatial_Intervention.png", dpi=300)
    print("Saved Fig4_Spatial_Intervention.png")
    plt.show()

if __name__ == "__main__":
    main()