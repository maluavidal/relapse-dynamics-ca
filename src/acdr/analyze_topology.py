import numpy as np
from scipy.ndimage import label
from rules import evolve, S, E, R, C, T, P
from run_experiments import init_grid, EXP_B, EXP_C_CONFIG, EXP_C_AFTER_PARAMS

# --- 1. FRACTAL DIMENSION (The Core Metric) ---
def fractal_dimension(grid, state=R):
    """
    Calculates Box-Counting Dimension (Minkowski-Bouligand).
    Slope of log(N_boxes) vs log(1/box_size).
    """
    # Create binary image of the state
    Z = (grid == state)
    if not np.any(Z):
        return np.nan # Topologically empty

    p = min(Z.shape)

    # Largest power of 2 less than image size
    n = 2**np.floor(np.log(p)/np.log(2))
    n = int(np.log(n)/np.log(2))

    # Build box sizes: 64, 32, 16, 8, 4, 2
    sizes = 2**np.arange(n, 1, -1)

    counts = []
    for size in sizes:
        # Split grid into boxes of 'size'
        # Count how many boxes contain at least one True
        coords = np.array(np.where(Z)).T
        hashed = np.floor(coords / size)
        unique_boxes = len(np.unique(hashed, axis=0))
        counts.append(unique_boxes)

    # Linear Regression to find D
    # log(N) = D * log(1/s) + c
    coeffs = np.polyfit(np.log(1/sizes), np.log(counts), 1)
    return coeffs[0]

# --- 2. CLUSTER ANALYSIS (The "Collapse" Metric) ---
def analyze_clusters(grid, state=R):
    """Returns number of clusters and size of the largest cluster."""
    # Scipy handles the graph traversal for us
    labeled_array, num_features = label(grid == state)
    if num_features == 0:
        return 0, 0

    unique, counts = np.unique(labeled_array, return_counts=True)
    # unique[0] is usually the background (0), so skip it
    cluster_sizes = counts[1:] 
    if len(cluster_sizes) == 0:
        return 0, 0

    return num_features, np.max(cluster_sizes)

# --- 3. RUN CHECK ---
if __name__ == "__main__":
    N = 100
    STEPS = 300

    print(f"--- ANALYZING TOPOLOGY (t={STEPS}) ---")

    # 1. Run Exp B (Chronic)
    print("\nRunning Experiment B (Chronic)...")
    grid_B, k, susc = init_grid(N, EXP_B["p_supp"], 0.02)
    for _ in range(STEPS):
        grid_B, k = evolve(grid_B, k, EXP_B["params"], susc)

    d_B = fractal_dimension(grid_B, state=R)
    n_clusters_B, max_size_B = analyze_clusters(grid_B, state=R)

    print(f"  Fractal Dimension (D): {d_B:.4f} (Target: ~1.63)")
    print(f"  Number of Clusters:    {n_clusters_B}")
    print(f"  Largest Cluster Size:  {max_size_B}")

    # 2. Run Exp C (Intervention)
    print("\nRunning Experiment C (Intervention at t=150)...")
    grid_C, k, susc = init_grid(N, 0.01, 0.02)
    params = EXP_B["params"].copy() # Start as B

    for t in range(STEPS):
        if t == 150:
            # Apply Intervention Logic
            params = EXP_C_AFTER_PARAMS.copy()
            # Add supports (simplified logic for this script)
            needed = int((0.06 * N * N) - np.sum(grid_C == P))
            avail = np.where(grid_C == S)[0]
            if len(avail) > 0:
                idx = np.random.choice(avail, min(len(avail), needed), replace=False)
                grid_C.flat[idx] = P

        grid_C, k = evolve(grid_C, k, params, susc)

    d_C = fractal_dimension(grid_C, state=R)
    n_clusters_C, max_size_C = analyze_clusters(grid_C, state=R)

    print(f"  Fractal Dimension (D): {d_C} (Target: nan)")
    print(f"  Number of Clusters:    {n_clusters_C}")
    print(f"  Largest Cluster Size:  {max_size_C}")