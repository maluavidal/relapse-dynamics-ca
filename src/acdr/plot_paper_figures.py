import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# --- 1. CONFIGURATION ---
N = 100
STEPS = 350 # Run slightly past 300 to show stability
S, E, R, C, T, P = 0, 1, 2, 3, 4, 5
STATES = [S, E, R, C, T, P]
COLORS = ['#f0f0f0', '#ffcc00', '#d62728', '#2ca02c', '#1f77b4', '#000000']
# Colors: White(S), Yellow(E), Red(R), Green(C), Blue(T), Black(P)
CMAP = mcolors.ListedColormap(COLORS)

# --- 2. EXACT PARAMETERS (MATCHING TABLE 1) ---
PARAMS = {
    # Exp A: Ideal (High Access, Low Dropout)
    'A': {'beta': 0.25, 'gamma': 0.25, 'sigma': 0.60, 'tau': 0.15, 'theta': 0.30, 'rho': 0.02, 'alpha': 0.05, 'eta': 0.10},

    # Exp B: Chronic (Low Access, High Exhaustion)
    'B': {'beta': 0.35, 'gamma': 0.25, 'sigma': 0.10, 'tau': 0.02, 'theta': 0.10, 'rho': 0.10, 'alpha': 0.015, 'eta': 0.50},

    # Exp C (Post-Intervention): High Access, High Support
    'C_Post': {'beta': 0.35, 'gamma': 0.25, 'sigma': 0.60, 'tau': 0.20, 'theta': 0.30, 'rho': 0.05, 'alpha': 0.05, 'eta': 0.10}
}

# --- 3. CORE SIMULATION LOGIC ---
def get_neighbors(grid, r=1):
    padded = np.pad(grid, r, mode='wrap')
    neighs = np.zeros_like(grid, dtype=int)
    for i in range(2*r + 1):
        for j in range(2*r + 1):
            if i == r and j == r: continue
            neighs += (padded[i:i+N, j:j+N] == R) | (padded[i:i+N, j:j+N] == E)
    return neighs

def run_simulation(mode='A'):
    # Init Grid
    grid = np.zeros((N, N), dtype=int)
    support_mask = np.zeros((N, N), dtype=bool)
    exhaustion = np.zeros((N, N), dtype=int)

    # Init Support Nodes (P)
    p_supp = 0.05 if mode == 'A' else 0.01
    num_supp = int(N*N * p_supp)
    supp_idx = np.random.choice(N*N, num_supp, replace=False)
    grid.flat[supp_idx] = P
    support_mask.flat[supp_idx] = True

    # Init Exposure (E)
    avail = np.where(grid.flat == S)[0]
    grid.flat[np.random.choice(avail, int(N*N * 0.02), replace=False)] = E

    history = np.zeros((STEPS, 6))
    current_grid = grid.copy()

    for t in range(STEPS):
        # Handle Parameter Switching for Exp C
        if mode == 'C':
            if t < 150:
                p = PARAMS['B'] # Start as Chronic
            elif t == 150:
                p = PARAMS['C_Post'] # Switch params
                # Add extra support nodes (Intervention)
                candidates = np.where((current_grid != P) & (current_grid != T))[0]
                new_supp = np.random.choice(candidates, int(N*N * 0.05), replace=False)
                current_grid.flat[new_supp] = P
                support_mask.flat[new_supp] = True
            else:
                p = PARAMS['C_Post']
        else:
            p = PARAMS[mode]

        # Record History
        counts = [np.sum(current_grid == s) for s in STATES]
        history[t] = counts

        # --- CALC PROBABILITIES ---
        neigh_risk = get_neighbors(current_grid, r=1)
        prob_S_E = 1 - (1 - p['beta'])**neigh_risk

        # Mitigation Mask (Simplified efficient check)
        is_protected = np.zeros((N, N), dtype=bool)
        # Check P (r=2)
        padded_P = np.pad(current_grid == P, 2, mode='wrap')
        for i in range(5):
            for j in range(5):
                is_protected |= padded_P[i:i+N, j:j+N]
        # Check C (r=1)
        padded_C = np.pad(current_grid == C, 1, mode='wrap')
        for i in range(3):
            for j in range(3):
                is_protected |= padded_C[i:i+N, j:j+N]

        prob_S_E[is_protected] *= (1 - p['sigma'])

        # --- TRANSITIONS ---
        rand = np.random.random((N, N))
        next_grid = current_grid.copy()

        # S -> E
        mask_S = (current_grid == S)
        next_grid[mask_S & (rand < prob_S_E)] = E

        # E -> R
        mask_E = (current_grid == E)
        next_grid[mask_E & (rand < p['gamma'])] = R

        # R -> T or C
        mask_R = (current_grid == R)
        to_treatment = mask_R & (rand < p['tau'])
        next_grid[to_treatment] = T

        remaining_R = mask_R & ~to_treatment
        prob_R_C = p['alpha'] / (1 + p['eta'] * exhaustion)
        to_recovery = remaining_R & (rand < prob_R_C)
        next_grid[to_recovery] = C

        exhaustion[mask_R] += 1 # Increment exhaustion

        # T -> C or R
        mask_T = (current_grid == T)
        success = mask_T & (rand < p['theta'])
        next_grid[success] = C
        dropout = mask_T & (rand > (1 - p['rho'])) & ~success
        next_grid[dropout] = R

        # C -> S (Relapse)
        mask_C = (current_grid == C)
        next_grid[mask_C & (rand < 0.005)] = S

        # Support stays Support
        next_grid[support_mask] = P

        current_grid = next_grid

    return history, current_grid

# --- 4. EXECUTION ---
print("Simulating Exp A...")
hist_A, grid_A = run_simulation('A')
print("Simulating Exp B...")
hist_B, grid_B = run_simulation('B')
print("Simulating Exp C...")
hist_C, grid_C = run_simulation('C')

# --- 5. PLOTTING ---

# FIGURE 1: DYNAMICS
plt.figure(figsize=(10, 5))
plt.plot(hist_A[:, R]/100, label='Experiment A (Ideal)', color='green', linewidth=2)
plt.plot(hist_B[:, R]/100, label='Experiment B (Chronic)', color='red', linewidth=2)
plt.axhline(y=30, color='gray', linestyle='--', alpha=0.5, label='Endemic Threshold')
plt.title('Figure 4: Relapse Prevalence Over Time')
plt.xlabel('Time Steps')
plt.ylabel('% Population in Relapse')
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig('Fig4_Dynamics.png', dpi=300)
print("Saved Fig4_Dynamics.png")

# FIGURE 2: INTERVENTION STACKPLOT
plt.figure(figsize=(10, 5))
x = np.arange(STEPS)
y = [hist_C[:, R], hist_C[:, T], hist_C[:, C]]
labels = ['Relapsed (R)', 'Treatment (T)', 'Controlled (C)']
colors = ['#d62728', '#1f77b4', '#2ca02c'] 
plt.stackplot(x, y, labels=labels, colors=colors, alpha=0.85)
plt.axvline(x=150, color='black', linestyle='--', linewidth=2, label='Intervention (t=150)')
plt.title('Figure 5: The "Push-Pull" Intervention Effect')
plt.xlabel('Time Steps')
plt.ylabel('Number of Agents')
plt.legend(loc='upper left')
plt.savefig('Fig5_Intervention.png', dpi=300)
print("Saved Fig5_Intervention.png")

# FIGURE 3: SPATIAL TOPOLOGY
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

# Plot Exp B (Chronic)
cmap = mcolors.ListedColormap(['white', 'orange', 'red', 'green', 'blue', 'black'])
ax1.imshow(grid_B, cmap=cmap, vmin=0, vmax=5)
ax1.set_title('Experiment B (Chronic)\nFractal Dust (D ≈ 1.58)')
ax1.axis('off')

# Plot Exp C (Intervention)
ax2.imshow(grid_C, cmap=cmap, vmin=0, vmax=5)
ax2.set_title('Experiment C (Intervention)\nTopological Collapse (Clusters)')
ax2.axis('off')

# Legend
patches = [plt.Rectangle((0,0),1,1, color=COLORS[i]) for i in [2,3,4,5]]
labels = ['Relapsed', 'Controlled', 'Treatment', 'Support']
fig.legend(patches, labels, loc='lower center', ncol=4, bbox_to_anchor=(0.5, 0.02))

plt.suptitle('Figure 6: Spatial Morphology of Addiction Dynamics')
plt.tight_layout(rect=[0, 0.05, 1, 1])
plt.savefig('Fig6_Spatial.png', dpi=300)
print("Saved Fig6_Spatial.png")