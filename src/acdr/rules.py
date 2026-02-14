import numpy as np

# States
S = 0  # Susceptible
E = 1  # Exposed
R = 2  # Relapsed
C = 3  # Controlled (Recovered)
T = 4  # Treatment (Transient state, optional)
P = 5  # Support (Permanent Mentor)

def count_neighbors(grid, state):
    """Count Moore neighbors (8 surrounding cells) with a specific state."""
    N = grid.shape[0]
    padded = np.pad(grid, 1, mode='wrap')
    count = np.zeros_like(grid)

    # Check all 8 directions
    for dx in [-1, 0, 1]:
        for dy in [-1, 0, 1]:
            if dx == 0 and dy == 0:
                continue
            count += (padded[1+dx:N+1+dx, 1+dy:N+1+dy] == state)
    return count

def evolve(grid, relapse_counts, p, susceptibility=None):
    """
    Evolve the grid by one time step based on paper parameters p.
    susceptibility: optional (N,N) array in [0,1]; cells with 0 (e.g. non-drinkers) never transition S->E.
    """
    N = grid.shape[0]
    new_grid = grid.copy()

    # Pre-calculate neighbor counts to vectorize the process
    n_R = count_neighbors(grid, R)
    n_E = count_neighbors(grid, E)
    n_P = count_neighbors(grid, P)
    n_C = count_neighbors(grid, C)

    # 1. Susceptible (S) -> Exposed (E)
    # Social pressure from R and E neighbors (paper); mitigated by P or C (sigma)
    risk_neighbors = n_R + n_E
    pressure = np.where(risk_neighbors > 0, 1 - (1 - p['beta']) ** risk_neighbors, 0.0)
    has_support = (n_P > 0) | (n_C > 0)
    mitigation = 1 - p['sigma'] * has_support
    prob_S_to_E = pressure * mitigation
    if susceptibility is not None:
        prob_S_to_E = prob_S_to_E * susceptibility

    mask_S = (grid == S)
    random_S = np.random.random((N, N))
    new_grid[mask_S & (random_S < prob_S_to_E)] = E

    # 2. Exposed (E) -> Relapsed (R)
    # Fixed progression rate gamma; increment relapse count k when entering R (for exhaustion)
    mask_E = (grid == E)
    random_E = np.random.random((N, N))
    transition_E_R = mask_E & (random_E < p['gamma'])
    new_grid[transition_E_R] = R
    relapse_counts[transition_E_R] += 1

    # 3a. Relapsed (R) -> Treatment (T) with probability tau (optional)
    # 3b. Relapsed (R) -> Controlled (C) for R that did not go to T; P(R->C) = alpha/(1+eta*k)
    mask_R = (grid == R)
    tau = p.get('tau', 0)
    random_R1 = np.random.random((N, N))
    transition_R_T = mask_R & (random_R1 < tau)
    new_grid[transition_R_T] = T

    k = np.maximum(relapse_counts, 0)
    prob_R_to_C = p['alpha'] / (1 + p['eta'] * k)
    random_R2 = np.random.random((N, N))
    transition_R_C = mask_R & (~transition_R_T) & (random_R2 < prob_R_to_C)
    new_grid[transition_R_C] = C

    # 3c. Treatment (T) -> C (theta) or T -> R (rho); else stay T
    mask_T = (grid == T)
    theta = p.get('theta', 0)
    rho = p.get('rho', 0)
    if theta > 0 or rho > 0:
        random_T = np.random.random((N, N))
        transition_T_C = mask_T & (random_T < theta)
        transition_T_R = mask_T & (random_T >= theta) & (random_T < theta + rho)
        new_grid[transition_T_C] = C
        new_grid[transition_T_R] = R
        relapse_counts[transition_T_R] += 1  # dropout counts as relapse for exhaustion

    # 4. Controlled (C) -> Susceptible (S)
    # Relapse risk delta (loss of sobriety)
    mask_C = (grid == C)
    random_C = np.random.random((N, N))
    new_grid[mask_C & (random_C < p['delta'])] = S

    # Ensure Support nodes (P) never change
    new_grid[grid == P] = P

    return new_grid, relapse_counts
