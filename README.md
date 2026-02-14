# Alcohol Relapse Epidemiology — Cellular Automata Model

A **cellular automata (CA)** model used to study **alcohol relapse** from an epidemiology perspective. The grid represents a population where each cell can be in one of several states; local interactions (neighborhood contagion and support) drive transitions.

## States

| Code | State      | Meaning |
|------|------------|--------|
| S    | Susceptible | Not currently exposed to relapse triggers |
| E    | Exposed     | At risk of relapsing (social pressure) |
| R    | Relapsed    | In relapse |
| C    | Controlled  | In recovery; can provide support to neighbors |
| T    | Treatment   | Under active clinical intervention (transient) |
| P    | Support     | Mentor / peer support node (reduces contagion) |

## Dynamics

- **S → E**: Social contagion — probability increases with nearby R/E; **mitigated by support** (C or P in neighborhood, strength σ). Part of S has **minimal susceptibility** (non-drinkers): a fraction `P_LOW_SUSCEPTIBILITY` (default 25%) never transitions S→E.
- **E → R**: Relapse progression (γ); relapse count *k* is incremented when entering R (for exhaustion).
- **R → T**: With probability τ, relapsed individuals enter Treatment (alternative to natural recovery).
- **R → C**: For R that did not go to T: probability α/(1 + η·*k*) (exhaustion: harder after repeated relapses).
- **T → C**: With probability θ, treatment completes and the individual enters Controlled recovery.
- **T → R**: With probability ρ, treatment dropout (return to Relapsed); *k* is incremented.
- **C → S**: Degradation (δ) — recovery can fade.

Support nodes (P) and recovered individuals (C) reduce the chance that susceptible individuals become exposed. Initial condition: random placement of P and E; same logic in experiments and spatial plots.

## Experiments (paper parameter table)

| Parameter        | Symbol | Exp A (Ideal) | Exp B (Chronic) | Exp C (Intervention)      |
|-----------------|--------|---------------|-----------------|---------------------------|
| Social Pressure | β      | 0.25          | 0.35            | 0.35 (as B)               |
| Support Strength| σ      | 0.60          | 0.10            | 0.10 → 0.60               |
| Recovery Rate   | α      | 0.05          | 0.015           | 0.015 → 0.05              |
| Exhaustion      | η      | 0.10          | 0.50            | 0.50 → 0.10               |
| Relapse Risk    | δ      | 10⁻⁴          | 0.005           | 0.005                     |
| Support Density | p_supp | 5%            | 1%              | 1% → 6%                   |

- **A (Ideal)** — High support, low exhaustion; recovery spreads.
- **B (Chronic)** — Low support, high exhaustion; relapse persists.
- **C (Intervention)** — Start as B; at t=150 implement σ=0.60, α=0.05, η=0.10, and p_supp=6%.

## Project structure

```
relapse-dynamics-ca/
├── README.md
├── src/acdr/           # All Python scripts and generated data/figures (when run from here)
│   ├── rules.py
│   ├── run_experiments.py
│   ├── plot_figures.py
│   ├── plot_spatial.py
│   ├── plot_paper_figures.py
│   ├── compute_fractal_dimension.py
│   ├── analyze_topology.py
│   ├── experiment_*_results.npy   # produced by run_experiments.py
│   └── Fig*.png                   # produced by plotting scripts
└── docs/                # Paper, discrepancies note, optional copy of figures
```

## Files (`src/acdr/`)

- **`rules.py`** — CA transition rules, neighbor counting (periodic boundaries), optional per-cell susceptibility for S→E, and Treatment (T) transitions: R→T (τ), T→C (θ), T→R (ρ).
- **`run_experiments.py`** — Runs A, B, C; builds susceptibility at init; saves time series for R, E, C, T in `experiment_*_results.npy`.
- **`plot_figures.py`** — Time-series: Fig1_Success (A), Fig2_Chronic (B), Fig3_Intervention (C vs B). Reads `.npy` from current working directory (or `docs/` if not present).
- **`plot_spatial.py`** — Spatial snapshots: Fig3_Spatial (A), Fig4_Spatial_Intervention (C, snapshots around t=150). Uses `rules` and `run_experiments` for init/params.
- **`plot_paper_figures.py`** — Standalone script that runs its own simulations and produces paper-style figures: Fig4_Dynamics (prevalence over time), Fig5_Intervention (stack plot R/T/C with intervention at t=150), Fig6_Spatial (morphology of B vs C at end of run). Does not use the `.npy` files.
- **`compute_fractal_dimension.py`** — Runs the model to t=300 (and t=100) and reports the box-counting fractal dimension of the R-set for experiments A, B, C. **Requires** a `spatial_clusters` module providing `fractal_dimension_boxcounting` and `R_STATE`.
- **`analyze_topology.py`** — Runs one replicate of B and C to a fixed time, then reports fractal dimension (box-counting), number of R-clusters, and largest cluster size. Uses `scipy.ndimage.label` for cluster analysis. Useful for topology/collapse metrics.

## How to execute the workflow

All commands below are run from the **`src/acdr/`** directory. Scripts read and write files in the current working directory.

### Prerequisites

- **Python 3** with **NumPy** and **Matplotlib** (and **SciPy** for `analyze_topology.py`):
  ```bash
  pip install numpy matplotlib scipy
  ```

### Step 1: Run experiments

This runs the three experiments (A, B, C), each for 100 runs × 500 steps, and saves the time-series data:

```bash
cd src/acdr
python3 run_experiments.py
```

**Output:** `experiment_a_results.npy`, `experiment_b_results.npy`, `experiment_c_results.npy`

### Step 2: Plot time-series figures

Uses the `.npy` files from Step 1 to produce the main time-series plots:

```bash
python3 plot_figures.py
```

**Output:** `Fig1_Success.png`, `Fig2_Chronic.png`, `Fig3_Intervention.png`

### Step 3: Plot spatial snapshots

Generates spatial visualizations (independent of the `.npy` files; runs its own simulations):

```bash
python3 plot_spatial.py
```

**Output:** `Fig3_Spatial.png`, `Fig4_Spatial_Intervention.png`

### Additional scripts (optional)

- **Paper-style figures (no `.npy` needed):** From `src/acdr/`, run  
  `python3 plot_paper_figures.py`  
  to generate **Fig4_Dynamics.png**, **Fig5_Intervention.png**, **Fig6_Spatial.png** (standalone simulation + plotting).

- **Fractal dimension report:**  
  `python3 compute_fractal_dimension.py`  
  prints box-counting fractal dimension of the R-set for A, B, C at t=300 and t=100. Requires a `spatial_clusters` module.

- **Topology / cluster stats:**  
  `python3 analyze_topology.py`  
  prints fractal dimension, number of R-clusters, and largest cluster size for one run of B and C at a fixed time.

### Run the full pipeline

From the project root:

```bash
cd src/acdr && python3 run_experiments.py && python3 plot_figures.py && python3 plot_spatial.py
```

Optional: run `python3 plot_paper_figures.py` for Fig4–Fig6, and/or `python3 analyze_topology.py` for topology metrics. All figures and outputs are written in `src/acdr/`. Move or copy them (e.g. to `docs/`) if you want them in a different folder.
