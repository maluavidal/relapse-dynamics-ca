import numpy as np
import matplotlib.pyplot as plt

# Load Data (same filenames as run_experiments.py)
data_A = np.load("experiment_a_results.npy")
data_B = np.load("experiment_b_results.npy")
data_C = np.load("experiment_c_results.npy")

steps = np.arange(data_A.shape[1])
N2 = 100 * 100

def get_mean_std(data, idx):
    return np.mean(data[:,:,idx], axis=0)/N2*100, np.std(data[:,:,idx], axis=0)/N2*100

# --- FIGURE 1: Experiment A (Success) ---
plt.figure(figsize=(8, 5))
mR, sR = get_mean_std(data_A, 0)
mC, sC = get_mean_std(data_A, 2)
mT, sT = get_mean_std(data_A, 3) if data_A.shape[2] > 3 else (np.zeros_like(steps, dtype=float), np.zeros_like(steps, dtype=float))

plt.plot(steps, mR, 'r-', lw=2, label='Relapsed (R)')
plt.fill_between(steps, mR-sR, mR+sR, color='r', alpha=0.1)
plt.plot(steps, mC, 'g-', lw=2, label='Controlled (C)')
plt.fill_between(steps, mC-sC, mC+sC, color='g', alpha=0.1)
if data_A.shape[2] > 3:
    plt.plot(steps, mT, color='magenta', lw=1.5, linestyle='--', label='Treatment (T)')
    plt.fill_between(steps, mT-sT, mT+sT, color='magenta', alpha=0.08)

plt.title("Figure 1: Baseline Dynamics (Experiment A)")
plt.xlabel("Time Steps")
plt.ylabel("% Population")
plt.legend()
plt.grid(alpha=0.3)
plt.savefig("Fig1_Success.png", dpi=300)
print("Saved Fig1_Success.png")

# --- FIGURE 2: Experiment B (Chronic) ---
plt.figure(figsize=(8, 5))
mR, sR = get_mean_std(data_B, 0)
mC, sC = get_mean_std(data_B, 2)
mT, sT = get_mean_std(data_B, 3) if data_B.shape[2] > 3 else (np.zeros_like(steps, dtype=float), np.zeros_like(steps, dtype=float))

plt.plot(steps, mR, 'r-', lw=2, label='Relapsed (R)')
plt.fill_between(steps, mR-sR, mR+sR, color='r', alpha=0.1)
plt.plot(steps, mC, 'g-', lw=2, label='Controlled (C)')
plt.fill_between(steps, mC-sC, mC+sC, color='g', alpha=0.1)
if data_B.shape[2] > 3:
    plt.plot(steps, mT, color='magenta', lw=1.5, linestyle='--', label='Treatment (T)')
    plt.fill_between(steps, mT-sT, mT+sT, color='magenta', alpha=0.08)

plt.title("Figure 2: Chronic Dynamics (Experiment B)")
plt.xlabel("Time Steps")
plt.ylabel("% Population")
plt.legend()
plt.grid(alpha=0.3)
plt.savefig("Fig2_Chronic.png", dpi=300)
print("Saved Fig2_Chronic.png")

# --- FIGURE 3: Experiment C (Intervention) ---
plt.figure(figsize=(9, 5))
mR_B, _ = get_mean_std(data_B, 0)
mR_C, sR_C = get_mean_std(data_C, 0)

plt.plot(steps, mR_B, color='gray', lw=2, label='Baseline (No Intervention, Exp B)')
plt.plot(steps, mR_C, color='steelblue', lw=3, label='With Intervention (Support Added)')
plt.fill_between(steps, mR_C-sR_C, mR_C+sR_C, color='steelblue', alpha=0.2)

plt.axvline(x=150, color='k', linestyle='--', label='Intervention (t=150)')

plt.title("Figure 3: Impact of Support Injection (Experiment C)")
plt.xlabel("Time Steps")
plt.ylabel("% Relapsed Population")
plt.legend()
plt.grid(alpha=0.3)
plt.savefig("Fig3_Intervention.png", dpi=300)
print("Saved Fig3_Intervention.png")