import uproot
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import os

# --- Config ---
input_dir = "/eos/user/h/hfatehi/JetTagger/TESTPIPELINE/"
output_dir = "/eos/user/h/hfatehi/JetTagger/Plots_2D_Grid_HiRes/"
os.makedirs(output_dir, exist_ok=True)
tree_name = "events"

# --- Sample mapping ---
sample_map = {
    "_gg_": "gg",
    "_bb_": "bb",
    "_cc_": "cc",
    "_thadthad_": "thad",
    "_tleptlep_": "tlept",
    "_whadwhad_": "whad",
    "_zhadzhad_": "zhad",
    "_uuddss_": "uuddss"
}

# --- Load all samples ---
samples_data = {}
for fname in os.listdir(input_dir):
    if not fname.endswith(".root"):
        continue

    file_path = os.path.join(input_dir, fname)
    sample_label = next((label for pat, label in sample_map.items() if pat in fname), "unknown")
    print(f"Reading {file_path} -> {sample_label}")

    leading_pts = []
    subleading_pts = []

    with uproot.open(file_path) as f:
        tree = f[tree_name]
        jet_pts = tree["jet_pT"].array(library="np")
        for jets in jet_pts:
            if len(jets) >= 2:
                leading_pts.append(jets[0])
                subleading_pts.append(jets[1])

    leading_pts = np.array(leading_pts)
    subleading_pts = np.array(subleading_pts)
    sum_pts = leading_pts + subleading_pts

    samples_data[sample_label] = (sum_pts, leading_pts)

# --- Determine grid size ---
n_samples = len(samples_data)
ncols = 3
nrows = int(np.ceil(n_samples / ncols))

# --- Determine global axis limits ---
all_x = np.concatenate([v[0] for v in samples_data.values()])
all_y = np.concatenate([v[1] for v in samples_data.values()])
x_min, x_max = all_x.min(), all_x.max()
y_min, y_max = all_y.min(), all_y.max()

# --- Create figure with subplots (larger) ---
fig, axes = plt.subplots(nrows, ncols, figsize=(7*ncols,5*nrows), constrained_layout=True)  # increased size
axes = axes.flatten()

for ax, (sample, (x, y)) in zip(axes, samples_data.items()):
    h = ax.hist2d(
        x, y, bins=100,
        range=[[x_min, x_max], [y_min, y_max]],
        norm=mcolors.LogNorm(),
        cmap='viridis'
    )
    ax.set_title(sample)
    ax.set_xlabel("Sum of Leading + Subleading Jet $p_T$ [GeV]")
    ax.set_ylabel("Leading Jet $p_T$ [GeV]")

# Hide any unused axes
for i in range(len(samples_data), len(axes)):
    fig.delaxes(axes[i])

# Add a single colorbar for all subplots
cbar = fig.colorbar(h[3], ax=axes, orientation='vertical', fraction=0.02, pad=0.02)
cbar.set_label("Counts")

plt.suptitle("2D Histograms: Leading vs Sum of Jet $p_T$ per Sample", fontsize=16)

# --- Save figure with double resolution ---
plt.savefig(os.path.join(output_dir, "all_samples_2Dgrid_hi_res.png"), dpi=400)  # doubled dpi
plt.close()

print(f"✅ High-resolution 2D histogram grid saved to {output_dir}")
