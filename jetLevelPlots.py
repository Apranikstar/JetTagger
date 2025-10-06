import os
import torch
import numpy as np
import matplotlib.pyplot as plt

# List of jet-level branches to plot
jet_branches = [
    'jet_mass', 'jet_p', 'jet_e', 'jet_phi', 'jet_theta', 'jet_pT',
    'jet_nnhad', 'jet_ngamma', 'jet_nchad', 'jet_nel', 'jet_nmu', 'jet_nconst'
]

def plot_with_ratios(pt_files, reference_keyword="thadthad", outdir="plots", bins=100):
    os.makedirs(outdir, exist_ok=True)

    # Load all files into memory
    data = {}
    for fname in pt_files:
        save_dict = torch.load(fname, map_location="cpu")
        base = os.path.splitext(os.path.basename(fname))[0]
        data[base] = {branch: save_dict[branch].numpy() for branch in jet_branches if branch in save_dict}

    # Pick reference sample (must contain reference_keyword in name)
    ref_name = None
    for name in data:
        if reference_keyword in name:
            ref_name = name
            break
    if ref_name is None:
        raise RuntimeError(f"No file contains '{reference_keyword}' for ratio reference!")

    print(f"Using {ref_name} as reference for ratios")

    # Now plot each variable
    for branch in jet_branches:
        plt.figure(figsize=(16,16),dpi=300)
        
        # Main histo (upper panel)
        ax1 = plt.subplot(2,1,1)
        histograms = {}
        bin_edges = None

        # Build reference histogram first
        ref_vals = data[ref_name][branch].flatten()
        ref_hist, bin_edges = np.histogram(ref_vals, bins=bins, density=True)
        histograms[ref_name] = ref_hist
        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

        # Plot reference
        ax1.step(bin_centers, ref_hist, where="mid", label=ref_name, linewidth=2)

        # Other samples
        for name, branch_dict in data.items():
            if name == ref_name or branch not in branch_dict:
                continue
            vals = branch_dict[branch].flatten()
            hist, _ = np.histogram(vals, bins=bin_edges, density=True)
            histograms[name] = hist
            ax1.step(bin_centers, hist, where="mid", label=name)

        ax1.set_ylabel("Normalized entries")
        ax1.set_title(branch)
        ax1.legend()
        ax1.grid(True, alpha=0.3)

        # Ratio panel
        ax2 = plt.subplot(2,1,2, sharex=ax1)
        ref_hist = histograms[ref_name]
        for name, hist in histograms.items():
            if name == ref_name:
                continue
            ratio = np.divide(hist, ref_hist, out=np.ones_like(hist), where=ref_hist>0)
            ax2.step(bin_centers, ratio, where="mid", label=name)
        ax2.axhline(1.0, color="black", linestyle="--")
        ax2.set_ylabel("Ratio to " + ref_name)
        ax2.set_xlabel(branch)
        ax2.grid(True, alpha=0.3)
        plt.yscale('log')
        plt.tight_layout()
        outpath = os.path.join(outdir, f"{branch}_comparison.png")
        plt.savefig(outpath)
        #plt.show()
        plt.close()
        print(f"Saved {outpath}")


if __name__ == "__main__":
    pt_files = [f for f in os.listdir(".") if f.endswith(".pt")]
    plot_with_ratios(pt_files, reference_keyword="thadthad", outdir="plots", bins=20)
