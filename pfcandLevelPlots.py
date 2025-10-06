import os
import torch
import numpy as np
import matplotlib.pyplot as plt

# List of per-constituent variables
pfcand_keys = [
    'pfcand_isMu', 'pfcand_isEl', 'pfcand_isChargedHad', 'pfcand_isGamma',
    'pfcand_isNeutralHad', 'pfcand_e', 'pfcand_p', 'pfcand_theta', 'pfcand_phi',
    'pfcand_charge', 'pfcand_type', 'pfcand_erel', 'pfcand_erel_log',
    'pfcand_thetarel', 'pfcand_phirel', 'pfcand_dxy', 'pfcand_dz', 'pfcand_phi0',
    'pfcand_C', 'pfcand_ct', 'pfcand_dptdpt', 'pfcand_dxydxy', 'pfcand_dzdz',
    'pfcand_dphidphi', 'pfcand_detadeta', 'pfcand_dxydz', 'pfcand_dphidxy',
    'pfcand_phidz', 'pfcand_phictgtheta', 'pfcand_dxyctgtheta', 'pfcand_dlambdadz',
    'pfcand_cctgtheta', 'pfcand_phic', 'pfcand_dxyc', 'pfcand_cdz',
    'pfcand_btagSip2dVal', 'pfcand_btagSip2dSig', 'pfcand_btagSip3dVal',
    'pfcand_btagSip3dSig', 'pfcand_btagJetDistVal', 'pfcand_btagJetDistSig'
]

def plot_pfcands_with_ratios(pt_files, reference_keyword="thadthad", outdir="plots_pfcands", bins=100):
    os.makedirs(outdir, exist_ok=True)

    # Load all data
    data = {}
    for fname in pt_files:
        save_dict = torch.load(fname, map_location="cpu")
        base = os.path.splitext(os.path.basename(fname))[0]
        data[base] = {}
        for key in pfcand_keys:
            if key in save_dict:
                arr = save_dict[key]
                # Flatten (num_jets, n_const) → 1D
                vals = arr.flatten().numpy()
                # Filter out padding zeros (if needed)
                if np.issubdtype(vals.dtype, np.floating):
                    vals = vals[vals != 0]
                data[base][key] = vals

    # Find reference sample
    ref_name = None
    for name in data:
        if reference_keyword in name:
            ref_name = name
            break
    if ref_name is None:
        raise RuntimeError(f"No file contains '{reference_keyword}' for ratio reference!")

    print(f"Using {ref_name} as reference for ratios")

    # Loop over variables
    for key in pfcand_keys:
        plt.figure(figsize=(16,16),dpi=300)

        ax1 = plt.subplot(2,1,1)
        histograms = {}
        bin_edges = None

        # Reference hist
        ref_vals = data[ref_name][key]
        ref_hist, bin_edges = np.histogram(ref_vals, bins=bins, density=True)
        bin_centers = 0.5*(bin_edges[:-1] + bin_edges[1:])
        histograms[ref_name] = ref_hist
        ax1.step(bin_centers, ref_hist, where="mid", label=ref_name, linewidth=2)

        # Other samples
        for name, branch_dict in data.items():
            if name == ref_name or key not in branch_dict:
                continue
            vals = branch_dict[key]
            hist, _ = np.histogram(vals, bins=bin_edges, density=True)
            histograms[name] = hist
            ax1.step(bin_centers, hist, where="mid", label=name)

        ax1.set_ylabel("Normalized entries")
        ax1.set_title(key)
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
        ax2.set_xlabel(key)
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()
        outpath = os.path.join(outdir, f"{key}_comparison.png")
        plt.yscale('log')
        plt.savefig(outpath)
        plt.close()
        print(f"Saved {outpath}")


if __name__ == "__main__":
    pt_files = [f for f in os.listdir(".") if f.endswith(".pt")]
    plot_pfcands_with_ratios(pt_files, reference_keyword="thadthad", outdir="plots_pfcands", bins=20)
