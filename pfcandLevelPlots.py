import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from itertools import cycle
from matplotlib.colors import to_rgb
import os
from matplotlib.backends.backend_pdf import PdfPages

# Variables of interest
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

# Short labels for legend
file_label_map = {
    "mgp8_pp_tleptlep": "tlept",
    "mgp8_pp_gg": "gg",
    "mgp8_pp_zhadzhad": "zhad",
    "mgp8_pp_cc": "cc",
    "mgp8_pp_whadwhad": "whad",
    "mgp8_pp_bb": "bb",
    "mgp8_pp_thadthad": "thad"
}

def flatten_column(col):
    if len(col) == 0:
        return np.array([])
    first = col.iloc[0]
    if isinstance(first, (list, np.ndarray)):
        return np.concatenate([np.array(x).ravel() for x in col if isinstance(x, (list, np.ndarray))])
    return col.to_numpy()

def darken_color(color, factor=0.6):
    rgb = np.array(to_rgb(color))
    return tuple(rgb * factor)

def auto_integer_bins(vals, max_bins=50):
    vmin, vmax = np.nanmin(vals), np.nanmax(vals)
    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmin == vmax:
        return None
    vmin_int, vmax_int = int(np.floor(vmin)), int(np.ceil(vmax))
    n_bins = vmax_int - vmin_int
    step = 1 if n_bins <= max_bins else int(np.ceil(n_bins / max_bins))
    bins = np.arange(vmin_int, vmax_int + step, step)
    return bins

def plot_pfcands_with_ratios_pdf(parquet_files, reference_keyword="thad", outdir="plots_pfcands"):
    os.makedirs(outdir, exist_ok=True)
    data = {}

    # --- Load data ---
    for fname in parquet_files:
        print(f"📂 Loading {fname} ...")
        df = pd.read_parquet(fname)
        base = os.path.splitext(os.path.basename(fname))[0]

        short_label = next((v for k,v in file_label_map.items() if k in base), base)
        data[short_label] = {}
        for key in pfcand_keys:
            if key in df.columns:
                vals = flatten_column(df[key])
                if np.issubdtype(vals.dtype, np.floating):
                    vals = vals[np.isfinite(vals)]
                if len(vals) > 0:
                    data[short_label][key] = vals

    # --- Reference ---
    ref_name = next((name for name in data if reference_keyword in name), None)
    if ref_name is None:
        raise RuntimeError(f"No file contains '{reference_keyword}' for ratio reference!")
    print(f"✅ Using '{ref_name}' as reference for ratios")

    color_cycle = cycle(plt.cm.tab10.colors + plt.cm.Set2.colors)

    pdf_path = os.path.join(outdir, "pfcands_all.pdf")
    with PdfPages(pdf_path) as pdf:
        for key in pfcand_keys:
            valid_arrays = [v[key] for v in data.values() if key in v and len(v[key]) > 0]
            if not valid_arrays:
                continue

            all_vals = np.concatenate(valid_arrays)
            bin_edges = auto_integer_bins(all_vals, max_bins=50)
            if bin_edges is None:
                continue
            vmin, vmax = bin_edges[0], bin_edges[-1]

            fig, (ax1, ax2) = plt.subplots(2,1, figsize=(12,10), dpi=300, gridspec_kw={'height_ratios':[3,1]})
            histograms = {}
            colors = {name: next(color_cycle) for name in data.keys()}

            # --- Top: filled histograms with edges ---
            for name, branch_dict in data.items():
                if key not in branch_dict:
                    continue
                vals = branch_dict[key]
                hist, _ = np.histogram(vals, bins=bin_edges, density=True)
                histograms[name] = hist

                color = colors[name]
                edgecolor = 'black' if name == ref_name else darken_color(color, 0.6)
                alpha = 0.5
                linewidth = 1.8 if name == ref_name else 1.2
                label = f"{name} (ref)" if name==ref_name else name

                ax1.hist(vals, bins=bin_edges, density=True,
                         histtype='stepfilled', facecolor=color, edgecolor=edgecolor,
                         alpha=alpha, linewidth=linewidth, label=label)

            ax1.set_ylabel("Normalized entries")
            ax1.set_title(f"{key} (range = [{vmin}, {vmax}])")
            ax1.grid(True, alpha=0.3)
            ax1.set_yscale("log")
            ax1.legend(fontsize=9)

            # --- Bottom: ratio plots ---
            ref_hist = histograms[ref_name]
            for name, hist in histograms.items():
                if name == ref_name:
                    continue
                ratio = np.divide(hist, ref_hist, out=np.ones_like(hist), where=ref_hist>0)
                bin_centers = 0.5*(bin_edges[:-1] + bin_edges[1:])
                ax2.plot(bin_centers, ratio, label=name, color=colors[name], linewidth=1.5)

            ax2.axhline(1.0, color='black', linestyle='--')
            ax2.set_ylabel(f"Ratio to {ref_name}")
            ax2.set_xlabel(key)
            ax2.grid(True, alpha=0.3)
            ax2.legend(fontsize=9)

            plt.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

            print(f"💾 Added {key} to PDF (bins={len(bin_edges)-1}, range={vmin}–{vmax}, step={bin_edges[1]-bin_edges[0]})")

    print(f"✅ All plots saved to {pdf_path}")


if __name__ == "__main__":
    parquet_files = [f for f in os.listdir(".") if f.endswith(".parquet")]
    plot_pfcands_with_ratios_pdf(parquet_files, reference_keyword="thad", outdir="plots_pfcands")
