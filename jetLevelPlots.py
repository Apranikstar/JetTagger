import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyarrow.parquet as pq
from matplotlib.backends.backend_pdf import PdfPages

# --- Jet-level branches ---
jet_branches = [
    'jet_mass', 'jet_p', 'jet_e', 'jet_phi', 'jet_theta', 'jet_pT',
    'jet_nnhad', 'jet_ngamma', 'jet_nchad', 'jet_nel', 'jet_nmu', 'jet_nconst'
]

# --- Short labels for legend ---
file_label_map = {
    "mgp8_pp_tleptlep": "tlept",
    "mgp8_pp_gg": "gg",
    "mgp8_pp_zhadzhad": "zhad",
    "mgp8_pp_cc": "cc",
    "mgp8_pp_whadwhad": "whad",
    "mgp8_pp_bb": "bb",
    "mgp8_pp_thadthad": "thad"
}

# --- Distinguishable colors ---
color_list = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"
]


def flatten_column(series):
    flat = []
    for entry in series:
        if isinstance(entry, (list, np.ndarray)):
            flat.extend(entry)
        elif pd.notna(entry):
            flat.append(entry)
    return np.array(flat, dtype=float) if flat else np.array([])


def auto_integer_bins(vals, max_bins=50):
    """
    Automatically determine integer bins:
    - Step = 1 if range <= max_bins
    - Otherwise, step = 2, 3,... to keep number of bins <= max_bins
    """
    vmin, vmax = np.nanmin(vals), np.nanmax(vals)
    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmin == vmax:
        return None
    vmin_int, vmax_int = int(np.floor(vmin)), int(np.ceil(vmax))
    n_bins = vmax_int - vmin_int
    if n_bins <= max_bins:
        step = 1
    else:
        step = int(np.ceil(n_bins / max_bins))
    bins = np.arange(vmin_int, vmax_int + step, step)
    return bins


def plot_jet_variables_with_ratios(parquet_files, reference_keyword="thadthad", outdir="plots_jets"):
    os.makedirs(outdir, exist_ok=True)
    data = {}

    # --- Load datasets ---
    for fname in parquet_files:
        print(f"📂 Loading {fname} ...")
        try:
            df = pq.read_table(fname).to_pandas()
        except Exception as e:
            print(f"❌ Failed to read {fname}: {e}")
            continue

        base = os.path.splitext(os.path.basename(fname))[0]
        data[base] = {}
        for branch in jet_branches:
            if branch in df.columns:
                vals = flatten_column(df[branch])
                if np.issubdtype(vals.dtype, np.floating):
                    vals = vals[np.isfinite(vals)]
                if len(vals) > 0:
                    data[base][branch] = vals

    if not data:
        raise RuntimeError("❌ No valid jet data loaded!")

    # --- Reference sample ---
    ref_name = next((name for name in data if reference_keyword in name), None)
    if ref_name is None:
        ref_name = list(data.keys())[0]
        print(f"⚠️ No '{reference_keyword}' found, using '{ref_name}' as reference instead.")
    else:
        print(f"✅ Using '{ref_name}' as reference for ratios.")

    # --- Prepare PDF output ---
    pdf_path = os.path.join(outdir, "jet_plots_all.pdf")
    with PdfPages(pdf_path) as pdf:
        print(f"📄 Creating PDF: {pdf_path}")

        # --- Loop over variables ---
        for branch in jet_branches:
            valid_arrays = [v[branch] for v in data.values() if branch in v and len(v[branch]) > 0]
            if not valid_arrays:
                print(f"⚠️ Skipping {branch}: no data available.")
                continue

            all_vals = np.concatenate(valid_arrays)
            bin_edges = auto_integer_bins(all_vals, max_bins=50)
            if bin_edges is None:
                print(f"⚠️ Skipping {branch}: invalid range or NaNs.")
                continue

            vmin, vmax = bin_edges[0], bin_edges[-1]

            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 12), dpi=300, gridspec_kw={'height_ratios': [3, 1]})
            histograms = {}

            # --- Top: filled histograms with 25% opacity and visible borders ---
            for i, (name, branch_dict) in enumerate(data.items()):
                if branch not in branch_dict:
                    continue
                vals = branch_dict[branch]
                hist, _ = np.histogram(vals, bins=bin_edges, density=True)
                histograms[name] = hist

                label = file_label_map.get(name.split("_PTmin")[0], name)
                color = color_list[i % len(color_list)]

                ax1.hist(
                    vals, bins=bin_edges, density=True,
                    facecolor=color, edgecolor=color, linewidth=1.4, alpha=0.25,
                    label=label
                )

            ax1.set_ylabel("Normalized entries")
            ax1.set_title(f"{branch} (range = [{vmin}, {vmax}])")
            ax1.legend()
            ax1.grid(True, alpha=0.3)
            ax1.set_yscale("log")

            # --- Bottom: ratio step plots ---
            ref_hist = histograms[ref_name]
            for i, (name, hist) in enumerate(histograms.items()):
                if name == ref_name:
                    continue
                ratio = np.divide(hist, ref_hist, out=np.ones_like(hist), where=ref_hist > 0)
                label = file_label_map.get(name.split("_PTmin")[0], name)
                color = color_list[i % len(color_list)]
                ax2.step(bin_edges[:-1], ratio, where='post', label=label, color=color, linewidth=1.8)

            ax2.axhline(1.0, color='black', linestyle='--')
            ax2.set_ylabel(f"Ratio to {file_label_map.get(ref_name.split('_PTmin')[0], ref_name)}")
            ax2.set_xlabel(branch)
            ax2.grid(True, alpha=0.3)
            ax2.legend()

            plt.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

            print(f"💾 Added {branch} to PDF (bins={len(bin_edges)-1}, range={vmin}–{vmax}, step={bin_edges[1]-bin_edges[0]})")

    print(f"✅ All plots saved to {pdf_path}")


if __name__ == "__main__":
    parquet_files = [f for f in os.listdir(".") if f.endswith(".parquet")]
    if not parquet_files:
        print("❌ No .parquet files found in current directory!")
    else:
        plot_jet_variables_with_ratios(
            parquet_files,
            reference_keyword="thadthad",
            outdir="plots_jets"
        )
