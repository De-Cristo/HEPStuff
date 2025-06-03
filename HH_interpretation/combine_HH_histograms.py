import os
import uproot
import numpy as np
import matplotlib.pyplot as plt
from hist import Hist
import HH_Combine  # Make sure HH_Combine.py is in the same directory or PYTHONPATH

# --- Configuration ---
ROOT_FILE_DIR = "/eos/user/l/lichengz/HHbbWW/outputDir_boosted_DL_0603/results/"
FILE_NAME_TEMPLATE = "ggH_bbww_2L2Nu_kl-{kl_str}_kt-1p00_c2-0p00_2022.root"
TREE_NAME = "DL_resolved_2b_emu_ml_vars"
PHYSICS_BRANCH = "InvM_met_diLep_diB"
WEIGHT_BRANCH = "weight"
BINS = 50
LOW = 0
HIGH = 1000

BENCHMARKS = [
    {"kl": 0.00, "kt": 1.0, "kl_str": "0p00", "label": "kl = 0.00", "xs": 0.069725},
    {"kl": 2.45, "kt": 1.0, "kl_str": "2p45", "label": "kl = 2.45", "xs": 0.013124},
    {"kl": 5.00, "kt": 1.0, "kl_str": "5p00", "label": "kl = 5.00", "xs": 0.091172},
]

def load_weighted_histogram(file_path, tree_name, x_branch, weight_branch, bins, low, high, label):
    try:
        with uproot.open(file_path) as file:
            tree = file[tree_name]
            arrays = tree.arrays([x_branch, weight_branch], library="np")
            x = arrays[x_branch]
            w = arrays[weight_branch]

            if len(x) == 0 or len(w) == 0 or len(x) != len(w):
                print(f"Warning: Empty or mismatched data in {file_path}")
                return None

            h = Hist.new.Reg(bins, low, high, name="x", label=x_branch).Double()
            h.label = label
            h.fill(x, weight=w)
            return h
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return None

def build_combined_histogram(target_kl, benchmarks, tree_name, x_branch, weight_branch, bins, low, high):
    samples_for_weights = [{"kl": b["kl"], "kt": b["kt"], "xs": b["xs"], "label": b["label"]} for b in benchmarks]
    coeffs, _, _, _ = HH_Combine.calculate_ggf_parameters(target_kl=target_kl, samples_data=samples_for_weights, target_kt=1.0)

    print("\nWeights for combination:")
    for b, w in zip(benchmarks, coeffs):
        print(f"  {b['label']}: weight = {w:.4f}")

    histograms = []
    for i, b in enumerate(benchmarks):
        file_path = os.path.join(ROOT_FILE_DIR, FILE_NAME_TEMPLATE.format(kl_str=b["kl_str"]))
        h = load_weighted_histogram(file_path, tree_name, x_branch, weight_branch, bins, low, high, label=b["label"])
        if h:
            histograms.append((h, coeffs[i]))

    if not histograms:
        raise RuntimeError("No valid histograms loaded for combination.")

    # Use axis from first histogram
    ref_axis = histograms[0][0].axes[0]
    combined = Hist(ref_axis)  # This creates the histogram
    combined.label = f"Combined kl = {target_kl:.2f}"

    for h, w in histograms:
        combined += h * w

    return combined, histograms

def main():
    try:
        target_kl = float(input("Enter target kl value (e.g. 1.0): "))
    except ValueError:
        print("Invalid input.")
        return


    # --- Plot and save individual benchmark histograms ---
    plt.figure(figsize=(10, 7))
    
    for benchmark in BENCHMARKS:
        file_path = os.path.join(ROOT_FILE_DIR, FILE_NAME_TEMPLATE.format(kl_str=benchmark["kl_str"]))
        h = load_weighted_histogram(
            file_path=file_path,
            tree_name=TREE_NAME,
            x_branch=PHYSICS_BRANCH,
            weight_branch=WEIGHT_BRANCH,
            bins=BINS,
            low=LOW,
            high=HIGH,
            label=benchmark["label"]
        )
    
        if h is not None and h.sum() > 0:
            h.plot(histtype="step", lw=2, label=benchmark["label"])
        else:
            print(f"Skipping benchmark {benchmark['label']}: histogram not valid.")
    
    plt.xlabel(f"{PHYSICS_BRANCH} (GeV)")
    plt.ylabel("Events / Bin Width")
    plt.title("Benchmark Samples (Before Combination)")
    plt.legend()
    plt.grid(True, linestyle="--", alpha=0.7)
    plt.tight_layout()
    
    output_benchmark_plot = "benchmark_histograms.png"
    plt.savefig(output_benchmark_plot)
    print(f"Saved benchmark histograms plot: {output_benchmark_plot}")
    plt.close()

    combined_hist, components = build_combined_histogram(
        target_kl=target_kl,
        benchmarks=BENCHMARKS,
        tree_name=TREE_NAME,
        x_branch=PHYSICS_BRANCH,
        weight_branch=WEIGHT_BRANCH,
        bins=BINS,
        low=LOW,
        high=HIGH
    )

    # --- Plot ---
    plt.figure(figsize=(10, 7))
    for h, w in components:
        h_scaled = h * w
        h_scaled.plot(histtype="step", label=f"{h.label} × {w:.2f}")

    combined_hist.plot(histtype="step", lw=2, color="black", label=combined_hist.label)

    plt.xlabel(f"{PHYSICS_BRANCH} (GeV)")
    plt.ylabel("Events / Bin Width")
    plt.title(f"Synthetic HH Signal: kl = {target_kl:.2f}")
    plt.legend()
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.tight_layout()
    out_file = f"combined_hist_kl_{str(target_kl).replace('.', 'p')}.png"
    plt.savefig(out_file)
    print(f"Saved combined histogram plot: {out_file}")
    # plt.show()

if __name__ == "__main__":
    main()
