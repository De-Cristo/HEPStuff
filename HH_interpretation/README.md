# 🧪 HH Signal Combination & Benchmark Plotting Toolkit

This toolkit contains two Python scripts for visualizing and synthesizing Higgs boson pair production (`ggF → HH`) signal shapes at arbitrary values of the Higgs self-coupling modifier (κλ), using benchmark Monte Carlo samples and a 3-operator reweighting model.

---

## 📜 Included Scripts

### 1. `HH_Combine.py`

**Purpose**:  
Provides the function `calculate_ggf_parameters()` that computes weights to linearly combine benchmark samples to represent any target κλ (with fixed κt = 1.0).

**Features**:
- Uses symbolic algebra via `sympy` to solve the parameterization:
  \[
  \sigma(\kappa_\lambda, \kappa_t) = A \kappa_t^4 + B \kappa_\lambda^2 \kappa_t^2 + C \kappa_\lambda \kappa_t^3
  \]
- Computes coefficients (`weights`) and predicted cross-section.
- Has an interactive CLI to test combinations manually.

---

### 2. `combine_HH_histograms.py`

**Purpose**:  
Plots both the original benchmark histograms and a combined synthetic signal shape for a user-specified `κλ` value.

**Workflow**:
1. **Plot individual benchmark samples** (`kl = 0.0`, `2.45`, `5.0`) using ROOT data.
2. **Prompt user for target κλ** (e.g. 1.0).
3. **Use `HH_Combine.calculate_ggf_parameters()`** to get weights.
4. **Apply weights to histograms** to synthesize the target shape.
5. **Plot and save**:
   - `benchmark_histograms.png` for benchmarks.
   - `combined_hist_kl_<target>.png` for combined shape.

**Run it**:
```bash
python combine_HH_histograms.py
```