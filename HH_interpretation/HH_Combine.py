# HH_Combine.py

import sympy
from collections import OrderedDict # Or just use a list of dicts

# --- Input Data: ggF HH samples ---
# Based on your provided data
# Each sample has 'kl', 'kt', 'xs' (cross-section), and a 'label'
ggf_samples_data = [
    {"kl": 0.0, "kt": 1.0, "xs": 0.069725, "label": "ggHH_kl_0_kt_1"},
    {"kl": 1.0, "kt": 1.0, "xs": 0.031047, "label": "ggHH_kl_1_kt_1"},
    {"kl": 2.45, "kt": 1.0, "xs": 0.013124, "label": "ggHH_kl_2p45_kt_1"},
    {"kl": 5.0, "kt": 1.0, "xs": 0.091172, "label": "ggHH_kl_5_kt_1"},
]

vbf_samples_data = [
    
    {"CV": 1.0,     "C2V": 0.0,    "kl": 1.0,   "xs": 0.029320, "label": "CV_1_C2V_0_C3_1"},
    {"CV": 1.0,     "C2V": 1.0,    "kl": 1.0,   "xs": 0.001906, "label": "CV_1_C2V_1_C3_1"},
    {"CV": 1.74,    "C2V": 1.37,   "kl": 14.4,  "xs": 0.395400, "label": "CV-1p74_C2V-1p37_C3-14p4"},
    {"CV": -0.012,  "C2V": 0.03,   "kl": 10.2,  "xs": 0.00001257, "label": "CV-m0p012_C2V-0p03_C3-10p2"},
    {"CV": -0.758,  "C2V": 1.44,   "kl": -19.3,  "xs": 0.35500, "label": "CV-m0p758_C2V-1p44_C3-m19p3"},
    {"CV": -0.962,  "C2V": 0.959,  "kl": -1.43, "xs": 0.00111, "label": "CV-m0p962_C2V-0p959_C3-m1p43"},
    {"CV": -1.21,   "C2V": 1.94,  "kl": -0.94, "xs": 0.0033739, "label": "CV-m1p21_C2V-1p94_C3-m0p94"},
    {"CV": -1.60,   "C2V": 2.72,  "kl": -1.36, "xs": 0.0105109, "label": "CV-m1p60_C2V-2p72_C3-m1p36"},
    {"CV": -1.83,   "C2V": 3.57,  "kl": -3.39, "xs": 0.0149850, "label": "CV-m1p83_C2V-3p57_C3-m3p39"},
    {"CV": -2.12,   "C2V": 3.87,  "kl": -5.96, "xs": 0.6322811, "label": "CV-m2p12_C2V-3p87_C3-m5p96"},
    # above are what were generated in Run3 (2022)
    
    {"CV": 1.0,     "C2V": 0.0,    "kl": 2.9,   "xs": 0.0191330, "label": "CV_0p4_C2V_1_C3_2p9"},
    {"CV": 1.1,     "C2V": 1.0,    "kl": 0.2,   "xs": 0.0095458, "label": "CV_1p1_C2V_1_C3_0p2"},
    {"CV": 1.0, "C2V": 1.0, "kl": 0.0, "xs": 0.0046089, "label": "qqHH_CV_1_C2V_1_kl_0"},
    {"CV": 1.0, "C2V": 1.0, "kl": 2.0, "xs": 0.0014228, "label": "qqHH_CV_1_C2V_1_kl_2"},
    {"CV": 1.0, "C2V": 2.0, "kl": 1.0, "xs": 0.0142178, "label": "qqHH_CV_1_C2V_2_kl_1"},
    {"CV": 0.5, "C2V": 1.0, "kl": 1.0, "xs": 0.0108237, "label": "qqHH_CV_0p5_C2V_1_kl_1"},
    {"CV": 1.5, "C2V": 1.0, "kl": 1.0, "xs": 0.0660185, "label": "qqHH_CV_1p5_C2V_1_kl_1"},
    
]

# We need at least 3 samples for the 3-operator basis (box, triangle, interf)
MIN_SAMPLES = 3
if len(ggf_samples_data) < MIN_SAMPLES:
    raise ValueError(f"Need at least {MIN_SAMPLES} samples for the calculation.")


def calculate_ggf_parameters(target_kl, samples_data, target_kt=1.0):
    """
    Calculates the coefficients to reweight ggF HH samples for a target kl
    (assuming kt=1.0) and the predicted cross-section.

    The formula for ggF HH cross-section parametrization is often given as:
    xs(kl, kt) = A * kt^4 + B * kl^2 * kt^2 + C * kl * kt^3
    where A, B, C are coefficients related to box, triangle (top-loop, Higgs-tria),
    and interference terms respectively.

    This function determines weights for existing samples to match a target (kl, kt).

    Args:
        target_kl (float): The target kl coupling modifier.
        samples_data (list): A list of dictionaries, where each dictionary
                             contains 'kl', 'kt', and 'xs' for a sample.
        target_kt (float, optional): The target kt coupling modifier.
                                     Defaults to 1.0 as per user request.

    Returns:
        tuple: A tuple containing:
            -coeffs_values (list): A list of numerical coefficients for each input sample.
            -predicted_xs_val (float): The predicted cross-section for the target (kl, kt).
            -coeffs_expr (sympy.Matrix): The symbolic expression for the coefficients.
            -predicted_xs_expr (sympy.Expr): The symbolic expression for the cross-section.
    """
    n_samples = len(samples_data)

    # --- 1. Construct the matrix M from the existing samples ---
    # Each row in M corresponds to a sample and its coupling dependence.
    # The columns correspond to the basis terms: kt^4, kt^2*kl^2, kt^3*kl
    M_list = []
    for sample in samples_data:
        s_kt = sample['kt']
        s_kl = sample['kl']
        M_list.append([
            s_kt**4,
            s_kt**2 * s_kl**2,
            s_kt**3 * s_kl
        ])
    M = sympy.Matrix(M_list)

    # --- 2. Define the target coupling vector c ---
    # This vector represents the coupling dependence for the target (kl, kt).
    kl_sym, kt_sym = sympy.symbols("kl kt") # Symbolic kl, kt for the target
    c_expr = sympy.Matrix([
        [kt_sym**4],
        [kt_sym**2 * kl_sym**2],
        [kt_sym**3 * kl_sym]
    ])

    # Substitute the actual target kl and kt values
    c_target = c_expr.subs({kl_sym: target_kl, kt_sym: target_kt})

    # --- 3. Define the vector of symbolic sample cross sections s ---
    # This is used to calculate the combined cross-section later.
    s_list_sym = [sympy.Symbol(f"xs{i}") for i in range(n_samples)]
    s_sym = sympy.Matrix(s_list_sym)

    # Vector of actual cross-sections from input data
    s_values = sympy.Matrix([sample['xs'] for sample in samples_data])


    # --- 4. Calculate coefficients ---
    # coeffs = c^T * M_pseudo_inverse
    # M_pinv is the Moore-Penrose pseudo-inverse, useful for non-square matrices
    # or singular matrices, providing a least-squares solution.
    try:
        M_pinv = M.pinv()
    except Exception as e:
        print(f"Error calculating pseudo-inverse of M: {e}")
        print("Matrix M:")
        sympy.pprint(M)
        # Check for linear dependency if a true inverse is expected (for 3x3 M)
        if M.shape == (3,3):
            print(f"Determinant of M: {M.det()}")
        return None, None, None, None

    # Symbolic coefficients (coeffs_expr is a 1 x n_samples matrix)
    coeffs_expr_transposed = c_expr.transpose() * M_pinv # (1x3) * (3xN_samples) = 1xN_samples
    coeffs_expr = coeffs_expr_transposed.transpose() # N_samples x 1, for consistency if preferred

    # Numerical coefficients for the target kl, kt
    coeffs_target_transposed = c_target.transpose() * M_pinv
    # Convert sympy.Float to Python float
    coeffs_target_values = [float(val.evalf()) for val in coeffs_target_transposed]


    # --- 5. Calculate the predicted cross section ---
    # sigma = coeffs * s
    predicted_xs_expr = (coeffs_expr_transposed * s_sym)[0,0] # (1xN_samples) * (N_samplesx1) -> 1x1
    
    # Numerical predicted cross-section
    predicted_xs_val_matrix = coeffs_target_transposed * s_values
    predicted_xs_val = predicted_xs_val_matrix[0,0].evalf() if predicted_xs_val_matrix else None

    return coeffs_target_values, predicted_xs_val, coeffs_expr_transposed, predicted_xs_expr

def calculate_vbf_parameters(target_CV, target_C2V, target_kl, samples_data):
    """
    Calculates weights for qqHH (VBF) samples to match a target (CV, C2V, kl).
    Uses 6-term basis from VBFFormula.

    Args:
        target_CV (float): Target CV coupling.
        target_C2V (float): Target C2V coupling.
        target_kl (float): Target kl coupling.
        samples_data (list): List of dicts with 'CV', 'C2V', 'kl', 'xs', and 'label'.

    Returns:
        tuple: (weights, predicted_xs, symbolic_weights, symbolic_sigma)
    """
    n_samples = len(samples_data)
    if n_samples < 6:
        raise ValueError("Need at least 6 samples for VBF reweighting.")

    # --- 1. Build the sample matrix M ---
    M_list = []
    for s in samples_data:
        M_list.append([
            s["CV"]**2 * s["kl"]**2,
            s["CV"]**4,
            s["C2V"]**2,
            s["CV"]**3 * s["kl"],
            s["CV"] * s["C2V"] * s["kl"],
            s["CV"]**2 * s["C2V"],
        ])
    M = sympy.Matrix(M_list)

    # --- 2. Build target vector c ---
    CV_sym, C2V_sym, kl_sym = sympy.symbols("CV C2V kl")
    c_expr = sympy.Matrix([
        [CV_sym**2 * kl_sym**2],
        [CV_sym**4],
        [C2V_sym**2],
        [CV_sym**3 * kl_sym],
        [CV_sym * C2V_sym * kl_sym],
        [CV_sym**2 * C2V_sym],
    ])
    c_target = c_expr.subs({CV_sym: target_CV, C2V_sym: target_C2V, kl_sym: target_kl})

    # --- 3. Symbolic ¦Ò and weights ---
    s_syms = [sympy.Symbol(f"xs{i}") for i in range(n_samples)]
    s_vec = sympy.Matrix(s_syms)
    xs_values = sympy.Matrix([s["xs"] for s in samples_data])

    try:
        M_pinv = M.pinv()
    except Exception as e:
        print("Error in matrix inversion:", e)
        return None, None, None, None

    coeffs_expr_T = c_expr.transpose() * M_pinv
    coeffs_expr = coeffs_expr_T.transpose()

    coeffs_vals_T = c_target.transpose() * M_pinv
    coeffs_vals = [float(val.evalf()) for val in coeffs_vals_T]

    sigma_expr = (coeffs_expr_T * s_vec)[0, 0]
    sigma_val = float((coeffs_vals_T * xs_values)[0, 0].evalf())

    return coeffs_vals, sigma_val, coeffs_expr, sigma_expr

if __name__ == "__main__":
    print("=== HH Combine Reweighting Tool ===")
    mode = input("Choose production mode [ggf/qqhh]: ").strip().lower()

    if mode == "ggf":
        # Use GGF benchmarks
        print("\nAvailable ggF benchmark samples:")
        for i, sample in enumerate(ggf_samples_data):
            print(f"  Index {i}: kl={sample['kl']}, kt={sample['kt']}, xs={sample['xs']:.5f} ({sample['label']})")

        while True:
            try:
                selected_indices_str = input(f"\nEnter the indices of the samples to use (comma-separated, need at least {MIN_SAMPLES}): ")
                selected_indices = [int(idx.strip()) for idx in selected_indices_str.split(',')]
                if len(selected_indices) < MIN_SAMPLES:
                    print(f"Please select at least {MIN_SAMPLES} samples.")
                    continue

                active_samples = [ggf_samples_data[idx] for idx in selected_indices]
                break
            except Exception as e:
                print(f"Error: {e}")

        target_kl = float(input("Enter the target kl: "))
        target_kt = float(input("Enter the target kt [default = 1.0]: ") or 1.0)

        print(f"\nCalculating for target kl = {target_kl}, kt = {target_kt}")
        coefficients, predicted_xs, _, _ = calculate_ggf_parameters(
            target_kl=target_kl,
            samples_data=active_samples,
            target_kt=target_kt
        )

    elif mode == "qqhh":
        print("\nAvailable qqHH (VBF) benchmark samples:")
        for i, sample in enumerate(vbf_samples_data):
            print(f"  Index {i}: CV={sample['CV']}, C2V={sample['C2V']}, kl={sample['kl']}, xs={sample['xs']:.5f} ({sample['label']})")

        while True:
            try:
                selected_indices_str = input(f"\nEnter the indices of the VBF samples to use (at least 6): ")
                selected_indices = [int(idx.strip()) for idx in selected_indices_str.split(',')]
                if len(selected_indices) < 6:
                    print("Please select at least 6 samples for VBF.")
                    continue

                active_samples = [vbf_samples_data[idx] for idx in selected_indices]
                break
            except Exception as e:
                print(f"Error: {e}")

        target_CV = float(input("Enter the target CV: "))
        target_C2V = float(input("Enter the target C2V: "))
        target_kl = float(input("Enter the target kl: "))

        print(f"\nCalculating for target CV = {target_CV}, C2V = {target_C2V}, kl = {target_kl}")
        coefficients, predicted_xs, _, _ = calculate_vbf_parameters(
            target_CV=target_CV,
            target_C2V=target_C2V,
            target_kl=target_kl,
            samples_data=active_samples
        )

    else:
        print("Invalid mode. Please choose either 'ggf' or 'qqhh'.")
        exit()

    # --- Display Results ---
    if coefficients is not None and predicted_xs is not None:
        print("\n--- Results ---")
        for i, coeff in enumerate(coefficients):
            print(f"  Weight for {active_samples[i]['label']} (xs{i}): {coeff:.6f}")
        print(f"\nPredicted cross-section: {predicted_xs:.6f} pb")
    else:
        print("Calculation failed.")


