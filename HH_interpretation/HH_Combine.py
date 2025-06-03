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


if __name__ == "__main__":
    print("Available benchmark samples:")
    for i, sample in enumerate(ggf_samples_data):
        print(f"  Index {i}: kl={sample['kl']}, kt={sample['kt']}, xs={sample['xs']:.5f} ({sample['label']})")

    while True:
        try:
            selected_indices_str = input(f"\nEnter the indices of the samples to use (comma-separated, need at least {MIN_SAMPLES}): ")
            selected_indices = [int(idx.strip()) for idx in selected_indices_str.split(',')]

            if len(selected_indices) < MIN_SAMPLES:
                print(f"Please select at least {MIN_SAMPLES} samples.")
                continue

            active_samples = []
            valid_selection = True
            for idx in selected_indices:
                if 0 <= idx < len(ggf_samples_data):
                    active_samples.append(ggf_samples_data[idx])
                else:
                    print(f"Error: Index {idx} is out of range.")
                    valid_selection = False
                    break
            if not valid_selection:
                continue
            
            # Optional: Check if all selected samples have kt=1 if that's an assumption
            # for a particular interpretation of "fixed kt" mode.
            # For now, the calculate_ggf_parameters function handles kt as given.

            break # Exit loop if selection is valid

        except ValueError:
            print("Invalid input. Please enter comma-separated numbers (e.g., 0, 1, 2).")
        except Exception as e:
            print(f"An error occurred: {e}")

    # --- User Input for target kl ---
    try:
        input_kl_str = input("Enter the target kl value (e.g., 1.0, 2.0): ")
        target_kl_input = float(input_kl_str)
    except ValueError:
        print("Invalid input. Please enter a numeric value for kl.")
        exit()

    # Assuming target_kt is 1.0 as per previous discussions for this scenario
    target_kt_input = 1.0

    print(f"\nCalculating for target kl = {target_kl_input} (with target kt = {target_kt_input})")
    print("Using the following selected benchmark samples:")
    for i, sample in enumerate(active_samples): # Changed to active_samples
        print(f"  Sample {i} (Original Index {ggf_samples_data.index(sample)}): kl={sample['kl']}, kt={sample['kt']}, xs={sample['xs']:.5f} ({sample['label']})")
    print("-" * 30)

    # --- Perform Calculation ---
    coefficients, predicted_xs, sym_coeffs, sym_xs_formula = calculate_ggf_parameters(
        target_kl=target_kl_input,
        samples_data=active_samples, # Pass only the selected samples
        target_kt=target_kt_input
    )

    # --- Display Results ---
    if coefficients is not None and predicted_xs is not None:
        print("\n--- Results ---")
        print(f"Target: kl = {target_kl_input}, kt = {target_kt_input}")
        
        print("\nCalculated Coefficients (Weights) for each selected input sample:")
        for i, coeff in enumerate(coefficients):
            sample_label = active_samples[i]['label'] # Use active_samples here
            print(f"  Weight for {sample_label} (xs{i}): {coeff:.6f}")

        print(f"\nPredicted Cross-Section (sigma) for target kl={target_kl_input}, kt={target_kt_input}: {predicted_xs:.6f}")
    else:
        print("Calculation failed. Please check input data or matrix properties (e.g., selected samples might lead to a singular matrix if not distinct enough).")

