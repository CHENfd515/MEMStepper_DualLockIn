import numpy as np
import pandas as pd
import re
import io
from datetime import datetime
import matplotlib.pyplot as plt

# =============================================================================
# DATA LOADING: Parses the analy_report.txt file using Regex
# =============================================================================
def load_report(file_path):
    rows = []
    
    # Identify dynamic column names from the header (e.g., Mean_l100)
    with open(file_path, 'r', encoding='utf-8') as f:
        header_line = ""
        for line in f:
            if "Mean_l" in line:
                header_line = line.strip()
                break

        last_mean_col = None
        last_mae_col = None
        last_mse_col = None

        if header_line:
            cols = [c.strip() for c in header_line.split("|")]
            for c in cols:
                if c.startswith("Mean_l"): last_mean_col = c
                if c.startswith("MAE_l"):  last_mae_col = c
                if c.startswith("MSE_l"):  last_mse_col = c

    # Regex to capture Group ID, PID Gains, Global Stats, Window Stats, and Dominant factor
    pattern = re.compile(
        r"G(\d+)\s+\|\s+"
        r"([\deE\+\-\.]+)\s+\|\s+([\deE\+\-\.]+)\s+\|\s+([\deE\+\-\.]+)\s+\|\s+"  # Kp, Ki, Kd
        r"([\deE\+\-\.]+)\s+\|\s+([\deE\+\-\.]+)\s+\|\s+([\deE\+\-\.]+)\s+\|\s+"  # Mean, Std, MAE
        r"([\deE\+\-\.]+)\s+\|\s+([\deE\+\-\.]+)\s+\|\s+([\deE\+\-\.]+)\s+\|\s+"  # Window Mean, MAE, MSE
        r"([\d\.]+)\s+\|\s+([\d\.]+)\s+\|\s+([\d\.]+)\s+\|\s+([PID])"           # P/I/D %, Dominant
    )

    with open(file_path, 'r', encoding='utf-8') as f:
        for line in f:
            m = pattern.search(line)
            if m:
                rows.append({
                    'Group': int(m.group(1)),
                    'Kp': float(m.group(2)),
                    'Ki': float(m.group(3)),
                    'Kd': float(m.group(4)),
                    'Mean': float(m.group(5)),
                    'Std': float(m.group(6)),
                    'MAE': float(m.group(7)),
                    last_mean_col: float(m.group(8)) if last_mean_col else None,
                    last_mae_col: float(m.group(9)) if last_mae_col else None,
                    last_mse_col: float(m.group(10)) if last_mse_col else None,
                    'P_pct': float(m.group(11)),
                    'I_pct': float(m.group(12)),
                    'D_pct': float(m.group(13)),
                    'Dominant': m.group(14)
                })

    return pd.DataFrame(rows)

# =============================================================================
# PARETO CALCULATION: Finds non-dominated solutions (minimizing all objectives)
# =============================================================================
def pareto_front(df, objectives):
    """
    Identifies the Pareto optimal set.
    objectives: list of column names to MINIMIZE.
    """
    values = df[objectives].values
    is_pareto = np.ones(values.shape[0], dtype=bool)

    for i in range(values.shape[0]):
        for j in range(values.shape[0]):
            # Check if point j dominates point i
            if all(values[j] <= values[i]) and any(values[j] < values[i]):
                is_pareto[i] = False
                break

    return df.iloc[is_pareto]

# =============================================================================
# VISUALIZATION: Plots candidates and highlights the Pareto Front
# =============================================================================
def plot_pareto(df, pareto_df, opt_last=False):
    current_time = datetime.now().strftime("%m%d")
    save_path = f"pareto_front_{current_time}"

    last_mean_col = [c for c in df.columns if c.startswith("abs_Mean_l")]
    last_mae_col  = [c for c in df.columns if c.startswith("MAE_l")]
    last_mse_col  = [c for c in df.columns if c.startswith("MSE_l")]

    if opt_last and last_mean_col:
        x_col = last_mean_col[0]  # Absolute Mean (X-axis)
        y_col = last_mae_col[0]   # MAE (Y-axis)
        c_col = last_mse_col[0]   # MSE (Color)
        title = "PID Pareto Optimization (|Mean| vs MAE)"
    else:
        x_col = "Std"
        y_col = "MAE"
        c_col = "Mean"
        title = "Global PID Pareto Front (Std vs MAE)"

    plt.figure(figsize=(10, 8))

    # Scatter plot for all candidate configurations
    sc_all = plt.scatter(
        df[x_col], df[y_col], c=df[c_col],
        cmap='viridis', s=70, edgecolors='black', alpha=0.6
    )

    # --- Logic to avoid overlap ---
    # We use a staggered offset to help separate labels
    for i, (_, row) in enumerate(df.iterrows()):
        offset_x = (i % 3 - 1) * 0.00001  # Slight staggered x offset
        offset_y = (i % 2 * 2 - 1) * 0.0001 # Slight staggered y offset
        
        label = f"G{row['Group']}\nP={row['Kp']:.2f}\nI={row['Ki']:.4f}\nD={row['Kd']:.2f}"
        plt.text(
            row[x_col] + offset_x, 
            row[y_col] + offset_y, 
            label, 
            fontsize=2, 
            fontweight='bold',
            bbox=dict(facecolor='white', alpha=0.5, edgecolor='none', pad=0.5)
        )
        
    # Highlight Pareto Optimal points with red circles
    plt.scatter(
        pareto_df[x_col], pareto_df[y_col],
        s=200, facecolors='none', edgecolors='red', linewidths=2, label='Pareto Optimal'
    )

    # Plot Pareto Front line (Sorted by X-axis)
    pareto_sorted = pareto_df.sort_values(by=x_col)
    plt.plot(pareto_sorted[x_col], pareto_sorted[y_col], 'r--', linewidth=2, label='Frontier Line')

    plt.colorbar(sc_all).set_label(c_col)
    plt.xlabel(f"Objective 1: {x_col}")
    plt.ylabel(f"Objective 2: {y_col}")
    plt.title(title)
    plt.grid(True, linestyle='--', alpha=0.4)
    plt.legend()
    plt.tight_layout()
    plt.savefig(save_path + ".png", dpi=300)
    plt.savefig(save_path + ".svg", dpi=300, format='svg', bbox_inches='tight')
    print(f"Plot saved: {save_path}.png and {save_path}.svg")

# =============================================================================
# MAIN OPTIMIZER: Candidate selection and Pareto ranking
# =============================================================================
def optimize_pid(report_path, opt_last=True, selection_mode="quantile", thresholds=None, top_k=10):
    df = load_report(report_path)

    last_mean_col = [c for c in df.columns if c and c.startswith("Mean_l")]
    last_mae_col = [c for c in df.columns if c and c.startswith("MAE_l")]
    last_mse_col = [c for c in df.columns if c and c.startswith("MSE_l")]

    if opt_last and last_mean_col:
        m_col = last_mean_col[0]
        mae_col = last_mae_col[0]
        mse_col = last_mse_col[0]
        
        # --- KEY MODIFICATION: Create Absolute Mean Column for Optimization ---
        abs_mean_col = 'abs_' + m_col
        df[abs_mean_col] = df[m_col].abs()

        print(f"\n[INFO] Optimizing based on windowed stats. Objective: Minimize {abs_mean_col}")

        # --- Candidate Filtering ---
        if selection_mode == "quantile":
            mae_q  = thresholds.get("mae_q", 0.5) if thresholds else 0.5
            mean_q = thresholds.get("mean_q", 0.5) if thresholds else 0.5
            
            mae_th  = df[mae_col].quantile(mae_q)
            mean_th = df[abs_mean_col].quantile(mean_q)
            
            candidates = df[(df[mae_col] <= mae_th) & (df[abs_mean_col] <= mean_th)].copy()

        elif selection_mode == "topk":
            # Rank based on closeness to 0 (abs mean) and low MAE
            score = df[mae_col].rank() + df[abs_mean_col].rank()
            candidates = df.iloc[score.argsort()].head(top_k).copy()

        # Define objectives to minimize for Pareto Front
        objectives = [abs_mean_col, mae_col, mse_col]

    else:
        # Global stats optimization
        df['abs_Mean'] = df['Mean'].abs()
        candidates = df[df['MAE'] <= df['MAE'].quantile(0.4)].copy()
        objectives = ['abs_Mean', 'Std']

    # Calculate Pareto Optimal set
    pareto = pareto_front(candidates, objectives)
    
    # Generate visualization
    plot_pareto(candidates, pareto, opt_last=opt_last)

    # --- Sorting the Pareto results by Absolute Mean (Ascending) ---
    # Sort Pareto results by Absolute Mean (closest to 0 first)
    pareto_sorted = pareto.sort_values(by=abs_mean_col, ascending=True)

    print("\n" + "="*85)
    print("FINAL PARETO OPTIMAL PID CONFIGURATIONS (Sorted by |Mean|)")
    print("="*85)
    # Header with MSE included
    print(f"{'GROUP':<6} {'Kp':<7} {'Ki':<8} {'Kd':<7} {'|Mean|':<10} {'MAE':<10} {'MSE':<12}")
    print("-" * 85)
    
    for _, row in pareto_sorted.iterrows():
        print(f"{int(row['Group']):<6} "
              f"{row['Kp']:<7.3f} "
              f"{row['Ki']:<8.4f} "
              f"{row['Kd']:<7.3f} "
              f"{row[abs_mean_col]:<10.5f} "
              f"{row[mae_col]:<10.4f} "
              f"{row[mse_col]:<5.5f}")
    print("="*85)
    
    return pareto_sorted

if __name__ == "__main__":
    # Settings for Optimization
    # thresholds["mean_q"] controls how strictly we demand Mean to be close to 0
    # thresholds["mae_q"] controls how strictly we demand low overall error
    optimize_pid(
        "analy_report.txt",
        opt_last=True,
        selection_mode="quantile",
        thresholds={"mae_q": 0.5, "mean_q": 0.2, "mse_q": 0.5}
    )