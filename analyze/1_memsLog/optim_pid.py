import numpy as np
import pandas as pd
import re
from datetime import datetime
import matplotlib.pyplot as plt



# =========================
# read analy_report.txt
# =========================
def load_report(file_path):
    rows = []

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
                if c.startswith("Mean_l"):
                    last_mean_col = c
                if c.startswith("MAE_l"):
                    last_mae_col = c
                if c.startswith("MSE_l"):
                    last_mse_col = c

    pattern = re.compile(
        r"G(\d+)\s+\|\s+"
        r"([\deE\+\-\.]+)\s+\|\s+([\deE\+\-\.]+)\s+\|\s+([\deE\+\-\.]+)\s+\|\s+"
        r"([\deE\+\-\.]+)\s+\|\s+([\deE\+\-\.]+)\s+\|\s+([\deE\+\-\.]+)\s+\|\s+"
        r"([\deE\+\-\.]+)\s+\|\s+([\deE\+\-\.]+)\s+\|\s+([\deE\+\-\.]+)\s+\|\s+"
        r"([\d\.]+)\s+\|\s+([\d\.]+)\s+\|\s+([\d\.]+)\s+\|\s+([PID])"
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



# =========================
# Pareto Front Visualization
# =========================
def plot_pareto(df, pareto_df, opt_last=False):

    current_time = datetime.now().strftime("%m%d")
    save_path = f"pareto_front_{current_time}.png"

    last_mean_col = [c for c in df.columns if c.startswith("Mean_l")]
    last_mae_col  = [c for c in df.columns if c.startswith("MAE_l")]
    last_mse_col  = [c for c in df.columns if c.startswith("MSE_l")]

    if opt_last and last_mean_col:
        x_col = last_mean_col[0]
        y_col = last_mae_col[0]
        c_col = last_mse_col[0]
        title = "Pareto Front (Last Window Optimization)"
    else:
        x_col = "Std"
        y_col = "MAE"
        c_col = "P2P"
        title = "Pareto Front (Global Optimization)"

    plt.figure(figsize=(10, 8))

    # =========================
    sc_all = plt.scatter(
        df[x_col],
        df[y_col],
        c=df[c_col],
        cmap='viridis',
        s=70,
        edgecolors='black'
    )

    for _, row in df.iterrows():
        label = f"G{row['Group']}\nP={row['Kp']:.1e}\nI={row['Ki']:.1e}\nD={row['Kd']:.1e}"
        plt.text(
            row[x_col],
            row[y_col],
            label,
            fontsize=9,
            ha='left',
            va='bottom'
        )

    # =========================
    plt.scatter(
        pareto_df[x_col],
        pareto_df[y_col],
        s=180,
        facecolors='none',
        edgecolors='red',
        linewidths=2,
        label='Pareto Optimal'
    )

    # =========================
    pareto_sorted = pareto_df.sort_values(by=x_col)
    x_p = pareto_sorted[x_col].values
    y_p = pareto_sorted[y_col].values

    if len(x_p) >= 3:
        coeffs = np.polyfit(x_p, y_p, 2)
        poly = np.poly1d(coeffs)

        x_fit = np.linspace(min(x_p), max(x_p), 200)
        y_fit = poly(x_fit)

        plt.plot(
            x_fit,
            y_fit,
            linestyle='--',
            linewidth=2,
            label='Pareto Front Fit'
        )
    else:
        plt.plot(x_p, y_p, linestyle='--', linewidth=2, label='Pareto Front')

    # =========================
    # Colorbar
    # =========================
    cbar = plt.colorbar(sc_all)
    cbar.set_label(c_col)

    plt.xlabel(x_col)
    plt.ylabel(y_col)
    plt.title(title)
    plt.grid(True, linestyle='--', alpha=0.4)
    plt.legend()
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()

    print(f"Pareto front plot saved as {save_path}")
      
def pareto_front(df, objectives):
    """
    df: pandas DataFrame
    objectives: list of column names to MINIMIZE
    """
    values = df[objectives].values
    is_pareto = np.ones(values.shape[0], dtype=bool)

    for i in range(values.shape[0]):
        for j in range(values.shape[0]):
            if all(values[j] <= values[i]) and any(values[j] < values[i]):
                is_pareto[i] = False
                break

    return df.iloc[is_pareto]

# =========================
# optimize PID parameters based on analy_report.txt
# =========================
def optimize_pid(report_path,
                 opt_last=False,
                 selection_mode="quantile",
                 thresholds=None,
                 top_k=10):
    df = load_report(report_path)

    last_mean_col = [c for c in df.columns if c and c.startswith("Mean_l")]
    last_mae_col = [c for c in df.columns if c and c.startswith("MAE_l")]
    last_mse_col = [c for c in df.columns if c and c.startswith("MSE_l")]

    if opt_last and last_mean_col:
        mean_col = last_mean_col[0]
        mae_col = last_mae_col[0]
        mse_col = last_mse_col[0]

        print(f"\nUsing LAST window optimization based on {mean_col}")

        # =========================
        # Candidate Selection
        # =========================
        if selection_mode == "quantile":
            mae_q  = thresholds.get("mae_q", 0.4)  if thresholds else 0.4
            mean_q = thresholds.get("mean_q", 0.4) if thresholds else 0.4
            mse_q  = thresholds.get("mse_q", 0.6)  if thresholds else 0.6

            mae_th  = df[mae_col].quantile(mae_q)
            mean_th = df[mean_col].abs().quantile(mean_q)
            mse_th  = df[mse_col].quantile(mse_q)

            candidates = df[
                (df[mae_col] <= mae_th) &
                (df[mean_col].abs() <= mean_th) &
                (df[mse_col] <= mse_th)
            ]


        elif selection_mode == "fixed":
            candidates = df[
                (df[mae_col] <= thresholds["mae"]) &
                (df[mean_col].abs() <= thresholds["mean"]) &
                (df[mse_col] <= thresholds["mse"])
            ]


        elif selection_mode == "topk":
            score = (
                df[mae_col].rank() +
                df[mean_col].abs().rank() +
                df[mse_col].rank()
            )

            df_sorted = df.iloc[score.argsort()]
            candidates = df_sorted.head(top_k)

        else:
            raise ValueError("Invalid selection_mode")


        objectives = [mae_col, mean_col, mse_col]

    else:
        print("\nUsing GLOBAL optimization")

        mae_th = df['MAE'].quantile(0.4)
        std_th = df['Std'].quantile(0.6)

        candidates = df[
            (df['MAE'] <= mae_th) &
            (df['Std'] <= std_th)
        ]

        objectives = ['MAE', 'Std']

    print("\n========== Valuable Candidates ==========")
    print(candidates[['Group', 'Kp', 'Ki', 'Kd'] + (objectives if opt_last else [])])

    pareto = pareto_front(candidates, objectives)
    plot_pareto(candidates, pareto, opt_last=opt_last)


    print("\n========== Pareto Front ==========")
    print(pareto[['Group', 'Kp', 'Ki', 'Kd']+ (objectives if opt_last else [])])

    return pareto


if __name__ == "__main__":
    # optimize_pid("analy_report.txt", opt_last=True)
    # optimize_pid("analy_report.txt",
    #             opt_last=True,
    #             selection_mode="topk",
    #             top_k=8)
    optimize_pid("analy_report.txt",
             opt_last=True,
             selection_mode="quantile",
             thresholds={"mae_q":0.9, "mean_q":0.3, "mse_q":0.9})
    # optimize_pid("analy_report.txt",
    #          opt_last=True,
    #          selection_mode="fixed",
    #          thresholds={"mae":1e-3,
    #                      "mean":5e-4,
    #                      "mse":1e-6})

