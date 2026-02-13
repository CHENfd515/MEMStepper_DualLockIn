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

    pattern = re.compile(
        r"G(\d+)\s+\|\s+"
        r"([\deE\+\-\.]+)\s+\|\s+([\deE\+\-\.]+)\s+\|\s+([\deE\+\-\.]+)\s+\|\s+"
        r"([\deE\+\-\.]+)\s+\|\s+([\deE\+\-\.]+)\s+\|\s+([\deE\+\-\.]+)\s+\|\s+([\deE\+\-\.]+)\s+\|\s+"
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
                    'P2P': float(m.group(7)),
                    'MAE': float(m.group(8)),
                    'P_pct': float(m.group(9)),
                    'I_pct': float(m.group(10)),
                    'D_pct': float(m.group(11)),
                    'Dominant': m.group(12)
                })

    return pd.DataFrame(rows)


# =========================
# Pareto Front Visualization
# =========================
def plot_pareto(df, pareto_df):
    current_time = datetime.now().strftime("%m%d")
    save_path = f"pareto_front_{current_time}.png"
    """
    df         : all candidate PID sets
    pareto_df  : Pareto-optimal subset
    """

    plt.figure(figsize=(9, 7))

    # candidates
    plt.scatter(
        df['Std'],
        df['MAE'],
        c='lightgray',
        s=40,
        label='Candidates'
    )

    # pareto optimal
    sc = plt.scatter(
        pareto_df['Std'],
        pareto_df['MAE'],
        c=pareto_df['P2P'],
        cmap='viridis',
        s=120,
        edgecolors='black',
        label='Pareto Optimal'
    )

    for _, row in pareto_df.iterrows():
        plt.text(
            row['Std'],
            row['MAE'],
            row['Group'],
            fontsize=9,
            ha='left',
            va='bottom'
        )

    cbar = plt.colorbar(sc)
    cbar.set_label('P2P (Peak-to-Peak)')

    plt.xlabel('Std (Stability)')
    plt.ylabel('MAE (Accuracy)')
    plt.title('Pareto Front of PID Parameter Sets')
    plt.grid(True, linestyle='--', alpha=0.4)
    plt.legend()
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
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
def optimize_pid(report_path):
    df = load_report(report_path)

    print("\n========== Loaded PID Report ==========")
    print(df[['Group', 'Kp', 'Ki', 'Kd', 'MAE', 'Std', 'P2P', 'Dominant']])

    # =========================
    # 1-Valuable Candidates Selection
    # =========================
    mae_th = df['MAE'].quantile(0.4)
    std_th = df['Std'].quantile(0.6)
    p2p_th = df['P2P'].quantile(0.6)

    candidates = df[
        (df['MAE'] <= mae_th) &
        (df['Std'] <= std_th) &
        (df['P2P'] <= p2p_th)
    ]

    print("\n========== Valuable Candidates ==========")
    print(candidates[['Group', 'Kp', 'Ki', 'Kd', 'MAE', 'Std', 'P2P']])

    # =========================
    # 2-Pareto Front Identification
    # =========================
    pareto = pareto_front(
        candidates,
        objectives=['MAE', 'Std', 'P2P']
    )
    plot_pareto(candidates, pareto)


    print("\n========== Pareto Front ==========")
    print(pareto[['Group', 'Kp', 'Ki', 'Kd', 'MAE', 'Std', 'P2P']])

    # =========================
    # 3-Recommended Search Space Expansion
    # =========================
    def expand(values, ratio=0.3):
        center = np.median(values)
        return sorted(set([
            round(center * (1 - ratio), 6),
            round(center, 6),
            round(center * (1 + ratio), 6)
        ]))

    Kp_list = expand(pareto['Kp'].values)
    Ki_list = expand(pareto['Ki'].values)
    Kd_list = expand(pareto['Kd'].values)

    print("\n========== Recommended PID Search Space ==========")
    print(f"Kp_list = {Kp_list};")
    print(f"Ki_list = {Ki_list};")
    print(f"Kd_list = {Kd_list};")

    return pareto, Kp_list, Ki_list, Kd_list


if __name__ == "__main__":
    optimize_pid("analy_report.txt")
