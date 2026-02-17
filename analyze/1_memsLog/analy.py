import re
import pandas as pd
import numpy as np
import time
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.patches import Patch
import datetime


cutPoint = 100  # Last 100 steps for focused analysis


def parse_pid_config(conf_str):
    parts = conf_str.replace(' ', '').split(',')
    kp = float(parts[0].split('=')[1])
    ki = float(parts[1].split('=')[1])
    kd = float(parts[2].split('=')[1])
    return kp, ki, kd

def analyze_mems_log(file_path):
    # Regular expression matching logic
    pattern = (
        r"Step Loss \(.*?\): \[\s*([\d\.-]+)\]"
        r".*?Gains\(P,I,D\): \[\s*([\d\.e\+-]+),\s*([\d\.e\+-]+),\s*([\d\.e\+-]+)\]"
        r".*?PID_Out\(P\+I\+D\): \[\s*([\d\.-]+)\s*\+\s*([\d\.-]+)\s*\+\s*([\d\.-]+)\s*="
    )

    data = []
    last_config = None
    group_id = 0
    
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            for line in f:
                match = re.search(pattern, line)
                if match:
                    step_loss = float(match.group(1))
                    kp, ki, kd = map(float, match.group(2, 3, 4))
                    p_out, i_out, d_out = map(float, match.group(5, 6, 7))
                    config_str = f"Kp={kp:.2e}, Ki={ki:.2e}, Kd={kd:.2e}"

                    if config_str != last_config:
                        group_id += 1
                        last_config = config_str

                    data.append({
                        'Group_ID': group_id,
                        'Step_Loss': step_loss,
                        'Abs_Step_Loss': abs(step_loss),
                        'PID_Config': config_str,
                        'P_Out': p_out, 'I_Out': i_out, 'D_Out': d_out,
                        'is_zero': (kp == 0 and ki == 0 and kd == 0)
                    })
    except Exception as e:
        print(f"Error reading file: {e}"); return

    if not data:
        print("No data found."); return

    df = pd.DataFrame(data)
    unique_group_ids = df['Group_ID'].unique()
    timestamp = time.strftime("%m%d")

    # ==========================================================
    # Synchronization Protocol: Generate analy_report.txt that meets load_report() regex requirements
    # ==========================================================
    report_lines = []
    header = (f"{'Group':<8} | {'Kp':>9} | {'Ki':>9} | {'Kd':>9} | "
          f"{'Mean':>10} | {'Std':>10} | {'MAE':>10} | "
          f"{'Mean_l100':>10} | {'MAE_l100':>10} | {'MSE_l100':>10} | "
          f"{'P%':>6} | {'I%':>6} | {'D%':>6} | {'Dom':>4}")
    report_lines.append(header)
    report_lines.append("-" * len(header))

    for gid in unique_group_ids:
        subset = df[df['Group_ID'] == gid]
        conf = subset['PID_Config'].iloc[0]
        kp, ki, kd = parse_pid_config(conf)
        
        # Basic statistical metrics
        v_mean = subset['Step_Loss'].mean()
        v_std = subset['Step_Loss'].std()
        v_p2p = subset['Step_Loss'].max() - subset['Step_Loss'].min()
        v_mae = subset['Abs_Step_Loss'].mean()
        
        # Last cutPoint steps analysis (or all if less than cutPoint)
        last_subset = subset.tail(cutPoint)
        v_mae_last = last_subset['Abs_Step_Loss'].mean()
        v_mse_last = (last_subset['Step_Loss']**2).mean()
        v_mean_last = last_subset['Step_Loss'].mean()
        
        # Contribution ratio calculation
        p_c, i_c, d_c = np.abs(subset[['P_Out', 'I_Out', 'D_Out']]).mean()
        total_c = p_c + i_c + d_c
        p_r, i_r, d_r = (p_c/total_c, i_c/total_c, d_c/total_c) if total_c > 0 else (0, 0, 0)
        
        # Dominant term determination (must comply with [PID] character set in regex, set to P when all zeros to prevent regex failure)
        dom = ['P', 'I', 'D'][np.argmax([p_r, i_r, d_r])] if total_c > 0 else "P"

        # Strictly construct row according to load_report pattern
        row = (f"G{gid:<7} | {kp:9.2e} | {ki:9.2e} | {kd:9.2e} | "
            f"{v_mean:10.2e} | {v_std:10.2e} | {v_mae:10.2e} | "
            f"{v_mean_last:10.2e} | {v_mae_last:10.2e} | {v_mse_last:10.2e} | "
            f"{p_r*100:6.1f} | {i_r*100:6.1f} | {d_r*100:6.1f} | {dom:>4}")

        report_lines.append(row)

    with open('analy_report.txt', 'w', encoding='utf-8') as f_rep:
        f_rep.write("\n".join(report_lines))
    print("Report 'analy_report.txt' has been generated.")

    # ==========================================================
    # Plotting logic
    # ==========================================================
    fig = plt.figure(figsize=(40, 14))
    ax1 = fig.add_subplot(2, 2, 1); ax2 = fig.add_subplot(2, 2, 2)
    ax3 = fig.add_subplot(2, 2, 3); ax4 = fig.add_subplot(2, 2, 4, polar=True)

    color_cycle = plt.cm.get_cmap('tab20', len(unique_group_ids))
    gid_to_color = {gid: (color_cycle(i) if not df[df['Group_ID'] == gid]['is_zero'].iloc[0] else (0,0,0,1)) 
                    for i, gid in enumerate(unique_group_ids)}

    df['Color'] = df['Group_ID'].map(gid_to_color)

    # Subplot 1: MAE
    mae_analysis = df.groupby('Group_ID')['Abs_Step_Loss'].mean()
    bars1 = ax1.bar(mae_analysis.index, mae_analysis.values, color=[gid_to_color[g] for g in mae_analysis.index], edgecolor='black', alpha=0.7)
    ax1.set_xticks(unique_group_ids)
    ax1.set_xticklabels([f"G{g}" for g in unique_group_ids])
    for gid in unique_group_ids:
        val = mae_analysis[gid]
        conf = df[df['Group_ID'] == gid]['PID_Config'].iloc[0]
        ax1.text(gid, val, f'{val:.2e}', va='bottom', ha='center', fontsize=8)
        ax1.text(gid, val/2, conf.replace(', ', '\n'), va='center', ha='center', fontsize=7, bbox=dict(facecolor='white', alpha=0.6, edgecolor='none'))
    ax1.set_ylabel('Step Loss MAE, where Step = CapReadOut / k')
    ax1.set_title('Mean Absolute Error by PID Group\nmean( abs(error) )')   
    ax1.grid(True, linestyle='--', alpha=0.5)
    
    # Subplot 2: Trend
    ax2.bar(range(len(df)), df['Step_Loss'], color=df['Color'], width=1.0)
    ax2.axhline(0, color='red', lw=1, alpha=0.5)
    ax2.set_ylabel('Step Loss')
    ax2.set_title('Step Loss Trend Across All Groups')
    ax2.grid(True, linestyle='--', alpha=0.5)
    
    # Subplot 3: Boxplot
    group_data = [df[df['Group_ID'] == gid]['Step_Loss'] for gid in unique_group_ids]
    bplot = ax3.boxplot(group_data, patch_artist=True, labels=[f"G{g}" for g in unique_group_ids])
    for patch, gid in zip(bplot['boxes'], unique_group_ids):
        patch.set_facecolor(gid_to_color[gid])
        patch.set_alpha(0.5)
    
    ymin, ymax = ax3.get_ylim()
    for i, gid in enumerate(unique_group_ids):
        std_v = df[df['Group_ID'] == gid]['Step_Loss'].std()
        ax3.text(i+1, ymin + (ymax-ymin)*0.02, f'Std\n{std_v:.1e}', 
                 color=gid_to_color[gid], ha='center', fontweight='bold', fontsize=8)
    ax3.axhline(0, color='grey', lw=1, alpha=0.5, ls='--')
    ax3.set_ylabel('Step Loss')
    ax3.set_title('Boxplot of Step Loss by Group')
    ax3.grid(True, linestyle='--', alpha=0.5)
    
    # Subplot 4: Radar
    angles = np.array([0, 2*np.pi/3, 4*np.pi/3, 0])
    for gid in unique_group_ids:
        sub = df[df['Group_ID'] == gid]
        if sub['is_zero'].iloc[0]: continue
        vals = np.abs(sub[['P_Out', 'I_Out', 'D_Out']]).mean()
        if vals.sum() == 0: continue
        norm_vals = np.append(vals.values / vals.sum(), vals.values[0] / vals.sum())
        ax4.plot(angles, norm_vals, label=f"G{gid}", color=gid_to_color[gid], lw=2)
        ax4.fill(angles, norm_vals, color=gid_to_color[gid], alpha=0.1)
    ax4.set_title('PID Contribution Radar Chart')
    ax4.set_xticks(angles[:-1])
    ax4.set_xticklabels(['P', 'I', 'D'])
    ax4.legend(loc='upper right', bbox_to_anchor=(1.3, 1.1))
    
    plt.tight_layout()
    plt.savefig(f'mems_pid_analysis_{timestamp}.png', dpi=300)

def curve_mems_log_byGroup(file_path, fine_plot=False):
    """
    Plot Step Loss curve subplots grouped by Group_ID
    Save as high-resolution PDF: mems_step_loss_curves_MMDD.pdf

    Parameters
    ----------
    file_path : str
        Log file path
    fine_plot : bool
        Whether to enable fine mode
        True  -> Fixed y-axis range [-0.12, 0.12]
        False -> Auto reasonable scaling
    """
    
    pattern = (
        r"Step Loss \(.*?\): \[\s*([\d\.-]+)\]"
        r".*?Gains\(P,I,D\): \[\s*([\d\.e\+-]+),\s*([\d\.e\+-]+),\s*([\d\.e\+-]+)\]"
        r".*?PID_Out\(P\+I\+D\): \[\s*([\d\.-]+)\s*\+\s*([\d\.-]+)\s*\+\s*([\d\.-]+)\s*="
    )

    data = []
    last_config = None
    group_id = 0

    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            for line in f:
                match = re.search(pattern, line)
                if match:
                    step_loss = float(match.group(1))
                    kp, ki, kd = map(float, match.group(2, 3, 4))

                    config_str = f"Kp={kp:.2e}, Ki={ki:.2e}, Kd={kd:.2e}"

                    if config_str != last_config:
                        group_id += 1
                        last_config = config_str

                    data.append({
                        'Group_ID': group_id,
                        'Step_Loss': step_loss,
                        'PID_Config': config_str
                    })

    except Exception as e:
        print(f"Error reading file: {e}")
        return

    if not data:
        print("No data found.")
        return

    df = pd.DataFrame(data)
    unique_group_ids = df['Group_ID'].unique()
    n_groups = len(unique_group_ids)

    n_cols = int(np.ceil(np.sqrt(n_groups)))
    n_rows = int(np.ceil(n_groups / n_cols))

    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(6*n_cols, 3.5*n_rows)
    )

    axes = np.array(axes).flatten()
    color_map = plt.cm.get_cmap('tab20', n_groups)

    ref_lines = [-0.1, -0.05, -0.01, 0.01, 0.05, 0.1]

    for i, gid in enumerate(unique_group_ids):
        ax = axes[i]
        subset = df[df['Group_ID'] == gid]

        y = subset['Step_Loss'].values
        x = np.arange(len(y))

        y_max = np.max(y)
        y_min = np.min(y)
        mae = np.mean(np.abs(y))
        if len(y) >= cutPoint:
            mse_100 = np.mean(y[-cutPoint:]**2)
        else:
            mse_100 = np.mean(y**2)


        # =============================
        # Y-axis control logic
        # =============================
        if fine_plot:
            lower, upper = -0.12, 0.12
        else:
            y_range = y_max - y_min
            margin = 0.1 * y_range if y_range != 0 else 0.01
            lower = min(y_min - margin, min(ref_lines) * 1.05)
            upper = max(y_max + margin, max(ref_lines) * 1.05)

        ax.set_ylim(lower, upper)

        # Main curve
        ax.plot(x, y, linewidth=1.5, color=color_map(i))

        # Zero line
        ax.axhline(0, color='red', lw=1, alpha=0.5)

        # Fixed reference lines
        for val in ref_lines:
            ax.axhline(val, linestyle='--', linewidth=0.8, alpha=0.4, color='gray')

        pid_label = subset['PID_Config'].iloc[0]

        # ax.set_title(
        #     f"G{gid} | {pid_label}\n"
        #     f"max={y_max:.2e}  min={y_min:.2e}  MAE={mae:.2e}",
        #     fontsize=9
        # )
        ax.set_title(
            f"G{gid} | {pid_label}\n"
            f"max={y_max:.2e}  min={y_min:.2e}  "
            f"MAE={mae:.2e}  MSE(Last{cutPoint:d})={mse_100:.2e}",
            fontsize=9
        )


        ax.set_xlabel("Step Index")
        ax.set_ylabel("Step Loss")
        ax.grid(True, linestyle=':', alpha=0.3)

    # Remove empty subplots
    for j in range(i+1, len(axes)):
        fig.delaxes(axes[j])

    timestamp = time.strftime("%m%d")
    mode = "fine" if fine_plot else "auto"
    filename = f"mems_step_loss_curves_{mode}_{timestamp}"

    plt.suptitle(
        f"Step Loss Trend by PID Group ({'Fine' if fine_plot else 'Auto'} Mode)",
        fontsize=16
    )

    plt.tight_layout(rect=[0, 0, 0.97, 0.95])
    plt.savefig(filename+".pdf", format='pdf')
    plt.savefig(filename+".png", format='png', dpi=300)
    plt.close()

    print(f"Saved: {filename}")

def curve_mems_log_byStep(file_path, debug=True, color_type="discrete"):
    """
    Parses MEMS log file, organizes data by theoretical step groups,
    and plots 4 subplots: Step, Capacitance, Loss over time.
    """
    
    # Regex pattern to parse valid log lines
    # Groups: timestamp, theor_cap, act_cap, cap_loss, step_loss, act_step, theor_step
    pattern = re.compile(
        r"Cap\(theor-act\): \[\s*(?P<theor_cap>[\d\.e\+-]+)\s*-\s*(?P<act_cap>[\d\.e\+-]+)\s*=\s*(?P<cap_loss>[\d\.e\+-]+)\]"
        r".*?Step Loss \(=Cap/k_fit\): \[\s*(?P<step_loss>[\d\.e\+-]+)\]"
        r".*?@\s*(?P<timestamp>\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}\.\d{6})\s*@"
        r".*?Step \[(?P<act_step>[\d\.e\+-]+)/(?P<theor_step>[\d\.e\+-]+)\]"
    )

    data = []
    match_count = 0

    try:
        line_count = 0
        with open(file_path, 'r', encoding='utf-8') as f:
            for line in f:
                line_count += 1
                line = line.strip()
                if not line:
                    continue

                # Skip lines that do not contain capacitance data
                if "Cap(theor-act)" not in line:
                    continue

                # Try to match the log pattern
                match = pattern.search(line)
                if match:
                    match_count += 1
                    try:
                        d = match.groupdict()
                        # Convert parsed strings to proper data types
                        d['timestamp']   = datetime.datetime.strptime(d['timestamp'], '%Y-%m-%d %H:%M:%S.%f')
                        d['theor_cap']   = float(d['theor_cap'])
                        d['act_cap']     = float(d['act_cap'])
                        d['cap_loss']    = float(d['cap_loss'])
                        d['step_loss']   = float(d['step_loss'])
                        d['theor_step']  = float(d['theor_step'])
                        d['act_step']    = float(d['act_step'])

                        data.append(d)
                    except ValueError as ve:
                        if debug:
                            print(f"Value error at line {line_count}: {ve}")
                else:
                    if debug:
                        print(f"Regex failed at line {line_count}: {line}")

    except Exception as e:
        print(f"File read error: {e}")
        return

    if not data:
        print("No valid data parsed.")
        return

    print(f"Done. Lines read: {line_count}, Matched: {match_count}")
    if debug:
        print(f"Sample parsed entry: {data[0]}")
    # ===== Convert timestamp to relative seconds (start from 0) =====
    df = pd.DataFrame(data)
    df = df.sort_values('timestamp') 
    t0 = df['timestamp'].iloc[0]
    df['time_sec'] = (df['timestamp'] - t0).dt.total_seconds()


    # Create 4 subplots
    fig, axs = plt.subplots(4, 1, figsize=(12, 16))
    plt.subplots_adjust(hspace=0.4)

    # Group data by theoretical step for coloring
    groups = df.groupby('theor_step')
    
    
    # # Use 32 fixed colors and cycle through them if groups exceed 32
    # colors = plt.cm.get_cmap('tab20', 32)
    # color_list = [colors(i % 32) for i in range(len(groups))]
    
    if color_type == "continus_eachStep" or color_type == "continus":
        if color_type == "continus_eachStep":
            continus_color = True
        else:   
            continus_color = False
        step = df['theor_step'].values
        # Integer part: 0~9
        step_int = np.floor(step).astype(int)
        # Fraction inside each block (0~1)
        step_frac = step - step_int
        # Normalize block index globally (0~1)
        block_norm = step_int / 10.0   # since range is 0~10
        # Base colormap (scientific-friendly)
        base_cmap = plt.cm.viridis   # green → blue → purple

        colors = []
        for i in range(len(step)):
            # Shift along global colormap by block index
            base_color_value = 0.15 + 0.7 * block_norm[i]
            # Add local gradient within block
            if continus_color:
                local_shift = 0.15 * step_frac[i]
                final_value = min(base_color_value + local_shift, 1.0)
            else:
                final_value = base_color_value

            colors.append(base_cmap(final_value))
            color_list = np.array(colors)
            
    elif color_type == "discrete":
        colors = plt.cm.get_cmap('tab20', 32)
        color_list = [colors(i % 32) for i in range(len(groups))]
            

    
    print(f"Unique theoretical steps: {len(groups)}")
    # Plot 1: Actual vs Theoretical Step over time
    for i, (step_val, group) in enumerate(groups):
        c = color_list[i]
        axs[0].plot(group['time_sec'], group['theor_step'], '--', color=c, alpha=0.6)
        axs[0].plot(group['time_sec'], group['act_step'], '-', color=c, label=f'{step_val:.2e}')
    axs[0].set_title('Step Position vs Time')
    axs[0].set_ylabel('Step')
    # axs[0].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    axs[0].grid(True)

    # Plot 2: Actual vs Theoretical Capacitance
    for i, (step_val, group) in enumerate(groups):
        c = color_list[i]
        axs[1].plot(group['time_sec'], group['theor_cap'], '--', color=c, alpha=0.6)
        axs[1].plot(group['time_sec'], group['act_cap'], '-', color=c)
    axs[1].set_title('Capacitance vs Time')
    axs[1].set_ylabel('Capacitance')
    # axs[1].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    axs[1].grid(True)
    
    # Plot 3: Capacitance Loss (Theor - Act)
    for i, (step_val, group) in enumerate(groups):
        c = color_list[i]
        axs[2].bar(group['time_sec'], group['cap_loss'], width=0.5, color=c, alpha=0.7)
    axs[2].set_title('Capacitance Loss (Theor - Actual)')
    axs[2].set_ylabel('Cap Loss')
    # axs[2].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    axs[2].grid(True)

    # Plot 4: Step Loss over time
    for i, (step_val, group) in enumerate(groups):
        c = color_list[i]
        axs[3].bar(group['time_sec'], group['step_loss'], width=0.5, color=c, alpha=0.7)
    axs[3].set_title('Step Loss vs Time')
    axs[3].set_ylabel('Step Loss')
    axs[3].set_xlabel('Time (sec)')
    # axs[3].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
    axs[3].grid(True)

    # plt.show()
    plt.savefig(f'mems_step_loss_curves_byStep_{time.strftime("%m%d")}.png', dpi=300)

if __name__ == "__main__":
    analyze_mems_log('mems.log')
    # curve_mems_log_byGroup('mems.log', fine_plot=False)
    # curve_mems_log_byGroup('mems.log', fine_plot=True)
    # curve_mems_log_byStep('mems.log', debug=False)