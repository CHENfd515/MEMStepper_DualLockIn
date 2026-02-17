import time
import re
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import numpy as np
from scipy.stats import linregress


def calc_linearity(x, y):
    slope, intercept, _, _, _ = linregress(x, y)
    max_dev = np.max(np.abs(y - (slope * x + intercept)))
    y_range = np.max(y) - np.min(y)
    return (max_dev / y_range) * 100 if y_range != 0 else 0


def analyze_mems_wait_time(file_path):

    with open(file_path, 'r', encoding='utf-8') as f:
        content = f.read()

    sessions = re.split(r'\[LOG\] goTo_fit', content)[1:]
    print(f"Total sessions found: {len(sessions)}")

    experiment_data = []
    all_k_zero = True
    for session in sessions:

        wt_match = re.search(r'wait_time=([\d\.]+)', session)
        if not wt_match:
            print("[Warning] Could not find wait_time, skipping session.")
            continue
        wait_time = float(wt_match.group(1))

        k_match = re.search(r'k_total = ([\d\.e\+-]+)', session)
        k_total = float(k_match.group(1)) if k_match else 1.0
        if abs(k_total) > 1e-12:
            all_k_zero = False  # k_total is not zero, we have valid data to analyze

        dx_match = re.search(
            r'Mean Delta_Step \(dx\) = ([\d\.]+) steps', session)
        dx_val = float(dx_match.group(1)) if dx_match else 0.0

        raw_pts = re.findall(r'\(([\d\.]+),\s*([\d\.e\+-]+)\)', session)

        if not raw_pts:
            print(
                f"[Error] No raw (x, y) points found for {wait_time}s session. [session content: {session[:200]}...]")
            continue

        x_vals = np.array([float(p[0]) for p in raw_pts])
        y_vals = np.array([float(p[1]) for p in raw_pts])

        total_pts = len(x_vals)

        # ===== 2N points div =====
        if total_pts % 2 != 0:
            print(
                f"[Warning] Odd number of points ({total_pts}) at {wait_time}s, skipping.")
            continue

        half = total_pts // 2
        xf, yf = x_vals[:half], y_vals[:half]
        xb, yb = x_vals[half:], y_vals[half:]

        # =========================
        # fit all points to get the overall slope and intercept for theoretical line
        slope_all, intercept_all, _, _, _ = linregress(x_vals, y_vals)

        # theoretical line values
        y_theory_f = slope_all * xf + intercept_all
        y_theory_b = slope_all * xb + intercept_all
        # loss
        fw_loss = (yf - y_theory_f) / \
            abs(k_total) if abs(k_total) > 1e-12 else np.zeros_like(yf)
        fb_loss = (yb - y_theory_b) / \
            abs(k_total) if abs(k_total) > 1e-12 else np.zeros_like(yb)

        experiment_data.append({
            'wait_time': wait_time,
            'k_total': k_total,
            'dx': dx_val,
            'xf': xf,
            'yf': yf,
            'xb': xb,
            'yb': yb,
            'dy_k': (np.flip(yb) - yf) / abs(k_total) if abs(k_total) > 1e-12 else np.zeros_like(yf),
            'lin_f': calc_linearity(xf, yf),
            'lin_b': calc_linearity(xb, yb),
            'fw_loss': fw_loss,
            'fb_loss': fb_loss,
            'slope': slope_all,
            'intercept': intercept_all
        })
    
    if len(experiment_data) == 0:
        print("No valid sessions found.")
        return
    # sort by wait_time
    experiment_data.sort(key=lambda x: x['wait_time'], reverse=True)

    num_sessions = len(experiment_data)
    
    # colors = plt.cm.viridis(np.linspace(0.1, 0.8, num_sessions))
    # colors = plt.cm.viridis(np.linspace(0.5, 0.9, num_sessions))
    colors = plt.cm.viridis(np.linspace(0.25, 0.8, num_sessions))
    bw_color = 'red'  # Color for backward data points and bars
    
    
    if all_k_zero:
        print("All k_total are zero. Plotting only Subplot 1.")
        rows = 1
        figsize = (4 * num_sessions, 12)
    else:
        print("k_total found. Plotting full report.")
        rows = 5
        if num_sessions == 1:
            figsize = (70, 12)
        else:
            figsize = (4 * num_sessions, 12)
    fig = plt.figure(figsize=figsize)
    plt.rc('axes', axisbelow=True)
    fig.patch.set_facecolor('white')
    rows, cols = 5, num_sessions
    

    all_y = np.concatenate(
        [d['yf'] for d in experiment_data] + [d['yb'] for d in experiment_data])
    y_min_val, y_max_val = all_y.min() * 0.995, all_y.max() * 1.005

    # --- Row 1: Subplot 1 (Scatters) - Keep Minor Grid ---
    for i, data in enumerate(experiment_data):
        ax = plt.subplot2grid((rows, cols), (0, i))
        ax.scatter(data['xf'], data['yf'], color=colors[i],
                   marker='o', s=30, label='Forward')
        ax.scatter(data['xb'], data['yb'], color=bw_color,
                   marker='x', s=40, label='Backward')
        ax.set_ylim(y_min_val, y_max_val)

        formatter = ScalarFormatter(useMathText=True)
        formatter.set_scientific(True)
        formatter.set_powerlimits((0, 0))
        ax.yaxis.set_major_formatter(formatter)
        ax.yaxis.get_offset_text().set_fontsize(8)


        ax.set_title(f'Wait: {data["wait_time"]}s',
                     fontweight='bold', fontsize=10)
        ax.set_xlabel('Position (Step)')

        # ax.legend(fontsize=8)
        k = data['slope']
        b = data['intercept']
        from matplotlib.lines import Line2D

        dummy = Line2D([], [], linestyle='none', marker='')

        ax.legend(
            handles=[
                ax.collections[0],   # Forward scatter
                ax.collections[1],   # Backward scatter
                dummy
            ],
            labels=[
                'Forward',
                'Backward',
                f'k={k:.3e}\nb={b:.3e}'
            ],
            fontsize=8
        )

        # ax.legend(
        #     [
        #         f'Forward',
        #         f'Backward\n(k={k:.3e}, b={b:.3e})'
        #     ],
        #     fontsize=8
        # )

        if i == 0:
            ax.set_ylabel('V', fontweight='bold')
        ax.grid(True, which='both', linestyle='--', alpha=0.3)
        ax.minorticks_on()

    if not all_k_zero:
        # --- Row 2: Subplot 2 (Point-by-Point Abs Error %) - Only Major Grid + Horizontal Text ---
        ax_err = plt.subplot2grid((rows, cols), (1, 0), colspan=cols)
        x_pos, x_indices = experiment_data[0]['xf'], np.arange(
            len(experiment_data[0]['xf']))
        bar_width, log_floor = 0.8 / num_sessions, 0.1

        for i, data in enumerate(experiment_data):
            offset = (i - num_sessions/2) * bar_width + bar_width/2
            dy_k_pct = np.abs(data['dy_k']) * 100
            plot_vals = [max(v, log_floor) for v in dy_k_pct]

            # Plot bars for one set of wait_time
            bars = ax_err.bar(x_indices + offset, plot_vals,
                              width=bar_width, color=colors[i], alpha=0.8)
            # bars = ax_err.bar(
            #     x_indices + offset,
            #     plot_vals,
            #     width=bar_width,
            #     color=colors[i],
            #     alpha=0.8,
            #     label=f'{data["wait_time"]}s'
            # )


            for j, bar in enumerate(bars):
                current_h = bar.get_height()
                label_val = abs(data['dy_k'][j]) * 100

                ax_err.text(bar.get_x() + bar.get_width()/2.,
                            current_h * 1.25,
                            f'{label_val:.2f}%',
                            ha='center',
                            va='bottom',
                            fontsize=5,
                            rotation=90,
                            fontweight='bold')

        # ax_err.legend(fontsize=8, loc='upper right', title='Wait Time')

        ax_err.set_yscale('log')
        max_pct = max([np.max(np.abs(d['dy_k'])*100) for d in experiment_data])
        # Increase upper limit to fit horizontal text
        ax_err.set_ylim(bottom=log_floor, top=max_pct * 30)
        ax_err.set_title('abs (Backward - Forward) (Step %)',
                         fontsize=12, fontweight='bold')
        ax_err.set_ylabel('Step (%)', fontweight='bold')
        ax_err.set_xticks(x_indices)
        ax_err.set_xticklabels([f"{x:.1f}" for x in x_pos])
        ax_err.set_xlabel('Position (Step)')
        ax_err.grid(True, which='major', linestyle='-',
                    alpha=0.4)  # Only Major Grid

        # --- Row 3: Subplot 3 (Mean dx %) - Only Major Grid ---
        ax_trend = plt.subplot2grid((rows, cols), (2, 0), colspan=cols)
        wt_labels = [f"{d['wait_time']}s" for d in experiment_data]
        dx_pct_vals = [d['dx'] * 100 for d in experiment_data]
        bars_trend = ax_trend.bar(
            wt_labels, dx_pct_vals, color=colors, alpha=0.7, width=0.4)

        ax_trend.set_yscale('log')
        ax_trend.set_ylim(bottom=min(dx_pct_vals)*0.5, top=max(dx_pct_vals)*8)
        ax_trend.set_ylabel('Step (%)', fontweight='bold')
        ax_trend.set_xlabel('Wait Time')
        ax_trend.set_title(
            'mean (abs (Backward - Forward)) (Step %)', fontsize=12, fontweight='bold')
        ax_trend.grid(True, which='major', linestyle='-',
                      alpha=0.4)  # Only Major Grid
        for bar in bars_trend:
            ax_trend.text(bar.get_x() + bar.get_width()/2., bar.get_height() * 1.15,
                          f'{bar.get_height():.3f}%', ha='center', fontweight='bold', fontsize=10)

        # --- Row 4: Subplot 4 (Linearity) - Only Major Grid + Data on Top ---
        ax_lin = plt.subplot2grid((rows, cols), (3, 0), colspan=cols)
        x_wt = np.arange(len(wt_labels))
        # bars_f = ax_lin.bar(x_wt - 0.2, [d['lin_f'] for d in experiment_data], width=0.4, label='Forward', color='navy', alpha=0.7)
        # bars_b = ax_lin.bar(x_wt + 0.2, [d['lin_b'] for d in experiment_data], width=0.4, label='Backward', color=bw_color, alpha=0.7)

        bars_f = ax_lin.bar(
            x_wt - 0.2,
            [d['lin_f'] for d in experiment_data],
            width=0.4,
            label='Forward',
            color=colors,
            alpha=0.8
        )

        bars_b = ax_lin.bar(
            x_wt + 0.2,
            [d['lin_b'] for d in experiment_data],
            width=0.4,
            label='Backward',
            color=bw_color,
            alpha=0.8
        )

        # Write data on top of Subplot 4
        for bars in [bars_f, bars_b]:
            for bar in bars:
                height = bar.get_height()
                ax_lin.text(bar.get_x() + bar.get_width()/2., height + 0.01,
                            f'{height:.2f}%', ha='center', va='bottom', fontsize=9, fontweight='bold')

        ax_lin.set_title('Non-linearity Error (%)',
                         fontsize=12, fontweight='bold')
        ax_lin.set_xticks(x_wt)
        ax_lin.set_xticklabels(wt_labels)
        ax_lin.set_ylabel('Error (%)', fontweight='bold')
        ax_lin.set_xlabel('Wait Time')
        ax_lin.grid(True, which='major', linestyle='-',
                    alpha=0.3)  # Only Major Grid
        # ax_lin.legend(fontsize=8, loc='upper right')
        ax_lin.legend(
            handles=[bars_b],
            labels=['Backward'],
            fontsize=8,
            loc='upper right'
        )
        max_l = max([max(d['lin_f'], d['lin_b']) for d in experiment_data])
        min_l = min([min(d['lin_f'], d['lin_b']) for d in experiment_data])
        ax_lin.set_ylim(min_l * 0.8, max_l * 1.2)

        # --- Row 5: Theory-based FW / FB Loss (Split by wait_time) ---
        for i, data in enumerate(experiment_data):

            ax_loss = plt.subplot2grid((rows, cols), (4, i))
            if np.allclose(data['xf'], data['xb'][::-1]):
                fb_aligned = data['fb_loss'][::-1]
                x_aligned = data['xf']
            else:
                fb_aligned = data['fb_loss']
                x_aligned = data['xf']

            bar_width = 0.35
            x_indices = np.arange(len(x_aligned))

            # bars
            ax_loss.bar(
                x_indices - bar_width/2,
                data['fw_loss'],
                width=bar_width,
                alpha=0.8,
                color=colors[i],
                # label='Forward',
            )

            ax_loss.bar(
                x_indices + bar_width/2,
                fb_aligned,
                width=bar_width,
                alpha=0.8,
                color=bw_color,
                label='Backward',
            )


            ax_loss.axhline(0, linestyle='--')
            ax_loss.set_title(f'Wait: {data["wait_time"]}s',
                              fontweight='bold',
                              fontsize=10)

            ax_loss.set_xlabel('Position (Step)', fontweight='bold')
            ax_loss.set_ylabel('Loss (Step)', fontweight='bold')
            ax_loss.set_xticks(x_indices)
            ax_loss.set_xticklabels([f"{x:.1f}" for x in x_aligned])
            ax_loss.grid(True, which='major', linestyle='-', alpha=0.4)
            
            formatter = ScalarFormatter(useMathText=True)
            formatter.set_scientific(True)
            formatter.set_powerlimits((0, 0))

            ax_loss.xaxis.set_major_formatter(formatter)
            ax_loss.yaxis.set_major_formatter(formatter)

            ax_loss.ticklabel_format(axis='both', style='sci', scilimits=(0, 0))
            ax_loss.xaxis.get_offset_text().set_fontsize(8)
            ax_loss.yaxis.get_offset_text().set_fontsize(8)

            if i == 0:
                ax_loss.legend(fontsize=8)

    # plot layout adjustments
    plt.tight_layout(pad=3)
    if not all_k_zero:
        plt.subplots_adjust(hspace=0.66)

    # Save
    mmdd = time.strftime("%m%d")
    plt.savefig('scan_analy_{}.svg'.format(mmdd),
                format='svg', bbox_inches='tight')
    plt.savefig('scan_analy_{}.pdf'.format(mmdd),
                format='pdf', bbox_inches='tight')
    plt.savefig('scan_analy_{}.png'.format(mmdd),
                format='png', bbox_inches='tight', dpi=300)
    # plt.show()


if __name__ == "__main__":
    # ------------ Usage Instructions ---------------------------
    # Assuming 'scan.log' exists in the same directory
    # scan.log: processed log file containing multiple sessions of MEMS stepper data
    # scan.txt: original log file with raw data (not used in this script)
    # -----------------------------------------------------------
    
    analyze_mems_wait_time('scan.log')
