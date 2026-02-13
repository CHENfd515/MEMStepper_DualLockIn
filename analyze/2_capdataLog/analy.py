import time
import pandas as pd
import matplotlib.pyplot as plt
import re
import os

def parse_cap_log(file_path):
    """
    Parses the capdata.log file.
    Extracts timestamps, indices, values, and units.
    """
    data = []
    # Regex to match: [2026-02-04 14:30:00] Pt 000005: +1.23456e-01 uV
    pattern = r"\[(.*?)\] Pt (\d+): ([\d\.e\+-]+) (\w+)"
    
    if not os.path.exists(file_path):
        print(f"Error: {file_path} not found.")
        return None

    with open(file_path, 'r') as f:
        for line in f:
            match = re.search(pattern, line)
            if match:
                timestamp = match.group(1)
                idx = int(match.group(2))
                val = float(match.group(3))
                unit = match.group(4)
                data.append({
                    'Timestamp': timestamp,
                    'Index': idx,
                    'Value': val,
                    'Unit': unit
                })
    
    return pd.DataFrame(data)

def plot_data(df, start_idx=None, end_idx=None):
    """
    Plots the data within a specific index range.
    """
    if df is None or df.empty:
        print("No data to plot.")
        return

    # Filter data based on custom range
    if start_idx is None: start_idx = df['Index'].min()
    if end_idx is None: end_idx = df['Index'].max()
    
    mask = (df['Index'] >= start_idx) & (df['Index'] <= end_idx)
    filtered_df = df.loc[mask]

    if filtered_df.empty:
        print(f"No data found in range [{start_idx} : {end_idx}]")
        return

    # Create Plot
    plt.figure(figsize=(12, 6))
    plt.plot(filtered_df['Index'], filtered_df['Value'], 
             marker='o', markersize=3, linestyle='-', linewidth=1, color='#2c3e50')
    
    # Formatting
    unit = filtered_df['Unit'].iloc[0]
    plt.title(f'Capacitance Data Analysis (Index {start_idx} to {end_idx})', fontweight='bold')
    plt.xlabel('Sample Index')
    plt.ylabel(f'Value ({unit})')
    plt.grid(True, linestyle='--', alpha=0.7)
    
    # Annotate first and last timestamp of the selection
    plt.annotate(f"Start: {filtered_df['Timestamp'].iloc[0]}", 
                 xy=(0.02, 0.95), xycoords='axes fraction', fontsize=9, verticalalignment='top')
    plt.annotate(f"End: {filtered_df['Timestamp'].iloc[-1]}", 
                 xy=(0.02, 0.90), xycoords='axes fraction', fontsize=9, verticalalignment='top')

    plt.tight_layout()
    
    timestamp = time.strftime("%m%d")
    plt.savefig('capdata_analysis_{}.png'.format(timestamp), dpi=300)
    # plt.show()
    

if __name__ == "__main__":
    LOG_FILE = 'capdata.log'
    
    # 1. Load data
    df_logs = parse_cap_log(LOG_FILE)
    
    if df_logs is not None:
        print(f"Total points loaded: {len(df_logs)}")
        print(f"Index range: {df_logs['Index'].min()} to {df_logs['Index'].max()}")
        
        # 2. Set your custom range here
        # Example: start_idx=100, end_idx=500
        # Leave as None to plot everything
        plot_data(df_logs, start_idx=None, end_idx=None)