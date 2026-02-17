import re
import datetime
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def plot_capdata(file_path, debug=True):
    """
    Parses batch-structured MEMS log files, calculates average batch time intervals,
    generates evenly spaced time points for each batch's data points, and plots:
    1. Cap_Readout vs Time (time starts at 0 seconds)
    2. Step vs Time (time starts at 0 seconds)
    
    Args:
        file_path (str): Path to log file
        debug (bool): Whether to print debug information
    """
    # --------------------------
    # Step 1: Regex Patterns for Parsing
    # --------------------------
    # Pattern 1: Match LOG BATCH header (capture batch timestamp and count)
    batch_header_pattern = re.compile(
        r">> LOG BATCH: (?P<batch_ts>\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}\.\d{6}) \| Count: (?P<count>\d+) <<"
    )
    
    # Pattern 2: Match data points (capture Pt ID, Cap_Readout, step)
    data_point_pattern = re.compile(
        r"\[\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}\.\d{6}\] Pt (?P<pt_id>\d+): (?P<Cap_Readout>[\+\-]\d+\.\d+e[\+\-]\d+) V \[step=(?P<step>\s*\d+\.\d+)\]"
    )

    # --------------------------
    # Step 2: Parse Log File
    # --------------------------
    batch_data = []  # Store all batches: [{"batch_ts": datetime, "count": int, "points": [...]}, ...]
    current_batch = None

    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line:
                    continue

                # Match batch header
                batch_match = batch_header_pattern.match(line)
                if batch_match:
                    # Finalize previous batch if exists
                    if current_batch:
                        batch_data.append(current_batch)
                    
                    # Initialize new batch
                    batch_ts_str = batch_match.group('batch_ts')
                    count = int(batch_match.group('count'))
                    current_batch = {
                        "batch_ts": datetime.datetime.strptime(batch_ts_str, '%Y-%m-%d %H:%M:%S.%f'),
                        "count": count,
                        "points": []
                    }
                    if debug:
                        print(f"[DEBUG] Found Batch at line {line_num}: TS={batch_ts_str}, Count={count}")
                    continue

                # Match data point (only if in a batch)
                if current_batch:
                    point_match = data_point_pattern.match(line)
                    if point_match:
                        try:
                            pt_id = int(point_match.group('pt_id'))
                            Cap_Readout = float(point_match.group('Cap_Readout'))
                            step = float(point_match.group('step').strip())
                            current_batch["points"].append({
                                "pt_id": pt_id,
                                "Cap_Readout": Cap_Readout,
                                "step": step
                            })
                        except ValueError as ve:
                            if debug:
                                print(f"[DEBUG] Invalid data point at line {line_num}: {ve}")
                            continue

            # Add the last batch
            if current_batch:
                batch_data.append(current_batch)

    except Exception as e:
        print(f"Error reading/parsing file: {e}")
        return

    if not batch_data:
        print("No batch data found in file.")
        return
    print(f"Successfully parsed {len(batch_data)} batches")

    # --------------------------
    # Step 3: Calculate Time Intervals & Generate Time Points
    # --------------------------
    # Convert batch timestamps to epoch seconds (for calculation)
    batch_timestamps = [b["batch_ts"].timestamp() for b in batch_data]
    all_plot_data = []  # Final data for plotting: [(time_sec, Cap_Readout, step), ...]

    if len(batch_data) == 1:
        # Case 1: Only 1 batch - use arbitrary small interval (0.1s between points)
        batch = batch_data[0]
        num_points = len(batch["points"])
        time_interval = 0.1  # 100ms between points
        for i, point in enumerate(batch["points"]):
            time_sec = i * time_interval
            all_plot_data.append((time_sec, point["Cap_Readout"], point["step"]))
        print(f"Single batch detected - using {time_interval}s interval between points")

    else:
        # Case 2: Multiple batches - calculate average batch interval
        # Compute time differences between consecutive batches
        batch_intervals = []
        for i in range(1, len(batch_timestamps)):
            interval = batch_timestamps[i] - batch_timestamps[i-1]
            batch_intervals.append(interval)
        
        avg_batch_interval = np.mean(batch_intervals)
        print(f"Average batch time interval: {avg_batch_interval:.6f} seconds")

        # Generate evenly spaced time points for each batch's points
        cumulative_time = 0.0  # Start at 0 seconds
        for batch_idx, batch in enumerate(batch_data):
            num_points = len(batch["points"])
            if num_points == 0:
                continue
            
            # Evenly split batch interval into count points
            point_interval = avg_batch_interval / num_points
            
            # Assign time to each point in the batch
            for i, point in enumerate(batch["points"]):
                time_sec = cumulative_time + (i * point_interval)
                all_plot_data.append((time_sec, point["Cap_Readout"], point["step"]))
            
            # Update cumulative time for next batch
            cumulative_time += avg_batch_interval

    # Convert to DataFrame for easy plotting
    df = pd.DataFrame(all_plot_data, columns=["time_sec", "Cap_Readout", "step"])

    # --------------------------
    # Step 4: Generate Plots
    # --------------------------
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), sharex=True)
    plt.subplots_adjust(hspace=0.3)

    # Subplot 1: Cap_Readout vs Time (seconds)
    ax1.plot(df["time_sec"], df["Cap_Readout"], marker='o', linestyle='-', color='blue', alpha=0.7, markersize=2)
    ax1.set_title('Cap_Readout vs Time (Start at 0 Seconds)')
    ax1.set_ylabel('Cap_Readout (V)')
    ax1.grid(True, alpha=0.3)

    # Subplot 2: Step vs Time (seconds)
    ax2.plot(df["time_sec"], df["step"], marker='s', linestyle='-', color='red', alpha=0.7, markersize=2)
    ax2.set_title('Step vs Time (Start at 0 Seconds)')
    ax2.set_xlabel('Time (seconds)')
    ax2.set_ylabel('Step Value')
    ax2.grid(True, alpha=0.3)

    # plt.show()
    date_now = datetime.datetime.now().strftime("%m%d")
    plt.savefig("capdata_plots_{}.png".format(date_now))

# --------------------------
# Example Usage
# --------------------------
plot_capdata("capdata.log", debug=True)