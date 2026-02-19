# MEMStepper_DualLockIn

**MATLAB Application for MEMS Stepper Control with Dual Lock-in Amplifier Readout**

📫 [alexw25.us@berkeley.edu](mailto:alexw25.us@berkeley.edu) | [alexw25.us@utexas.edu](mailto:alexw25.us@utexas.edu)

> Developed by Alex · Feb 2026
> UC Berkeley–MIT Cao & Tang Lab
> © 2026 All Rights Reserved

---
![MEMS](./files/mems.png)

## Overview

`MEMStepper_DualLockIn` ([code](./MEMStepper_DualLockIn.m) | [app](./MEMStepper_DualLockIn.mlapp)) is a MATLAB App Designer–based control and characterization platform for MEMS stepper actuators.


The system integrates:

* Dual lock-in amplifier communication
* Capacitive sensing with real-time linear fitting
* Camera-based motion tracking
* Open-loop and closed-loop positioning
* PID feedback tuning and parameter sweep
* Automated logging and data visualization

The application enables **end-to-end drive, sensing, feedback control, and performance analysis** of MEMS stepper devices within a unified graphical interface.

See the [tutorial](tutorial.pdf) for a quick start guide.

---

# 1. Main Features

| Module                     | Capability                     | Description                                           |
| -------------------------- | ------------------------------ | ----------------------------------------------------- |
| **Drive Control**          | 3-Phase Output (V1/V2/V3)      | Generates programmable three-phase drive waveforms    |
|                            | Stepping / Microstepping Modes | Supports discrete stepping and fine microstep control |
|                            | Bias Control (Vz / Vd)         | Independent bias voltage tuning                       |
|                            | Output Protection              | Safe ON/OFF switching logic                           |
| **Dual Lock-in Interface** | VISA Auto-Detection            | Automatic instrument scanning                         |
|                            | SR830 (Drive)                  | Sensing signal generation                           |
|                            | SR844 (Sense)                  | Capacitive amplitude/phase measurement                |
|                            | Real-Time Plotting             | Continuous waveform visualization                     |
| **Capacitive Feedback**    | Live Capacitance Readout       | Streaming sensing data                                |
|                            | Step–Cap Linear Fit            | Extract k/b calibration parameters                    |
|                            | Loss Tracking                  | Theoretical vs. actual step deviation monitoring      |
|                            | Auto Logging                   | Automatic log file generation                         |
| **Vision Module**          | USB Camera Control             | Preview, capture, and recording                       |
|                            | Polar Transform                | Angular/radial displacement extraction                |
|                            | FFT Analysis                   | Frequency-domain displacement characterization        |
|                            | Reference Alignment            | Experimental alignment tool                           |
| **Positioning Control**    | Open-Loop Mode                 | Direct target stepping                                |
|                            | Closed-Loop Mode               | Capacitance-based feedback positioning                |
|                            | Scan & Auto-Zero               | Automated calibration routine                         |
|                            | PID Sweep & Tuning             | Stability evaluation and parameter optimization       |
| **Data & Analysis**        | Real-Time Plot Update          | Capacitance-step monitoring                           |
|                            | Motion Logging                 | Structured `mems.log` file                            |
|                            | Snapshot/Video Export          | Data archival and documentation                       |

---

# 2. Hardware Requirements

| Category    | Requirement                              |
| ----------- | ---------------------------------------- |
| Software    | MATLAB R2025b (App Designer recommended) |
| Instruments | 2× Lock-in                    |
|             | SR830 (sensing drive control)                    |
|             | SR844 (capacitive sensing)               |
| Interface   | VISA backend (e.g., NI-VISA)             |
| Device      | MEMS Stepper Actuator                    |

---

# 3. Core Internal Functions

| Function       | Description                                          |
| -------------- | ---------------------------------------------------- |
| `scanUSBPort`  | Auto-detect VISA instruments                         |
| `readLockin`   | Parse lock-in amplifier output data                  |
| `setVoltages`  | Generate three-phase drive voltages                  |
| `goTo`         | Open-loop positioning                                |
| `goTo_fit`     | Positioning to calibrate k/b linear model                     |
| `goTo_fb`      | Closed-loop feedback positioning                     |
| `goTo_fb_scan` | Closed-loop positioning (scan-based)  |
| PID Module     | Closed-loop control with logging and stability check |

---

# 4. File Structure

```
├─analyze
│  ├─0_fitscanLog  # Scripts for visualizing and analyzing fit process parameters
│  ├─1_memsLog     # Scripts for reading and analyzing mems.log
│  └─2_capdataLog  # Scripts for visualizing capdata.log
├─docs
│  ├─our_mems      # Article of our MEMS
│  └─ref_doc       # MEMS PID related reference
├─files
│  ├─MEMStepper_DAC_Approach.mlapp  # Original version of the MEMS control app code
│  ├─MEMS_phy.mp4                   # Video demonstrating MEMS rotation angle characterization
│  ├─mems_r_2um.png                 # Optical microscopy image of MEMS device (2μm scale)
│  ├─PID_demo.mp4                   # Demo video for PID operation
│  ├─ScanPID-1.mp4                  # Scan PID control demo video 1
│  └─ScanPID-2.mp4                  # Scan PID control demo video 2
└─latex
```

---

# 5. File Descriptions

| File                                          | Description                                         |
| --------------------------------------------- | --------------------------------------------------- |
| `MEMStepper_DualLockIn.mlapp`                 | Main App Designer class (MATLAB R2025b recommended) |
| `MEMStepper_DualLockIn.m`                     | Script-compatible class file (editable in R2024b)   |
| `mems.log`                                    | Runtime motion feedback log (auto-generated)        |
| Snapshots / Capacitance Data / Python Scripts | Exported experimental data and analysis scripts     |

---


