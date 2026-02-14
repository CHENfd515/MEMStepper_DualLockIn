# MEMStepper_DualLockIn
MATLAB App for MEMS Stepper Control with Dual Lock‑in Amplifier Readout

📫 alexw25.us@berkeley.edu | alexw25.us@utexas.edu
> By Alex @ Feb.2026

> UC Berkeley-MIT Cao & Tang Lab, All Copyright Reserved.
## Overview
`MEMStepper_DualLockIn` is a MATLAB App Designer application for **precision control and characterization of MEMS stepper actuators**.
It integrates:
- Dual Stanford Research lock‑in amplifier communication (SR830, SR844)
- Capacitive sensing & real‑time fitting
- Camera imaging & motion tracking (angular / radial displacement)
- Open‑loop & `closed‑loop` positioning
- `PID` feedback tuning
- Data logging & visualization

## 1. Main Features
- **MEMS Drive Control**
  - 3‑phase voltage output (V1, V2, V3)
  - Direct voltage / stepping / microstepping modes
  - Vz and Vd bias control
  - Output ON/OFF safety

- **Dual Lock‑in Amplifier Interface**
  - Automatic VISA instrument scanning
  - `SR830` for actuation control
  - `SR844` for capacitive readout
  - Real‑time capacitance / phase plotting

- **Capacitive Motion Detection**
  - `Live readout waveform`
  - Step‑capacitance `linear fitting`
  - `Feedback loss tracing`
  - Data logging to file

- **Camera & Image Processing**
  - USB camera preview, snapshot, video recording
  - Polar coordinate transformation (depolar view)
  - Angular (θ) and radial (R) motion extraction
  - FFT‑based displacement analysis
  - Reference design alignment (experimental)

- **Positioning Modes**
  - Open‑loop target stepping
  - `Closed‑loop capacitance feedback`
  - Scan & auto‑zero routines
  - `PID` parameter sweep & tuning

- **Data & Logging**
  - `Capacitance vs. step fitting`
  - Real‑time plot updates
  - Detailed motion log (`mems.log`)
  - Snapshot & video export

## 2. Hardware Requirements
- software: `MATLAB R2025b`
- 2× Lock‑in amplifiers:
  - SR830 (for driving)
  - SR844 (for sensing)
- VISA backend (e.g., NI‑VISA)
- MEMS stepper actuator

## 3. Core Internal Functions
- `scanUSBPort`: auto‑detect VISA instruments
- `readLockin`: read and parse lock‑in output
- `setVoltages`: 3‑phase drive waveform generation
- `goTo`, `goTo_fit`, `goTo_fb`: positioning engines
- PID control with logging & stability check

## 4. File Descriptions
- `MEMStepper_DualLockIn.mlapp` — main app class (must be edited in **MATLAB R2025b**)
- `MEMStepper_DualLockIn.m` — main app class (can be edited in **MATLAB R2024b**)
- `mems.log` — runtime feedback log (auto‑created)
- Saved snapshots / `capacitance data` / `python analysis script`
