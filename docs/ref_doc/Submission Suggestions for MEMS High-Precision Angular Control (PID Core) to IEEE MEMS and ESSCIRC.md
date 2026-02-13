# Submission Suggestions for MEMS High-Precision Angular Control (PID Core) to IEEE MEMS and ESSCIRC

The MEMS control software focuses on achieving high-precision angular control, with PID (Proportional-Integral-Derivative) control as its core technology. Combined with the characteristics of the two top conferences (IEEE MEMS and ESSCIRC) and the key points of MEMS interface circuit design, this document provides targeted submission suggestions to help highlight the research value and improve the acceptance rate. The paper submission deadline is April 3, 2026.

# 1. General Background

In MEMS interface circuits, PID control is a classic solution for implementing closed-loop systems, which is particularly common in the scanning control of MEMS micromirrors (key components for high-precision angular control) and the constant-temperature control of pressure sensors. For your research on high-precision angular control, the core is to highlight how the designed PID control scheme solves the key challenges of MEMS devices (such as low mechanical frequency, sensitivity to phase delay, and actuator voltage limitations) and achieves superior control performance, which is the core focus of both conferences.

# 2. Submission Suggestions for ESSCIRC (ESSERC)

ESSCIRC focuses on the efficiency of integrated circuit implementation rather than pure control theory. For your MEMS high-precision angular control system with PID core, the submission should shift the focus to the circuit-level implementation of PID control and its optimization for MEMS characteristics.

## 2.1 Key Focus on Two PID Implementation Routes

Reviewers at ESSCIRC pay close attention to the circuit implementation scheme, performance indicators, and engineering feasibility of PID controllers. You need to clarify the implementation route of your PID control and highlight the corresponding advantages:

### 2.1.1 Analog PID Implementation

Implementation Method: Using operational transconductance amplifiers (OTAs), capacitors, and resistors to implement addition, integration, and differentiation operations, which is suitable for scenarios requiring ultra-low latency (critical for high-precision angular control of MEMS micromirrors).

ESSCIRC Focus Points: The core competitiveness lies in **low power consumption** and **low drift**. Analog PID circuits are highly sensitive to temperature and process variations (PVT: Process, Voltage, Temperature). If your design can achieve an ultra-low-power analog PID circuit with high robustness to PVT variations, it will be highly favored by reviewers. It is recommended to provide specific test data (such as power consumption under typical operating conditions, drift values in different temperature ranges) to verify the performance advantages.

### 2.1.2 Digital/Mixed-Signal PID Implementation

Implementation Method: Composed of ADC (Analog-to-Digital Converter), digital PID logic (DSP/FPGA/Logic), and DAC (Digital-to-Analog Converter), which is the mainstream implementation route in current MEMS control systems, especially suitable for scenarios requiring flexible parameter adjustment and complex control logic.

ESSCIRC Focus Points: Emphasize **area efficiency** and **latency**. Although the mechanical frequency of MEMS devices is relatively low, closed-loop control for high-precision angular control is extremely sensitive to phase delay—excessive latency will lead to system oscillation and reduce angular control accuracy. In addition, the low-power design of ADC and DAC is also a key evaluation index. It is recommended to provide specific indicators such as chip area occupied by PID control logic, end-to-end latency of the control loop, and power consumption of ADC/DAC modules.

## 2.2 Value-Added Points in the Paper: MEMS-Specific Optimization

A standard PID equation implemented on a chip is not sufficient to stand out. You need to demonstrate customized design optimization for the characteristics of MEMS high-precision angular control, which are important innovation points to improve the paper quality:

1. **Anti-Windup Control**: The driving voltage of MEMS actuators (DAC output) has an upper limit. When the integral term accumulates excessively and causes saturation, the system will produce severe oscillation, which seriously affects the high-precision angular positioning. Adding anti-windup logic in the PID circuit (such as clipping the integral term or resetting the integral component when saturation is detected) is an important support for the engineering feasibility of the system. It is recommended to provide simulation and test data to verify that the anti-windup design can effectively suppress system oscillation.

2. **Nonlinear Compensation**: Many MEMS actuators (such as electrostatic drives) adopt an electrostatic driving mode, where the driving force has a square relationship with the voltage. This nonlinear characteristic will introduce errors into high-precision angular control. If your PID circuit has a built-in **linearization processing unit** (such as a square root circuit to compensate for the square relationship of electrostatic driving), it will be a significant innovation point. It is recommended to compare the angular control accuracy with and without nonlinear compensation to highlight the effect of the design.

3. **Online Parameter Self-Tuning**: The mechanical quality factor (Q value) of MEMS devices may change with working conditions (such as temperature, load), which will affect the control performance. If your chip can automatically adjust PID parameters (Kp, Ki, Kd) according to the variation of the MEMS device’s Q value to maintain stable high-precision angular control, it will be a top-conference-level contribution. It is recommended to demonstrate the self-tuning mechanism and test results under different Q value conditions.

## 2.3 Suggestions for Performance Comparison

In the Comparison Table of the paper (a necessary part for ESSCIRC submission), you must list the following key indicators to highlight the advantages of your design compared with the state-of-the-art (SOTA) works. All indicators should be closely related to high-precision angular control:

- **Steady-State Error**: The positioning error of the MEMS device after the closed-loop control reaches stability, which directly reflects the angular control precision (e.g., in μrad level).

- **Settling Time**: The time required for the MEMS device to reach a stable state after receiving a control command, which proves the response speed of the control loop (critical for high-speed and high-precision angular scanning scenarios).

- **Bandwidth**: The effective operating frequency range of the closed-loop control system, which reflects the adaptability of the system to different angular control speeds.

- **Energy Efficiency (Controller Power/Bandwidth)**: The power consumption of the PID control logic divided by the system bandwidth, which reflects the energy efficiency of the circuit implementation and is a key indicator for ESSCIRC.

## 2.4 Suggestions for Key Figures

Key figures are crucial for demonstrating the performance of the control system. ESSCIRC reviewers pay close attention to the stability and reliability of the circuit. It is recommended to include the following figures:

- **Step Response**: Show the comparison of ringing phenomena before and after closed-loop control. It can intuitively reflect the effect of PID control on suppressing system oscillation and improving settling time, which is closely related to angular control stability.

- **Bode Plot**: Focus on showing the **Phase Margin**. A sufficient phase margin (usually ≥45°) proves that the system will not be unstable under various working conditions (such as PVT variations), which is an important guarantee for the reliability of high-precision angular control.

# 3. Submission Suggestions for IEEE MEMS

IEEE MEMS is a top conference focusing on the entire field of MEMS, including device design, fabrication, packaging, and system integration. Compared with ESSCIRC, IEEE MEMS pays more attention to the integration of MEMS devices and control systems, as well as the application value of the technology in high-precision scenarios. For your high-precision angular control software with PID core, the submission should balance the circuit implementation of PID and the integration effect with MEMS devices.

## 3.1 Core Focus Points

1. **Integration of PID Control and MEMS Devices**: Highlight how the designed PID control system is tailored to the characteristics of MEMS micromirrors (or other MEMS devices for angular control), such as the matching between the control loop parameters and the mechanical characteristics of the MEMS device, and the interface design between the PID controller and the MEMS actuator/sensor. It is recommended to provide the integration scheme diagram and test results of the entire system (MEMS device + PID controller).

2. **Angular Control Precision and Application Prospects**: IEEE MEMS attaches great importance to the application value of research. You need to clearly demonstrate the high-precision angular control performance achieved by the system (such as specific steady-state error, angular resolution) and its application scenarios (such as optical communication, precision measurement, medical imaging). It is recommended to provide experimental results of actual angular control (such as the measured angular positioning data of the MEMS micromirror).

3. **Innovation of Control Scheme**: On the basis of PID control, highlight the innovative points combined with MEMS characteristics, which can be consistent with the MEMS-specific optimization points mentioned in ESSCIRC (anti-windup, nonlinear compensation, online parameter self-tuning). However, it is necessary to focus more on how these innovations improve the integration effect and control precision of the MEMS system, rather than just the circuit-level performance.

## 3.2 Supplementary Suggestions for Figures and Indicators

- **System Integration Diagram**: Clearly show the composition of the entire MEMS high-precision angular control system, including MEMS sensors (for angular detection), PID controllers (circuit implementation), MEMS actuators, and interface circuits, to reflect the integration capability.

- **Angular Control Precision Test Curve**: Show the variation curve of the MEMS device’s angular position over time, and mark key indicators such as steady-state error and settling time to intuitively demonstrate the high-precision control effect.

- **Environmental Robustness Test**: Supplement the test results of the control system under different environmental conditions (such as temperature, humidity), which reflects the reliability of the system in practical applications—an important evaluation point for IEEE MEMS.

## 3.3 Difference from ESSCIRC

Compared with ESSCIRC (which focuses on circuit implementation efficiency), IEEE MEMS allows a proper reduction in the detailed description of circuit-level design (such as the internal structure of OTA, ADC/DAC circuit parameters), but needs to strengthen the description of the interaction between the control system and MEMS devices, as well as the system-level performance and application value. It is not necessary to overemphasize circuit indicators such as chip area and power consumption, but it is necessary to highlight the improvement of angular control precision and system stability brought by the PID control scheme.

# 4. General Submission Tips for Both Conferences

1. **Clarify the Implementation Route**: Clearly state whether the PID control is implemented in analog, digital, or mixed-signal mode, and explain the reasons for choosing this route (combined with the requirements of high-precision angular control, such as low latency, flexible parameter adjustment).

2. **Highlight Innovation and Practicality**: Avoid simply repeating the basic principle of PID control; focus on the customized optimization for MEMS high-precision angular control and the engineering feasibility of the design (supported by sufficient simulation and test data).

3. **Standardize the Comparison with SOTA Works**: In the Comparison Table, ensure that the selected SOTA works are recent (within 3-5 years) and from top conferences/journals. Clearly point out the advantages of your design (such as higher angular control precision, lower energy efficiency, or better environmental robustness).

4. **Pay Attention to Deadline**: The paper submission deadline is April 3, 2026. It is recommended to reserve sufficient time for revision, including optimizing the experimental data, improving the figures, and modifying the paper language (ensure that the English expression is accurate and professional, avoiding grammatical errors).
