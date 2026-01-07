# 🚗 Active Suspension Control - Quarter-Car Model

<div align="center">

![MATLAB](https://img.shields.io/badge/MATLAB-R2024a-orange?style=for-the-badge&logo=mathworks)
![Simulink](https://img.shields.io/badge/Simulink-Modeling-blue?style=for-the-badge)
![License](https://img.shields.io/badge/License-Academic-green?style=for-the-badge)

**Complements of Automatic Controls - Final Project**

*Università degli Studi di Napoli Federico II*

</div>

---

## 📋 Overview

Design and comparison of different control strategies for an **active suspension system** using the Quarter-Car model. The goal is to control passenger comfort by reducing body accelerations while maintaining tire-road contact.

---

## 🎯 Control Objectives

- Zero steady-state error on body acceleration
- Suspension deflection constraint: |zₛ - zᵤ| < 0.08 m
- 1% settling time: tₐ ≤ 3 s
- Overshoot: S% ≤ 10%
- Peak acceleration < 3 m/s²

---

## 🔧 Implemented Controllers

| Controller | Description |
|------------|-------------|
| **State Feedback** | Pole placement with dominant poles method |
| **State Feedback + Integral Action** | Augmented system for zero steady-state error |
| **Observer-Based Control** | Luenberger observer for state estimation |
| **LQ Control** | Linear Quadratic Regulator with output weighting |
| **LQG Control** | LQ + Kalman Filter for noisy measurements |
| **H∞ Control** | Mixed Sensitivity Design for robust control |

---

## 📁 Repository Contents

| File | Description |
|------|-------------|
| `CC_Project_Report.pdf` | Complete project report |
| `LQ_Controller.slx` | LQ controller Simulink model |
| `LQ_Controller_schema.pdf` | LQ controller block diagram |
| `LQ_Kalman_LQG.slx` | LQG controller Simulink model |
| `LQ_Kalman_LQG_schema.pdf` | LQG controller block diagram |
| `LQ_Kalman_Filter_and_LQG.mlx` | MATLAB Live Script for LQ/LQG design |
| `Observer.slx` | Observer-based controller Simulink model |
| `Osservatore_sospensioni_schema.pdf` | Observer block diagram |
| `State_feedback_integrator.slx` | State feedback with integral action model |
| `Statefeedback_and_Observer.mlx` | MATLAB Live Script for state feedback design |
| `Sospensione_state_feedback_and_integrator.slx` | Complete state feedback model |
| `Hinfinity_controller.slx` | H∞ controller Simulink model |
| `Hinf_sim_schema.pdf` | H∞ controller block diagram |
| `Mixed_Sensitivity_Design_e_Hinfinity.mlx` | MATLAB Live Script for H∞ design |
| `Test_prestazioni_controllo_H_infinity.m` | H∞ performance testing script |

---

## 🛠️ Requirements

- MATLAB R2020a or later
- Simulink
- Control System Toolbox
- Robust Control Toolbox

---

## 🚀 Usage

1. Open MATLAB and navigate to the repository folder
2. Run the desired `.mlx` Live Script for controller design
3. Open the corresponding `.slx` Simulink model
4. Run the simulation to analyze closed-loop performance

---

## 👤 Author

**Emanuela Varone** - 

**Course**: Complements of Automatic Controls  
**Institution**: Università degli Studi di Napoli Federico II  
**Department**: DIETI

---

## 📄 License

This project is developed for academic purposes.
