# Comparative Analysis of Low-GWP Refrigerants Using CoolProp

Refer Inference.txt for detailed analysis

## Overview

This project performs a **thermodynamic comparison of modern low-Global Warming Potential (GWP) refrigerants** using the **CoolProp thermophysical property library** in Python.
CoolProp is a thermophysical property database and wrappers for a selection of programming environments. It offers similar functionality to REFPROP, but CoolProp is open-source and free. 

The simulation models an **ideal vapor compression refrigeration cycle** and evaluates multiple **Hydrofluoroolefin (HFO)** refrigerants against the conventional **HFC refrigerant R134a**.

The goal is to determine which refrigerant provides the best performance based on:

* Cooling capacity
* Compressor power consumption
* Mass flow rate requirement
* Coefficient of Performance (COP)

The project also generates **thermodynamic diagrams and comparative plots** to visualize refrigerant behavior and performance.

---

# Objectives

1. Simulate a **vapor compression refrigeration cycle** using real thermodynamic properties.
2. Compare multiple **next-generation HFO refrigerants**.
3. Evaluate refrigerant performance using **CoolProp property calculations**.
4. Generate visual graphs for comparative analysis.
5. Rank refrigerants automatically based on performance metrics.

---

# Refrigerants Studied

### Hydrofluoroolefins (HFO)

* R1234yf
* R1234ze(E)
* R1234ze(Z)
* R1233zd(E)
* R1336mzz(E)
* R1243zf

### Control Refrigerant (HFC)

* R134a

These refrigerants are selected due to their **low Global Warming Potential (GWP)** and potential to replace conventional HFC refrigerants.

---

# Simulation Model

The project simulates an **Ideal Vapor Compression Cycle (VCC)** consisting of four main components:

1. **Evaporator**
2. **Compressor**
3. **Condenser**
4. **Expansion Valve**

The thermodynamic states are calculated using **CoolProp's `PropsSI()` function**, which provides accurate fluid properties.

### Assumptions

* Steady state operation
* Isentropic compressor efficiency = **0.75**
* Condenser temperature = **40°C**
* Evaporator temperature range = **−40°C to 10°C**
* Ideal expansion valve (isenthalpic process)

---

# Cycle Process

| State | Process                       |
| ----- | ----------------------------- |
| 1 → 2 | Isentropic Compression        |
| 2 → 3 | Heat Rejection in Condenser   |
| 3 → 4 | Expansion through Valve       |
| 4 → 1 | Heat Absorption in Evaporator |

The performance is evaluated using:

**Cooling capacity**

Q_in = h1 − h4

**Compressor work**

W_comp = h2 − h1

**Coefficient of Performance**

COP = Q_in / W_comp

---

# Generated Plots

The project automatically generates the following graphs:

### Thermodynamic Property Plots

* Saturation **Pressure–Temperature (P–T) curves**
* **Enthalpy vs Temperature** for saturated liquid and vapor
* **Temperature–Entropy (T–s) diagrams**

### Performance Comparison Plots

* **COP vs Evaporator Temperature**
* **Cooling Capacity vs Evaporator Temperature**
* **Mass Flow Rate vs Evaporator Temperature**
* **Compressor Power vs Evaporator Temperature**


---

# Automatic Refrigerant Ranking

The project implements an **automatic ranking system** based on:

* Cooling capacity
* Compressor power requirement
* Mass flow rate

Each metric is normalized and combined into a **performance score**.

Final ranking example:

| Rank | Refrigerant |
| ---- | ----------- |
| 1    | R1234ze(Z)  |
| 2    | R1233zd(E)  |
| 3    | R1243zf     |
| 4    | R1234ze(E)  |
| 5    | R1336mzz(E) |
| 6    | R1234yf     |

This ranking indicates that **R1234ze(Z)** provides the best thermodynamic performance under the selected conditions.

---

# Project Structure

```
project/
│
├── scripts
|   ├── compare_refrigerants.py
|   ├── data/
│       ├──saturation_plots
|          ├── saturation_PT_logP.png
|          ├── saturation_h_T.png
│          ├── saturation_Ts.png
│          └── cop_vs_Tevap.png
│       ├──plots
|          ├── cop_vs_Tevap.png
│          ├── power_vs_Tevap.png
|          ├── hfo_cooling_vs_Te.png
│          ├── massflow_vs_Tevap.png
|          ├── hfo_comparative_raw.csv
|          ├── hfo_comparative_summary.csv
|          └── hfo_comparative_summary_extended.csv
├── Inference.txt
└── README.md
```

---

# Installation

### 1. Install Python

Python **3.12 or later** is recommended.

### 2. Install required libraries

```bash
pip install -r requirements.txt
```

---

# Running the Simulation

Run the main script:

```bash
python compare_refrigerants.py
```

The script will:

1. Calculate thermodynamic properties
2. Simulate the refrigeration cycle
3. Generate graphs
4. Export data to CSV files
5. Produce refrigerant performance rankings

---

# Example CoolProp Property Call

CoolProp calculates fluid properties using the `PropsSI()` function.

Example:

```python
from CoolProp.CoolProp import PropsSI

P = PropsSI('P','T',300,'Q',0,'R134a')
```

This returns the **saturation pressure of R134a at 300 K**.

---

# Key Results

Based on the combined thermodynamic scoring of cooling capacity, compressor power and mass flow rate, 
R1234ze(Z) emerged as the best performing refrigerant for the simulated vapor compression cycle, 
followed by R1233zd(E). 
R1234yf showed the lowest thermodynamic performance under the selected operating conditions.

---



---

# References

* CoolProp Documentation
  [https://coolprop.org/](https://coolprop.org/)


* Various research papers on **HFO refrigerants and vapor compression cycles**

---

# Author
Project developed as part of : Thermodynamic analysis of HydroFlouro Olefin(HFO) refrigerants using computational tools.


compare_refrigerants.py
Requires: CoolProp, numpy, pandas, matplotlib
refer docs : https://coolprop.org/coolprop/HighLevelAPI.html#propssi-function



HFO : R1132a (1,1-difluoroethylene) # NOT AVAILABLE