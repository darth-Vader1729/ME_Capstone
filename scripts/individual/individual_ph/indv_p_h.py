import sys
import os

# Get the path to the 'sem8' directory (two levels up from this file)
root_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))

# Add it to the system path so Python can find your module
if root_path not in sys.path:
    sys.path.insert(0, root_path)

from math import ceil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI
from compare_refrigerants import make_saturation_data


PH_PLOT_DIR = './'

# Operating conditions 
Tmin_C = -75.0
Tmax_C = 250.0
nT = 181
T_C = np.linspace(Tmin_C, Tmax_C, nT)
T_K = T_C + 273.15
T_cond_C = 40.0
superheat_C = 5.0
subcool_C = 3.0
eta_isentropic = 0.75
T_evap_range = np.linspace(-40, 10, 26)
fluids = ['R1234yf', 'R1234ze(E)','R1234ze(Z)', 'R1233zd(E)', 'R1336mzz(E)', 'R1243zf', 'R134a']

def safe_evap_temp_C(fluid, desired_Te_C, margin_C=2.0):
    """
    Return a valid evaporator saturation temperature for the fluid.
    Keeps Te above Tmin + margin.
    """
    try:
        Tmin_C = PropsSI('Tmin', fluid) - 273.15
    except Exception:
        Tmin_C = -100.0
    return max(desired_Te_C, Tmin_C + margin_C)


sat_data = {}
for fluid in fluids:
    sat_data[fluid] = make_saturation_data(fluid, T_K)
# print(sat_data)

def vcrs_ph_points(fluid, Te_C, Tc_C, sh_C=5.0, sc_C=3.0, eta_s=0.75):
    """
    Return P-h state points for a simple VCRS cycle.
    """
    try:
        Te_K = Te_C + 273.15
        Tc_K = Tc_C + 273.15

        # Saturation pressures
        P_evap = PropsSI('P', 'T', Te_K, 'Q', 1, fluid)
        P_cond = PropsSI('P', 'T', Tc_K, 'Q', 0, fluid)

        # State 1: Suction vapor (Superheated)
        T1_K = Te_C + sh_C + 273.15
        h1 = PropsSI('H', 'T', T1_K, 'P', P_evap, fluid)
        s1 = PropsSI('S', 'T', T1_K, 'P', P_evap, fluid)

        # State 2: Discharge (Actual compression)
        h2s = PropsSI('H', 'P', P_cond, 'S', s1, fluid)
        h2 = h1 + (h2s - h1) / eta_s
        
        # State 3: Condenser outlet (Subcooled)
        T3_K = Tc_C - sc_C + 273.15
        h3 = PropsSI('H', 'T', T3_K, 'P', P_cond, fluid)

        # State 4: Expansion outlet (Isenthalpic)
        h4 = h3

        return {
            '1': {'h': h1, 'P': P_evap},
            '2': {'h': h2, 'P': P_cond},
            '3': {'h': h3, 'P': P_cond},
            '4': {'h': h4, 'P': P_evap}
        }
    except Exception as e:
        return None

nfluids = len(fluids)
ncols = 2
nrows = ceil(nfluids / ncols)

fig, axes = plt.subplots(nrows, ncols, figsize=(12, 4 * nrows), dpi=140)
axes = np.array(axes).reshape(-1)

for ax, fluid in zip(axes, fluids):
    # Use saturation data from your make_saturation_data function
    # Assuming it provides 'h_l', 'h_v', and 'p' (in Pa)
    data = sat_data[fluid]

    # Plot Saturation Dome
    # Note: h is converted to kJ/kg, P to bar or kPa (using bar here for standard P-h)
    ax.plot(data['h_l'] * 1e-3, data['P'] * 1e-5, color='darkblue', lw=1.2)
    ax.plot(data['h_v'] * 1e-3, data['P'] * 1e-5, color='darkblue', lw=1.2)
        
    # VCRS cycle points
    Te_plot_C = safe_evap_temp_C(fluid, T_evap_range[len(T_evap_range)//2])
    pts = vcrs_ph_points(fluid, Te_plot_C, T_cond_C, 
                         sh_C=superheat_C, sc_C=subcool_C, eta_s=eta_isentropic)
    
    if pts:
        # Cycle sequence: 1 -> 2 -> 3 -> 4 -> 1
        h_cycle = np.array([pts['1']['h'], pts['2']['h'], pts['3']['h'], pts['4']['h'], pts['1']['h']]) * 1e-3
        p_cycle = np.array([pts['1']['P'], pts['2']['P'], pts['3']['P'], pts['4']['P'], pts['1']['P']]) * 1e-5

        ax.plot(h_cycle, p_cycle, color='red', lw=1.4)

        for label in ['1', '2', '3', '4']:
            x, y = pts[label]['h'] * 1e-3, pts[label]['P'] * 1e-5
            ax.plot(x, y, marker='o', color='red', ms=3)
            ax.annotate(label, (x, y), textcoords='offset points', xytext=(4, 4), fontsize=10, color='red')

    ax.set_title(fluid, fontsize=11)
    ax.set_xlabel(r'$h$ (kJ/kg)')
    ax.set_ylabel(r'$P$ (bar)')
    ax.set_yscale('log') # Pressure is typically logarithmic in P-h charts
    ax.grid(True, which='both', ls=':', alpha=0.5)
    ax.tick_params(direction='in', top=True, right=True)

# Clean up
for ax in axes[nfluids:]:
    ax.axis('off')

plt.suptitle('Individual P-h Saturation Curves with VCRS Cycle', fontsize=14)
plt.tight_layout()
ph_vcrs_png = os.path.join(PH_PLOT_DIR, 'ph_vcrs_comparison.png')
plt.savefig(ph_vcrs_png, bbox_inches='tight')
plt.show()

print("Saved:", ph_vcrs_png)



"""
# INTERACTIVE PH:-
import sys
import os

# Get the path to the 'sem8' directory (two levels up from this file)
root_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))

# Add it to the system path so Python can find your module
if root_path not in sys.path:
    sys.path.insert(0, root_path)


import plotly.graph_objects as go
from plotly.subplots import make_subplots
from math import ceil
import numpy as np
from CoolProp.CoolProp import PropsSI
from compare_refrigerants import make_saturation_data



PH_PLOT_DIR = './'

# Operating conditions 
Tmin_C = -75.0
Tmax_C = 250.0
nT = 181
T_C = np.linspace(Tmin_C, Tmax_C, nT)
T_K = T_C + 273.15
T_cond_C = 40.0
superheat_C = 5.0
subcool_C = 3.0
eta_isentropic = 0.75
T_evap_range = np.linspace(-40, 10, 26)
fluids = ['R1234yf', 'R1234ze(E)','R1234ze(Z)', 'R1233zd(E)', 'R1336mzz(E)', 'R1243zf', 'R134a']

T_K = np.linspace(-75, 250, 181) + 273.15
T_cond_C = 40.0
sh_C = 5.0
sc_C = 3.0
eta_s = 0.75
fluids = ['R1234yf', 'R1234ze(E)','R1234ze(Z)', 'R1233zd(E)', 'R1336mzz(E)', 'R1243zf', 'R134a']


def safe_evap_temp_C(fluid, desired_Te_C, margin_C=2.0):

    # Return a valid evaporator saturation temperature for the fluid.
    # Keeps Te above Tmin + margin.

    try:
        Tmin_C = PropsSI('Tmin', fluid) - 273.15
    except Exception:
        Tmin_C = -100.0
    return max(desired_Te_C, Tmin_C + margin_C)


def vcrs_ph_points(fluid, Te_C, Tc_C, sh_C=5.0, sc_C=3.0, eta_s=0.75):
    try:
        Te_K, Tc_K = Te_C + 273.15, Tc_C + 273.15
        P_evap = PropsSI('P', 'T', Te_K, 'Q', 1, fluid)
        P_cond = PropsSI('P', 'T', Tc_K, 'Q', 0, fluid)
        T1_K = Te_C + sh_C + 273.15
        h1 = PropsSI('H', 'T', T1_K, 'P', P_evap, fluid)
        s1 = PropsSI('S', 'T', T1_K, 'P', P_evap, fluid)
        h2s = PropsSI('H', 'P', P_cond, 'S', s1, fluid)
        h2 = h1 + (h2s - h1) / eta_s
        h3 = PropsSI('H', 'T', Tc_C - sc_C + 273.15, 'P', P_cond, fluid)
        return {'h': [h1, h2, h3, h3, h1], 'P': [P_evap, P_cond, P_cond, P_evap, P_evap]}
    except Exception as e: 
        print(e)
        return None

# Create Subplots
ncols = 2
nrows = ceil(len(fluids) / ncols)
fig = make_subplots(rows=nrows, cols=ncols, subplot_titles=fluids, vertical_spacing=0.08)

for i, fluid in enumerate(fluids):
    row, col = (i // ncols) + 1, (i % ncols) + 1
    data = make_saturation_data(fluid, T_K)
    
    # 1. Saturation Dome
    fig.add_trace(go.Scatter(x=data['h_l']*1e-3, y=data['P']*1e-5, name=f'{fluid} Sat', 
                             line=dict(color='darkblue', width=1.5), showlegend=False), row=row, col=col)
    fig.add_trace(go.Scatter(x=data['h_v']*1e-3, y=data['P']*1e-5, name=f'{fluid} Sat', 
                             line=dict(color='darkblue', width=1.5), showlegend=False), row=row, col=col)
    
    # 2. VCRS Cycle
    Te_plot_C = safe_evap_temp_C(fluid, T_evap_range[len(T_evap_range)//2])
    pts = vcrs_ph_points(fluid, Te_plot_C, T_cond_C, 
                         sh_C=superheat_C, sc_C=subcool_C, eta_s=eta_isentropic)
    if pts:
        h_vals = np.array(pts['h']) * 1e-3
        p_vals = np.array(pts['P']) * 1e-5
        
        fig.add_trace(go.Scatter(x=h_vals, y=p_vals, name=f'{fluid} Cycle', 
                                 line=dict(color='red', width=2), mode='lines+markers+text',
                                 text=['1', '2', '3', '4', ''], textposition="top center",
                                 hovertemplate='h: %{x:.1f} kJ/kg<br>P: %{y:.2f} bar'), row=row, col=col)

    # Logarithmic scale for Pressure and axis titles
    fig.update_yaxes(type="log", title_text="P (bar)", row=row, col=col)
    fig.update_xaxes(title_text="h (kJ/kg)", row=row, col=col)

# Layout adjustments
fig.update_layout(height=450*nrows, width=1100, title_text="Interactive P-h Diagrams", 
                  template="plotly_white", showlegend=False)

# Save as HTML and Open
fig.write_html("ph_interactive.html")
fig.show()
"""