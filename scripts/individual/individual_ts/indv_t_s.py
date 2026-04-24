import sys
import os

# Get the path to the 'sem8' directory (two levels up from this file)
root_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))

# Add it to the system path so Python can find your module
if root_path not in sys.path:
    sys.path.insert(0, root_path)

from math import ceil
import numpy as np
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI
from compare_refrigerants import make_saturation_data

TS_PLOT_DIR = './'

fluids = ['R1234yf', 'R1234ze(E)','R1234ze(Z)', 'R1233zd(E)', 'R1336mzz(E)', 'R1243zf', 'R134a']

Tmin_C = -75.0
Tmax_C = 250.0
nT = 181
T_C = np.linspace(Tmin_C, Tmax_C, nT)
T_K = T_C + 273.15
T_evap_range = np.linspace(-40, 10, 26)


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


# ---------- T-s diagrams with VCRS cycle overlay ----------
# Operating conditions for the cycle overlay
T_cond_C = 40.0
superheat_C = 5.0
subcool_C = 3.0
eta_isentropic = 0.75


def vcrs_ts_points(fluid, Te_C, Tc_C, sh_C=5.0, sc_C=3.0, eta_s=0.75):
    """
    Return T-s state points for a simple VCRS cycle.
    State 1: compressor inlet, superheated vapor
    State 2: compressor outlet
    State 3: condenser outlet, subcooled liquid
    State 4: after throttling valve
    """
    try:
        Te_K = Te_C + 273.15
        Tc_K = Tc_C + 273.15

        # Saturation pressures
        P_evap = PropsSI('P', 'T', Te_K, 'Q', 1, fluid)
        P_cond = PropsSI('P', 'T', Tc_K, 'Q', 0, fluid)

        # State 1, suction vapor
        T1_K = Te_C + sh_C + 273.15
        h1 = PropsSI('H', 'T', T1_K, 'P', P_evap, fluid)
        s1 = PropsSI('S', 'T', T1_K, 'P', P_evap, fluid)

        # Isentropic and actual compression
        h2s = PropsSI('H', 'P', P_cond, 'S', s1, fluid)
        h2 = h1 + (h2s - h1) / eta_s
        T2_K = PropsSI('T', 'P', P_cond, 'H', h2, fluid)
        s2 = PropsSI('S', 'P', P_cond, 'H', h2, fluid)

        # State 3, condenser outlet with subcooling
        T3_K = Tc_C - sc_C + 273.15
        h3 = PropsSI('H', 'T', T3_K, 'P', P_cond, fluid)
        s3 = PropsSI('S', 'T', T3_K, 'P', P_cond, fluid)

        # State 4, throttling valve outlet
        h4 = h3
        T4_K = PropsSI('T', 'P', P_evap, 'H', h4, fluid)
        s4 = PropsSI('S', 'P', P_evap, 'H', h4, fluid)

        return {
            '1': {'T': T1_K, 's': s1},
            '2': {'T': T2_K, 's': s2},
            '3': {'T': T3_K, 's': s3},
            '4': {'T': T4_K, 's': s4},
            'P_evap': P_evap,
            'P_cond': P_cond
        }

    except Exception as e:
        print(f"{fluid} failed: {e}")
        return None


nfluids = len(fluids)
ncols = 2
nrows = ceil(nfluids / ncols)

fig, axes = plt.subplots(nrows, ncols, figsize=(12, 4 * nrows), dpi=140)
axes = np.array(axes).reshape(-1)

for ax, fluid in zip(axes, fluids):
    data = sat_data[fluid]

    # Saturation dome
    valid_l = np.isfinite(data['s_l']) & np.isfinite(data['T_K'])
    valid_v = np.isfinite(data['s_v']) & np.isfinite(data['T_K'])

    if np.any(valid_l):
        ax.plot(data['s_l'][valid_l] * 1e-3,
                data['T_K'][valid_l] - 273.15,
                color='darkblue', lw=1.2)

    if np.any(valid_v):
        ax.plot(data['s_v'][valid_v] * 1e-3,
                data['T_K'][valid_v] - 273.15,
                color='darkblue', lw=1.2)
        
    # VCRS cycle points
    Te_plot_C = safe_evap_temp_C(fluid, T_evap_range[len(T_evap_range)//2], margin_C=2.0)
    pts = vcrs_ts_points(fluid, Te_plot_C, T_cond_C,
                     sh_C=superheat_C, sc_C=subcool_C, eta_s=eta_isentropic)
    if pts is None:
        ax.text(0.5, 0.5, 'VCRS cycle not valid at chosen state', transform=ax.transAxes, ha='center', va='center', fontsize=9)

    if pts is not None:
        s_cycle = np.array([pts['1']['s'], pts['2']['s'], pts['3']['s'], pts['4']['s'], pts['1']['s']]) * 1e-3
        T_cycle = np.array([pts['1']['T'], pts['2']['T'], pts['3']['T'], pts['4']['T'], pts['1']['T']]) - 273.15

        ax.plot(s_cycle, T_cycle, color='red', lw=1.4)

        for label in ['1', '2', '3', '4']:
            x = pts[label]['s'] * 1e-3
            y = pts[label]['T'] - 273.15
            ax.plot(x, y, marker='o', color='red', ms=3)
            ax.annotate(label, (x, y), textcoords='offset points', xytext=(4, 4), fontsize=10, color='red')

    ax.set_title(fluid, fontsize=11)
    ax.set_xlabel(r'$s$ (kJ/kg$\cdot$K)')
    ax.set_ylabel(r'$t$ ($^\circ$C)')
    ax.grid(True, ls=':')
    ax.tick_params(direction='in', top=True, right=True)

# Hide extra axes if number of fluids is odd
for ax in axes[nfluids:]:
    ax.axis('off')

plt.suptitle('Individual T-s saturation curves',fontsize=14)
plt.tight_layout()
ts_vcrs_png = os.path.join(TS_PLOT_DIR, 'ts_vcrs_comparison.png')
plt.savefig(ts_vcrs_png, bbox_inches='tight')
plt.show()

print("Saved:", ts_vcrs_png)


"""
# INTERACTIVE T-S :-
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

TS_PLOT_DIR = './'

# Operating conditions 
Tmin_C = -75.0
Tmax_C = 250.0
nT = 181
T_C = np.linspace(Tmin_C, Tmax_C, nT)
T_K_range = T_C + 273.15
T_cond_C = 40.0
superheat_C = 5.0
subcool_C = 3.0
eta_isentropic = 0.75
T_evap_range = np.linspace(-40, 10, 26)
fluids = ['R1234yf', 'R1234ze(E)','R1234ze(Z)', 'R1233zd(E)', 'R1336mzz(E)', 'R1243zf', 'R134a']

T_K_range = np.linspace(-75, 250, 181) + 273.15
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


def vcrs_ts_points(fluid, Te_C, Tc_C, sh_C=5.0, sc_C=3.0, eta_s=0.75):
    try:
        Te_K, Tc_K = Te_C + 273.15, Tc_C + 273.15
        P_evap = PropsSI('P', 'T', Te_K, 'Q', 1, fluid)
        P_cond = PropsSI('P', 'T', Tc_K, 'Q', 0, fluid)
        
        # State 1: Suction
        T1_K = Te_C + sh_C + 273.15
        h1 = PropsSI('H', 'T', T1_K, 'P', P_evap, fluid)
        s1 = PropsSI('S', 'T', T1_K, 'P', P_evap, fluid)
        
        # State 2: Discharge
        h2s = PropsSI('H', 'P', P_cond, 'S', s1, fluid)
        h2 = h1 + (h2s - h1) / eta_s
        T2_K = PropsSI('T', 'P', P_cond, 'H', h2, fluid)
        s2 = PropsSI('S', 'P', P_cond, 'H', h2, fluid)
        
        # State 3: Condenser exit
        T3_K = Tc_C - sc_C + 273.15
        s3 = PropsSI('S', 'T', T3_K, 'P', P_cond, fluid)
        
        # State 4: Expansion exit
        h4 = PropsSI('H', 'T', T3_K, 'P', P_cond, fluid) # isenthalpic
        s4 = PropsSI('S', 'P', P_evap, 'H', h4, fluid)
        T4_K = PropsSI('T', 'P', P_evap, 'H', h4, fluid)

        return {
            's': [s1, s2, s3, s4, s1], 
            'T': [T1_K, T2_K, T3_K, T4_K, T1_K]
        }
    except: return None

# Create Subplots
ncols = 2
nrows = ceil(len(fluids) / ncols)
fig = make_subplots(rows=nrows, cols=ncols, subplot_titles=fluids, vertical_spacing=0.08)

for i, fluid in enumerate(fluids):
    row, col = (i // ncols) + 1, (i % ncols) + 1
    data = make_saturation_data(fluid, T_K_range)
    
    # 1. Saturation Dome
    fig.add_trace(go.Scatter(x=data['s_l']*1e-3, y=data['T_K']-273.15, name=f'{fluid} Sat', 
                             line=dict(color='darkblue', width=1.5), showlegend=False), row=row, col=col)
    fig.add_trace(go.Scatter(x=data['s_v']*1e-3, y=data['T_K']-273.15, name=f'{fluid} Sat', 
                             line=dict(color='darkblue', width=1.5), showlegend=False), row=row, col=col)
    
    # 2. VCRS Cycle
    Te_plot_C = safe_evap_temp_C(fluid, T_evap_range[len(T_evap_range)//2])
    pts = vcrs_ts_points(fluid, Te_plot_C, T_cond_C, 
                         sh_C=superheat_C, sc_C=subcool_C, eta_s=eta_isentropic)
    if pts:
        s_vals = np.array(pts['s']) * 1e-3
        t_vals = np.array(pts['T']) - 273.15
        
        fig.add_trace(go.Scatter(x=s_vals, y=t_vals, name=f'{fluid} Cycle', 
                                 line=dict(color='red', width=2), mode='lines+markers+text',
                                 text=['1', '2', '3', '4', ''], textposition="top right",
                                 hovertemplate='s: %{x:.3f} kJ/kgK<br>T: %{y:.1f} °C'), row=row, col=col)

    fig.update_xaxes(title_text="s (kJ/kg·K)", row=row, col=col)
    fig.update_yaxes(title_text="T (°C)", row=row, col=col)

# Layout adjustments
fig.update_layout(height=450*nrows, width=1100, title_text="Interactive T-s Diagrams", 
                  template="plotly_white", showlegend=False)

# Save and Show
fig.write_html("ts_interactive.html")
fig.show()
"""