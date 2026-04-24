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



HS_PLOT_DIR = './'

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
    try:
        Tmin_C = PropsSI('Tmin', fluid) - 273.15
    except Exception:
        Tmin_C = -100.0
    return max(desired_Te_C, Tmin_C + margin_C)

sat_data = {fluid: make_saturation_data(fluid, T_K) for fluid in fluids}

def vcrs_hs_points(fluid, Te_C, Tc_C, sh_C=5.0, sc_C=3.0, eta_s=0.75):
    """Return h-s state points for a simple VCRS cycle."""
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
        s2 = PropsSI('S', 'P', P_cond, 'H', h2, fluid)
        
        # State 3: Condenser outlet
        T3_K = Tc_C - sc_C + 273.15
        h3 = PropsSI('H', 'T', T3_K, 'P', P_cond, fluid)
        s3 = PropsSI('S', 'T', T3_K, 'P', P_cond, fluid)

        # State 4: Expansion outlet
        h4 = h3
        s4 = PropsSI('S', 'P', P_evap, 'H', h4, fluid)

        return {
            '1': {'s': s1, 'h': h1},
            '2': {'s': s2, 'h': h2},
            '3': {'s': s3, 'h': h3},
            '4': {'s': s4, 'h': h4}
        }
    except Exception:
        return None

nfluids = len(fluids)
ncols = 2
nrows = ceil(nfluids / ncols)

fig, axes = plt.subplots(nrows, ncols, figsize=(12, 4 * nrows), dpi=140)
axes = np.array(axes).reshape(-1)

for ax, fluid in zip(axes, fluids):
    data = sat_data[fluid]
    
    # Plot Saturation Dome (h on y-axis, s on x-axis)
    ax.plot(data['s_l'] * 1e-3, data['h_l'] * 1e-3, color='darkblue', lw=1.2)
    ax.plot(data['s_v'] * 1e-3, data['h_v'] * 1e-3, color='darkblue', lw=1.2)
        
    Te_plot_C = safe_evap_temp_C(fluid, T_evap_range[len(T_evap_range)//2])
    pts = vcrs_hs_points(fluid, Te_plot_C, T_cond_C, sh_C=superheat_C, sc_C=subcool_C, eta_s=eta_isentropic)
    
    if pts:
        s_cycle = np.array([pts['1']['s'], pts['2']['s'], pts['3']['s'], pts['4']['s'], pts['1']['s']]) * 1e-3
        h_cycle = np.array([pts['1']['h'], pts['2']['h'], pts['3']['h'], pts['4']['h'], pts['1']['h']]) * 1e-3

        ax.plot(s_cycle, h_cycle, color='red', lw=1.4)

        for label in ['1', '2', '3', '4']:
            x, y = pts[label]['s'] * 1e-3, pts[label]['h'] * 1e-3
            ax.plot(x, y, marker='o', color='red', ms=3)
            ax.annotate(label, (x, y), textcoords='offset points', xytext=(4, 4), fontsize=10, color='red')

    ax.set_title(fluid, fontsize=11)
    ax.set_xlabel(r'$s$ (kJ/kg$\cdot$K)')
    ax.set_ylabel(r'$h$ (kJ/kg)')
    ax.grid(True, ls=':', alpha=0.5)
    ax.tick_params(direction='in', top=True, right=True)

for ax in axes[nfluids:]:
    ax.axis('off')

plt.suptitle('Individual h-s Saturation Curves with VCRS Cycle', fontsize=14)
plt.tight_layout()
hs_vcrs_png = os.path.join(HS_PLOT_DIR, 'hs_vcrs_comparison.png')
plt.savefig(hs_vcrs_png, bbox_inches='tight')
plt.show()

print("Saved:", hs_vcrs_png)



"""FOR INTERACTIVE GRAPHS :-
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



HS_PLOT_DIR = './'

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

    # Return a valid evaporator saturation temperature for the fluid.
    # Keeps Te above Tmin + margin.

    try:
        Tmin_C = PropsSI('Tmin', fluid) - 273.15
    except Exception:
        Tmin_C = -100.0
    return max(desired_Te_C, Tmin_C + margin_C)


def vcrs_hs_points(fluid, Te_C, Tc_C, sh_C=5.0, sc_C=3.0, eta_s=0.75):
    try:
        Te_K, Tc_K = Te_C + 273.15, Tc_C + 273.15
        P_evap = PropsSI('P', 'T', Te_K, 'Q', 1, fluid)
        P_cond = PropsSI('P', 'T', Tc_K, 'Q', 0, fluid)
        T1_K = Te_C + sh_C + 273.15
        h1 = PropsSI('H', 'T', T1_K, 'P', P_evap, fluid)
        s1 = PropsSI('S', 'T', T1_K, 'P', P_evap, fluid)
        h2s = PropsSI('H', 'P', P_cond, 'S', s1, fluid)
        h2 = h1 + (h2s - h1) / eta_s
        s2 = PropsSI('S', 'P', P_cond, 'H', h2, fluid)
        T3_K = Tc_C - sc_C + 273.15
        h3 = PropsSI('H', 'T', T3_K, 'P', P_cond, fluid)
        s3 = PropsSI('S', 'T', T3_K, 'P', P_cond, fluid)
        h4 = h3
        s4 = PropsSI('S', 'P', P_evap, 'H', h4, fluid)
        return {'s': [s1, s2, s3, s4, s1], 'h': [h1, h2, h3, h4, h1]}
    except: return None

# Create Subplots
ncols = 2
nrows = ceil(len(fluids) / ncols)
fig = make_subplots(rows=nrows, cols=ncols, subplot_titles=fluids, vertical_spacing=0.1)

for i, fluid in enumerate(fluids):
    row, col = (i // ncols) + 1, (i % ncols) + 1
    data = make_saturation_data(fluid, T_K)
    
    # 1. Plot Saturation Dome
    fig.add_trace(go.Scatter(x=data['s_l']*1e-3, y=data['h_l']*1e-3, name=f'{fluid} Sat', 
                             line=dict(color='blue', width=1), showlegend=False), row=row, col=col)
    fig.add_trace(go.Scatter(x=data['s_v']*1e-3, y=data['h_v']*1e-3, name=f'{fluid} Sat', 
                             line=dict(color='blue', width=1), showlegend=False), row=row, col=col)
    
    # 2. Plot VCRS Cycle
    Te_plot_C = safe_evap_temp_C(fluid, T_evap_range[len(T_evap_range)//2])
    pts = vcrs_ph_points(fluid, Te_plot_C, T_cond_C, 
                         sh_C=superheat_C, sc_C=subcool_C, eta_s=eta_isentropic)
    if pts:
        s_vals = np.array(pts['s']) * 1e-3
        h_vals = np.array(pts['h']) * 1e-3
        
        # Cycle Path
        fig.add_trace(go.Scatter(x=s_vals, y=h_vals, name=f'{fluid} Cycle', 
                                 line=dict(color='red', width=2), mode='lines+markers+text',
                                 text=['1', '2', '3', '4', ''], textposition="top right",
                                 hovertemplate='s: %{x:.3f} kJ/kgK<br>h: %{y:.1f} kJ/kg'), row=row, col=col)

    # Axis labels
    fig.update_xaxes(title_text="s (kJ/kg·K)", row=row, col=col)
    fig.update_yaxes(title_text="h (kJ/kg)", row=row, col=col)

# General Layout Update
fig.update_layout(height=400*nrows, width=1000, title_text="Interactive h-s Diagrams", 
                  template="plotly_white", showlegend=True)

# Save and Show
fig.write_html("hs_interactive.html")
fig.show()
"""