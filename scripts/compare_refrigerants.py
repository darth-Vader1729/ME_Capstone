# compare_refrigerants.py
# Requires: CoolProp, numpy, pandas, matplotlib
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI

# ---------- USER SETTINGS ----------
"""
Hydro-flouro olefins used:-

R1234yf (2,3,3,3-tetrafluoropropene)
R1234ze(E) (trans-1,3,3,3-tetrafluoropropene)
R1234ze(Z) (cis-1,3,3,3-tetrafluoropropene)
R1233zd(E) (trans-1-chloro-3,3,3-trifluoropropene) - A commonly used HCFO
R1336mzz(E) (trans-1,1,1,4,4,4-hexafluoro-2-butene)
R1336mzz(Z) (cis-1,1,1,4,4,4-hexafluoro-2-butene) # NOT AVAILABLE
R1243zf (3,3,3-trifluoropropene)
R1224yd(Z) (2-chloro-3,3,3-trifluoropropene)

HydroFlouro carbon refrigerant as control :- R134a

HFO : R1132a (1,1-difluoroethylene) # NOT AVAILABLE
"""
HFOs = ['R1234yf', 'R1234ze(E)','R1234ze(Z)', 'R1233zd(E)', 'R1336mzz(E)', 'R1243zf']   # HFO candidates
HFCs = ['R134a'] # HFC candidates; others 'R404A', 'R410A'

# for saturation curves
Tmin_C = -75.0
Tmax_C = 250.0
nT = 181 # number of points to sample between Tmin and Tmax

T_cond_C = 40.0
T_evap_range = np.linspace(-40, 10, 26) 
# linspace : Returns evenly spaced numbers over a specified interval.
eta_isentropic = 0.75 # compressor isen efficiency

# Reference cooling load for mass flow calculations (user adjustable)
Q_ref_kW = 5.0  # kW of cooling to be provided by the evaporator


OUTDIR = './data'
PLOT_DIR = os.path.join(OUTDIR, 'plots')
os.makedirs(PLOT_DIR, exist_ok=True)
os.makedirs(OUTDIR, exist_ok=True)

# Pull-down simulation parameters (tune these to match experiments)
room_mass_kg = 1000.0   # effective mass of air+contents, kg
c_p_room = 1005.0       # J/kg-K
UA_room = 50.0          # overall U*A, W/K
T_ambient = 30.0        # ambient boundary, degC
T_initial = 30.0        # initial room temp, degC
T_setpoint = 5.0        # target pull-down temperature, degC
dt = 1.0                # s
t_max = 3600.0          # s


# ---------- HELPERS ----------
# refer docs : https://coolprop.org/coolprop/HighLevelAPI.html#propssi-function
# first parameter is the output property that will be returned from PropsSI
def is_fluid_available(fluid):
    try:
        PropsSI('P','T',273.15,'Q',0,fluid)
        return True
    except Exception:
        return False

# wrapper around PropsSI : to catch any exception
def safe_propsSI(output, input1_name, input1_val, input2_name, input2_val, fluid):
    try:
        return PropsSI(output, input1_name, input1_val, input2_name, input2_val, fluid)
    except Exception as e:
        # return NaN and optionally print a small debug message
        # (Don't spam; only called during array creation so this is OK)
        # print(f"Warning: PropsSI failed for {fluid} ({output}) at {input1_name}={input1_val}, {input2_name}={input2_val}: {e}")
        return np.nan

# ---------- CHECK FLUIDS ----------
all_fluids = HFOs + HFCs
available = {f: is_fluid_available(f) for f in all_fluids}
print("Availability check (True means CoolProp recognized the fluid):")
for f, ok in available.items():
    print(f"  {f:12s}: {ok}")
print()

HFOs_avail = [f for f in HFOs if available.get(f, False)]
HFCs_avail = [f for f in HFCs if available.get(f, False)]
if not (HFOs_avail or HFCs_avail):
    raise RuntimeError("No fluids available in CoolProp from the provided lists.")

# ---------- SATURATION CURVES ----------
T_C = np.linspace(Tmin_C, Tmax_C, nT)
T_K = T_C + 273.15

# ref docs : https://coolprop.org/coolprop/HighLevelAPI.html#vapor-liquid-and-saturation-states
# they take 'x' vapour mass fraction as 'Q'
def make_saturation_data(fluid, T_K_array):
    # compute arrays (list comprehensions keep safe_propsSI behavior)
    # calculating sat. values
    P_sat = np.array([safe_propsSI('P','T',T,'Q',0,fluid) for T in T_K_array])
    h_l = np.array([safe_propsSI('H','T',T,'Q',0,fluid) for T in T_K_array])
    h_v = np.array([safe_propsSI('H','T',T,'Q',1,fluid) for T in T_K_array])
    s_l = np.array([safe_propsSI('S','T',T,'Q',0,fluid) for T in T_K_array])
    s_v = np.array([safe_propsSI('S','T',T,'Q',1,fluid) for T in T_K_array])

    # ref docs : https://coolprop.org/coolprop/HighLevelAPI.html#trivial-inputs
    # try to get critical point (may fail for some named fluids)
    try:
        T_crit = PropsSI('Tcrit', fluid)
        P_crit = PropsSI('Pcrit', fluid)
    except Exception:
        T_crit = np.nan
        P_crit = np.nan
    return dict(T_K=T_K_array, P=P_sat, h_l=h_l, h_v=h_v, s_l=s_l, s_v=s_v,
                T_crit=T_crit, P_crit=P_crit)

sat_data = {}
for fluid in HFOs_avail + HFCs_avail:
    sat_data[fluid] = make_saturation_data(fluid, T_K)

# ---------- P vs T plot (log pressure) ----------
plt.figure(figsize=(10,5), dpi=120) # 10*5 inches
for fluid, data in sat_data.items():
    P_bar = data['P'] / 1e5  # Pa -> bar
    valid = np.isfinite(P_bar) & (P_bar > 0)
    if not np.any(valid):
        print(f"Skipping {fluid} on P-T plot: no valid saturation pressures.")
        continue
    # x - linear, y - log, marker
    plt.semilogy(T_C[valid], P_bar[valid], label=fluid)
    # mark critical point if available
    if np.isfinite(data['T_crit']) and np.isfinite(data['P_crit']):
        Tc_C = data['T_crit'] - 273.15
        Pc_bar = data['P_crit'] / 1e5
        plt.plot([Tc_C], [Pc_bar], marker='o', markersize=4)
        plt.text(Tc_C + 1.0, Pc_bar*1.05, 'Tc', fontsize=8)

plt.xlabel('Temperature (°C)')
plt.ylabel('Saturation pressure (bar, log scale)')
plt.title('Saturation P-T curves (log P)')
plt.grid(True, which='both', ls=':')
plt.legend(ncol=2, fontsize='small')
plt.tight_layout()
plt.savefig(os.path.join(PLOT_DIR, 'saturation_PT_logP.png'))
plt.show()

# ---------- Enthalpy vs T (saturated liquid & vapor) ----------
plt.figure(figsize=(10,5), dpi=120)
for fluid, data in sat_data.items():
    valid_liq = np.isfinite(data['h_l'])
    valid_vap = np.isfinite(data['h_v'])
    if np.any(valid_vap):
        plt.plot(T_C[valid_vap], data['h_v'][valid_vap]/1000.0, '--', label=f'{fluid} vap')
    if np.any(valid_liq):
        plt.plot(T_C[valid_liq], data['h_l'][valid_liq]/1000.0, '-', label=f'{fluid} liq')

plt.xlabel('Temperature (°C)')
plt.ylabel('Specific enthalpy (kJ/kg)')
plt.title('Saturated liquid & vapor enthalpies vs T')
plt.grid(True, ls=':')
# reduce legend clutter: show only unique labels in order
handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles))
plt.legend(by_label.values(), by_label.keys(), ncol=2, fontsize='small')
plt.tight_layout()
plt.savefig(os.path.join(PLOT_DIR, 'saturation_h_T.png'))
plt.show()

# ---------- T-s diagrams (saturation curves) ----------
plt.figure(figsize=(8,6), dpi=120)
for fluid, data in sat_data.items():
    valid_l = np.isfinite(data['s_l']) & np.isfinite(data['T_K'])
    valid_v = np.isfinite(data['s_v']) & np.isfinite(data['T_K'])
    if np.any(valid_l):
        plt.plot(data['s_l'][valid_l]*1e-3, data['T_K'][valid_l]-273.15, '-', label=f'{fluid} liq')
    if np.any(valid_v):
        plt.plot(data['s_v'][valid_v]*1e-3, data['T_K'][valid_v]-273.15, '--', label=f'{fluid} vap')

plt.xlabel('Specific entropy (kJ/kg·K)')
plt.ylabel('Temperature (°C)')
plt.title('T-s saturation curves')
plt.grid(True, ls=':')
handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles))
plt.legend(by_label.values(), by_label.keys(), ncol=2, fontsize='small')
plt.tight_layout()
plt.savefig(os.path.join(PLOT_DIR, 'saturation_Ts.png'))
plt.show()

# ---------- SIMPLE VAPOR-COMPRESSION CYCLE COP FUNCTION ----------
def ideal_vcc_cop(fluid, T_evap_C, T_cond_C, eta_isentropic=0.75):
    try:
        T_evap = T_evap_C + 273.15
        T_cond = T_cond_C + 273.15
        P_evap = PropsSI('P','T',T_evap,'Q',1,fluid)
        P_cond = PropsSI('P','T',T_cond,'Q',0,fluid)
        h1 = PropsSI('H','P',P_evap,'Q',1,fluid)
        s1 = PropsSI('S','P',P_evap,'Q',1,fluid)
        h2s = PropsSI('H','P',P_cond,'S',s1,fluid)
        # actual compressor work (isentropic eff)
        h2 = h1 + (h2s - h1) / eta_isentropic
        h3 = PropsSI('H','P',P_cond,'Q',0,fluid)
        h4 = h3 # ideal throttle
        q_in = h1 - h4
        w_comp = h2 - h1
        if w_comp <= 0 or not np.isfinite(q_in) or not np.isfinite(w_comp):
            return np.nan
        return q_in / w_comp
    except Exception:
        return np.nan

# ---------- Sweep COP vs Evaporator temperature ----------
plt.figure(figsize=(10,5), dpi=120)
for fluid in HFOs_avail + HFCs_avail:
    cops = [ideal_vcc_cop(fluid, Te, T_cond_C, eta_isentropic) for Te in T_evap_range]
    cops = np.array(cops)
    valid = np.isfinite(cops)
    if not np.any(valid):
        print(f"No valid COP points for {fluid} at T_cond={T_cond_C} °C.")
        continue
    plt.plot(T_evap_range[valid], cops[valid], label=fluid)

plt.xlabel('Evaporator Temperature (°C)')
plt.ylabel('Ideal VCC COP (dimensionless)')
plt.title(f'COP vs Evaporator Temp (T_cond = {T_cond_C} °C, η_isentropic={eta_isentropic})')
plt.grid(True, ls=':')
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(PLOT_DIR, 'cop_vs_Tevap.png'))
plt.show()

# ---------- Export CSVs ----------
for fluid, data in sat_data.items():
    df = pd.DataFrame({
        'T_C': T_C,
        'T_K': data['T_K'],
        'P_Pa': data['P'],
        'h_liq_Jkg': data['h_l'],
        'h_vap_Jkg': data['h_v'],
        's_liq_JkgK': data['s_l'],
        's_vap_JkgK': data['s_v'],
    })
    fname = os.path.join(OUTDIR, f'saturation_{fluid.replace("/","_")}.csv')
    df.to_csv(fname, index=False)
    print("Wrote:", fname)

print("Plots saved to:", PLOT_DIR)
print("Done.")
