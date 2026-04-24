# compare_refrigerants.py
# Requires: CoolProp, numpy, pandas, matplotlib
# refer docs : https://coolprop.org/coolprop/HighLevelAPI.html#propssi-function
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

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI


# ---------- USER SETTINGS ----------
HFOs = ['R1234yf', 'R1234ze(E)','R1234ze(Z)', 'R1233zd(E)', 'R1336mzz(E)', 'R1243zf']   # HFO candidates
HFCs = ['R134a'] # HFC candidates; others 'R404A', 'R410A'

# for saturation curves
Tmin_C = -75.0
Tmax_C = 250.0
nT = 181 # number of points to sample between Tmin and Tmax
T_cond_C = 40.0
T_evap_range = np.linspace(-40, 10, 26)  # linspace : Returns evenly spaced numbers(26 items) over a specified intervalSimulation[-40, 10].
eta_isentropic = 0.75 # compressor isen efficiency

# Reference cooling load for mass flow calculations (user adjustable)
Q_ref_kW = 5.0  # kW of cooling to be provided by the evaporator

# ----------- cooling capacity info -----------------

T_evap_range = np.linspace(-40, 10, 26)  # °C
T_cond_C = 40.0  # condenser temperature
superheat_C = 5.0
subcool_C = 3.0
eta_isentropic = 0.75
m_dot_fixed = 0.01  # kg/s fixed mass flow for cooling capacity plot

P_atm = 101325.0
T0 = 298.15

# ---------- HELPERS ----------
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


# they take (vapour mass fraction) 'x' as 'Q'
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
    return dict(T_K=T_K_array, P=P_sat, h_l=h_l, h_v=h_v, s_l=s_l, s_v=s_v, T_crit=T_crit, P_crit=P_crit)

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


def cycle_states(fluid, T_evap_C, T_cond_C, eta_isentropic=0.75):
    """Return primary state enthalpies and per-kg quantities for a simple 4-point VCC.
    All enthalpies in J/kg, temperatures in C.
    """
    try:
        Te = T_evap_C + 273.15
        Tc = T_cond_C + 273.15
        P_evap = PropsSI('P','T',Te,'Q',1,fluid)
        P_cond = PropsSI('P','T',Tc,'Q',0,fluid)
        # states
        h1 = PropsSI('H','P',P_evap,'Q',1,fluid)  # vap at evaporator exit
        s1 = PropsSI('S','P',P_evap,'Q',1,fluid)
        # isentropic compression to condenser pressure
        h2s = PropsSI('H','P',P_cond,'S',s1,fluid)
        # actual h2 with isentropic efficiency
        h2 = h1 + (h2s - h1) / eta_isentropic
        h3 = PropsSI('H','P',P_cond,'Q',0,fluid)  # condensate
        h4 = h3  # ideal throttle
        q_in_per_kg = h1 - h4  # refrigeration effect, J/kg
        w_comp_per_kg = h2 - h1  # compressor work per kg, J/kg
        return dict(
            h1=h1, h2=h2, h2s=h2s, h3=h3, h4=h4,
            P_evap=P_evap, P_cond=P_cond,
            q_in_per_kg=q_in_per_kg, w_comp_per_kg=w_comp_per_kg
        )
    except Exception as e:
        return None
    

# Approximate fallback property model (simple linear-ish)

def approx_props_simple(refrig, T_K, quality=0):
    base = {'R134a': 1.0, 'R1234yf': 0.95, 'R1234ze(E)':0.95, 'R1234ze(Z)':0.95, 'R1233zd(E)':0.9, 'R1336mzz(E)':0.9, 'R1243zf':0.95}
    b = base.get(refrig, 1.0)

    # provide different liquid/vapor baseline via quality flag
    if quality == 0:  # saturated liquid
        h = (100e3 * b) + (T_K - 273.15) * 50.0   # J/kg arbitrary baseline
        s = 800.0 * b + (T_K - 273.15) * 0.15

    else:  # saturated vapor
        h = (400e3 * b) + (T_K - 273.15) * 120.0  # J/kg arbitrary baseline
        s = 1200.0 * b + (T_K - 273.15) * 0.25

    p = P_atm
    return h, s, p


def get_props(refrig, T_C, quality=1):
    T_K = T_C + 273.15
    if PropsSI is not None:
        try:
            h = PropsSI('H','T',T_K,'Q',quality,refrig)   # J/kg
            s = PropsSI('S','T',T_K,'Q',quality,refrig)
            p = PropsSI('P','T',T_K,'Q',quality,refrig)
            return h, s, p

        except Exception:
            return approx_props_simple(refrig, T_K, quality)

    else:
        return approx_props_simple(refrig, T_K, quality)


def compute_cycle_metrics(refrig, Te_C, Tc_C, sh_C=superheat_C, sc_C=subcool_C, eta_s=eta_isentropic):

    # state 1: suction vapor at Te + superheat
    T1 = Te_C + sh_C
    h1, s1, p1 = get_props(refrig, T1, quality=1)

    # state 3: liquid after condenser and subcool
    T3 = Tc_C - sc_C
    h3, s3, p3 = get_props(refrig, T3, quality=0)

    # isentropic outlet enthalpy h2s approximate

    if PropsSI is not None:
        try:
            Te_K = Te_C + 273.15
            Tc_K = Tc_C + 273.15

            # get pressures
            p_evap = PropsSI('P','T',Te_K,'Q',0,refrig)
            p_cond = PropsSI('P','T',Tc_K,'Q',1,refrig)
            h2s = PropsSI('H','P',p_cond,'S',s1,refrig)

        except Exception:
            h2s = h1 * 1.5

    else:
        h2s = h1 * 1.5

    # actual h2
    h2 = h1 + (h2s - h1) / eta_s

    # q_e and compressor work per kg
    q_e = h1 - h3
    w_comp = h2 - h1
    COP = q_e / w_comp if (w_comp != 0 and np.isfinite(q_e) and np.isfinite(w_comp)) else np.nan

    return dict(COP=COP, q_e_Jkg=q_e, w_comp_Jkg=w_comp, h1=h1, h2=h2, h3=h3)

def normalize(series, higher_is_better=True):
    vals = series.to_numpy(dtype=float)
    # handle all-nan or constant series
    valid = np.isfinite(vals)
    if not np.any(valid):
        return np.full_like(vals, np.nan, dtype=float)
    vmin = np.nanmin(vals[valid])
    vmax = np.nanmax(vals[valid])
    if np.isclose(vmax, vmin):
        # no variation. give neutral score 0.5 for finite entries
        out = np.full_like(vals, 0.5, dtype=float)
        out[~valid] = np.nan
        return out
    norm = (vals - vmin) / (vmax - vmin)
    if higher_is_better:
        return norm
    else:
        return 1.0 - norm
    

if __name__ == "__main__":
    
    # plots in PLOT_DIR are -> cop/T_evap, m_dot/T_evap, power/T_evap, hfo_cooling/T_evap
    # saturation plots in -> SAT_PLOT_DIR

    OUTDIR = './data'
    PLOT_DIR = os.path.join(OUTDIR, 'plots')
    SAT_PLOT_DIR = os.path.join(OUTDIR, 'saturation_plots')

    os.makedirs(OUTDIR, exist_ok=True)
    os.makedirs(PLOT_DIR, exist_ok=True)
    os.makedirs(SAT_PLOT_DIR, exist_ok=True)
    # ---------- CHECK FLUIDS ----------
    HFOs_avail = []
    HFCs_avail = []

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
    # ref docs : https://coolprop.org/coolprop/HighLevelAPI.html#vapor-liquid-and-saturation-states

    T_C = np.linspace(Tmin_C, Tmax_C, nT)
    T_K = T_C + 273.15



    sat_data = {}
    for fluid in HFOs_avail + HFCs_avail:
        sat_data[fluid] = make_saturation_data(fluid, T_K)

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
    plt.savefig(os.path.join(SAT_PLOT_DIR, 'saturation_Ts.png'))
    plt.show()

    # ---------- P-h diagrams (saturation curves) ----------
    plt.figure(figsize=(8,6), dpi=120)

    for fluid, data in sat_data.items():
        valid_l = np.isfinite(data['h_l']) & np.isfinite(data['P']) & (data['P'] > 0)
        valid_v = np.isfinite(data['h_v']) & np.isfinite(data['P']) & (data['P'] > 0)

        if np.any(valid_l):
            plt.plot(data['h_l'][valid_l] * 1e-3, data['P'][valid_l] * 1e-5, '-', label=f'{fluid} liq')
        if np.any(valid_v):
            plt.plot(data['h_v'][valid_v] * 1e-3, data['P'][valid_v] * 1e-5, '--', label=f'{fluid} vap')

    plt.xlabel('Specific enthalpy (kJ/kg)')
    plt.ylabel('Pressure (bar)')
    plt.title('P-h saturation curves')
    plt.yscale('log')
    plt.grid(True, which='both', ls=':')
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys(), ncol=2, fontsize='small')
    plt.tight_layout()
    plt.savefig(os.path.join(SAT_PLOT_DIR, 'saturation_Ph.png'))
    plt.show()


    # ---------- h-s diagrams (saturation curves) ----------
    plt.figure(figsize=(8,6), dpi=120)

    for fluid, data in sat_data.items():
        valid_l = np.isfinite(data['h_l']) & np.isfinite(data['s_l'])
        valid_v = np.isfinite(data['h_v']) & np.isfinite(data['s_v'])

        if np.any(valid_l):
            plt.plot(data['h_l'][valid_l] * 1e-3, data['s_l'][valid_l] * 1e-3, '-', label=f'{fluid} liq')
        if np.any(valid_v):
            plt.plot(data['h_v'][valid_v] * 1e-3, data['s_v'][valid_v] * 1e-3, '--', label=f'{fluid} vap')

    plt.xlabel('Specific enthalpy (kJ/kg)')
    plt.ylabel('Specific entropy (kJ/kg·K)')
    plt.title('h-s saturation curves')
    plt.grid(True, ls=':')
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys(), ncol=2, fontsize='small')
    plt.tight_layout()
    plt.savefig(os.path.join(SAT_PLOT_DIR, 'saturation_hs.png'))
    plt.show()







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

    # ---------- mass flow, power ----------
    motor_efficiency = 0.90  # assume motor/electrical efficiency for electrical power estimate


    # Sweep and compute mass flow and powers
    Q_ref_W = Q_ref_kW * 1000.0
    massflow_results = {}
    for fluid in HFOs_avail + HFCs_avail:
        m_dot_list = []
        P_mech_list = []
        P_elec_list = []
        COP_list = []
        valid_Tev = []
        for Te in T_evap_range:
            st = cycle_states(fluid, Te, T_cond_C, eta_isentropic)
            if st is None:
                m_dot_list.append(np.nan)
                P_mech_list.append(np.nan)
                P_elec_list.append(np.nan)
                COP_list.append(np.nan)
                continue
            q_in = st['q_in_per_kg']  # J/kg
            wkg = st['w_comp_per_kg']  # J/kg
            # require positive refrigeration effect and compressor work
            if (not np.isfinite(q_in)) or q_in <= 1e-9:
                m_dot_list.append(np.nan)
                P_mech_list.append(np.nan)
                P_elec_list.append(np.nan)
                COP_list.append(np.nan)
                continue
            m_dot = Q_ref_W / q_in  # kg/s needed to meet Q_ref
            P_mech = m_dot * wkg  # W mechanical to compressor shaft
            P_elec = P_mech / motor_efficiency  # W electrical input estimate
            cop = q_in / wkg if wkg > 0 and np.isfinite(q_in/wkg) else np.nan
            m_dot_list.append(m_dot)
            P_mech_list.append(P_mech)
            P_elec_list.append(P_elec)
            COP_list.append(cop)
            valid_Tev.append(Te)
        massflow_results[fluid] = dict(
            Te=T_evap_range,
            m_dot=np.array(m_dot_list),
            P_mech=np.array(P_mech_list),
            P_elec=np.array(P_elec_list),
            COP=np.array(COP_list)
        )

    # Plot mass flow vs evaporator temperature
    plt.figure(figsize=(10,5), dpi=120)
    for fluid, d in massflow_results.items():
        valid = np.isfinite(d['m_dot']) & (d['m_dot'] > 0)
        if not np.any(valid):
            print(f"No valid mass flow points for {fluid}.")
            continue
        plt.plot(d['Te'][valid], d['m_dot'][valid], label=fluid)
    plt.xlabel('Evaporator Temperature (°C)')
    plt.ylabel('Mass flow rate (kg/s)')
    plt.title(f'Mass flow required for {Q_ref_kW} kW cooling vs Evaporator Temp')
    plt.grid(True, ls=':')
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(PLOT_DIR, 'massflow_vs_Tevap.png'))
    plt.show()

    # Plot mechanical and electrical power vs evaporator temperature
    plt.figure(figsize=(10,5), dpi=120)
    for fluid, d in massflow_results.items():
        valid = np.isfinite(d['P_mech']) & (d['P_mech'] > 0)
        if not np.any(valid):
            continue
        plt.plot(d['Te'][valid], d['P_mech'][valid]/1000.0, '--', label=f'{fluid} mech kW')
        plt.plot(d['Te'][valid], d['P_elec'][valid]/1000.0, '-', label=f'{fluid} elec kW')
    plt.xlabel('Evaporator Temperature (°C)')
    plt.ylabel('Power (kW)')
    plt.title('Compressor mechanical and electrical power vs Evap Temp')
    plt.grid(True, ls=':')
    plt.legend(ncol=2, fontsize='small')
    plt.tight_layout()
    plt.savefig(os.path.join(PLOT_DIR, 'power_vs_Tevap.png'))
    plt.show()


    # Prepare data frame to collect results

    rows = []
    for refrig in HFOs:
        for Te in T_evap_range:
            metrics = compute_cycle_metrics(refrig, Te, T_cond_C)

            # cooling capacity for fixed mass flow
            if np.isfinite(metrics['q_e_Jkg']):
                Qdot_kW = metrics['q_e_Jkg'] * m_dot_fixed / 1000.0

            else:
                Qdot_kW = np.nan

            rows.append({
                'refrigerant': refrig,
                'Te_C': Te,
                'COP': metrics['COP'],
                'q_e_Jkg': metrics['q_e_Jkg'],
                'Qdot_kW_fixed_mdot': Qdot_kW,
                'w_comp_Jkg': metrics['w_comp_Jkg']
            })

    df = pd.DataFrame(rows)


    # Save CSV summary

    csv_out = os.path.join(PLOT_DIR, 'hfo_comparative_raw.csv')
    df.to_csv(csv_out, index=False)


    # Plot 1: Cooling capacity vs Te (fixed mass flow)

    plt.figure(figsize=(7,4))
    for refrig, g in df.groupby('refrigerant'):
        valid = np.isfinite(g['Qdot_kW_fixed_mdot'])
        if not np.any(valid):
            continue
        plt.plot(g.loc[valid,'Te_C'], g.loc[valid,'Qdot_kW_fixed_mdot'], marker='o', label=refrig)

    plt.xlabel('Evaporator temp (°C)')
    plt.ylabel(f'Cooling capacity (kW) at ṁ={m_dot_fixed} kg/s')
    plt.title('Cooling capacity vs Evaporator Temperature (HFO comparison)')
    plt.grid(True)
    plt.legend(fontsize='small')
    plt.tight_layout()

    png1 = os.path.join(PLOT_DIR, 'hfo_cooling_vs_Te.png')

    plt.savefig(png1)
    plt.close()


    # Compute ranking by mean cooling capacity across Te range and mean COP

    summary = df.groupby('refrigerant').agg(
        mean_Q_kW = ('Qdot_kW_fixed_mdot', 'mean'),
        mean_COP = ('COP', 'mean')
    ).reset_index()

    # Rank primarily by mean cooling capacity, tie-break by mean COP

    summary['rank'] = summary[['mean_Q_kW','mean_COP']].rank(method='dense', ascending=False).mean(axis=1).rank(method='dense').astype(int)
    summary = summary.sort_values(['mean_Q_kW','mean_COP'], ascending=False).reset_index(drop=True)
    summary_csv = os.path.join(PLOT_DIR, 'hfo_comparative_summary.csv')
    summary.to_csv(summary_csv, index=False)



    # # Display the summary to the user in a small table
    # print(summary.head())


    # print("Saved raw data:", csv_out)
    # print("Saved plots:", png1)
    # print("Saved summary:", summary_csv)


    # ----------------- Extended ranking using multiple metrics -----------------
    # Requires: massflow_results dict computed earlier, and summary DataFrame created above

    # List fluids to score. Use fluids present in the summary table (keeps ordering consistent)
    fluids_to_score = summary['refrigerant'].tolist()

    # Prepare arrays to hold mean power (electrical) and mean mass flow
    mean_power_kW = []
    mean_mdot = []

    for refrig in fluids_to_score:
        # default nan
        mpow = np.nan
        mmdot = np.nan
        if refrig in massflow_results:
            arr = massflow_results[refrig]
            # P_elec stored in W
            P_elec = np.array(arr.get('P_elec', []), dtype=float)
            m_dot = np.array(arr.get('m_dot', []), dtype=float)
            # compute mean over finite entries
            if P_elec.size > 0:
                valid_pow = np.isfinite(P_elec)
                if np.any(valid_pow):
                    mpow = np.nanmean(P_elec[valid_pow]) / 1000.0  # convert to kW
            if m_dot.size > 0:
                valid_mdot = np.isfinite(m_dot)
                if np.any(valid_mdot):
                    mmdot = np.nanmean(m_dot[valid_mdot])
        mean_power_kW.append(mpow)
        mean_mdot.append(mmdot)

    # attach to summary DataFrame
    summary_ext = summary.copy()
    summary_ext['mean_power_kW'] = mean_power_kW
    summary_ext['mean_mdot_kg_s'] = mean_mdot

    # Normalise each metric to 0-1 where 1 is best
    # For mean_Q_kW and mean_COP higher is better
    # For mean_power_kW and mean_mdot_kg_s lower is better, so invert the normalisation


    summary_ext['score_Q'] = normalize(summary_ext['mean_Q_kW'], higher_is_better=True)
    summary_ext['score_COP'] = normalize(summary_ext['mean_COP'], higher_is_better=True)
    summary_ext['score_power'] = normalize(summary_ext['mean_power_kW'], higher_is_better=False)
    summary_ext['score_mdot'] = normalize(summary_ext['mean_mdot_kg_s'], higher_is_better=False)

    # Combined score: equal weights for all four metrics. Adjust weights if you want.
    weights = {'Q': 1.0, 'COP': 1.0, 'power': 1.0, 'mdot': 1.0}
    weight_sum = sum(weights.values())

    summary_ext['combined_score'] = (
        weights['Q'] * summary_ext['score_Q']
    + weights['COP'] * summary_ext['score_COP']
    + weights['power'] * summary_ext['score_power']
    + weights['mdot'] * summary_ext['score_mdot']
    ) / weight_sum

    # Final ranking (1 = best)
    summary_ext = summary_ext.sort_values('combined_score', ascending=False).reset_index(drop=True)
    summary_ext['final_rank'] = np.arange(1, len(summary_ext) + 1)

    # Reorder columns for readability
    cols_order = [
        'final_rank', 'refrigerant',
        'mean_Q_kW', 'mean_COP',
        'mean_power_kW', 'mean_mdot_kg_s',
        'score_Q', 'score_COP', 'score_power', 'score_mdot',
        'combined_score'
    ]
    # keep only available columns in case some are missing
    cols_order = [c for c in cols_order if c in summary_ext.columns]
    summary_ext = summary_ext[cols_order]

    # Save and print
    out_csv_ext = os.path.join(PLOT_DIR, 'hfo_comparative_summary_extended.csv')
    summary_ext.to_csv(out_csv_ext, index=False)

    print("Extended ranking (best first):")
    print(summary_ext.head(len(summary_ext)))
    print("Saved extended summary:", out_csv_ext)


    # # # ---------- Export CSVs ----------
    # for fluid, data in sat_data.items():
    #     df = pd.DataFrame({
    #         'T_C': T_C,
    #         'T_K': data['T_K'],
    #         'P_Pa': data['P'],
    #         'h_liq_Jkg': data['h_l'],
    #         'h_vap_Jkg': data['h_v'],
    #         's_liq_JkgK': data['s_l'],
    #         's_vap_JkgK': data['s_v'],
    #     })
    #     fname = os.path.join(OUTDIR, f'saturation_{fluid.replace("/","_")}.csv')
    #     df.to_csv(fname, index=False)
    #     print("Wrote:", fname)

    # print("Plots saved to:", PLOT_DIR)
    # print("Done.")    