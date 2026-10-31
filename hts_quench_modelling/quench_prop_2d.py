######################################################
# ----------------- Import libraries -----------------
######################################################
import copy
import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolours
from matplotlib.collections import PolyCollection
from matplotlib.animation import FuncAnimation
from scipy.integrate import solve_ivp, cumulative_trapezoid
from scipy.interpolate import interp1d
import scipy.sparse as sparse
from dataclasses import dataclass, field, replace, fields
from IPython.display import display, HTML
from itertools import product
from typing import Optional
from joblib import Parallel, delayed
from pathlib import Path


######################################################
# ----------------- Model parameters -----------------
######################################################
def default_materials():
    return { # density = kg m^(-3)
        'cu': {'conductor': True, 'normal_conductor': True, 'anisotropic': False, 'density': 8960, 'base_rrr': 150, 'depth': 40},
        'ag': {'conductor': True, 'normal_conductor': True, 'anisotropic': False, 'density': 10490, 'base_rrr': 30, 'depth': 2},
        'ha_c276': {'conductor': True, 'normal_conductor': True, 'anisotropic': False, 'density': 8890, 'base_rrr': 0, 'depth': 50},
        'ybco': {'conductor': True, 'normal_conductor': False, 'anisotropic': True, 'density': 5450, 'base_rrr': 0, 'depth': 2},
        'kapton': {'conductor': False, 'normal_conductor': False, 'anisotropic': False, 'density': 1420, 'base_rrr': 0, 'depth': 25},
    }


try:
    BASE_DIR = Path(__file__).resolve().parent
except NameError:          # interactive session: use the current working directory
    BASE_DIR = Path.cwd()

@dataclass
class CoreSolverParams:
    # ----------------- Generally static parameters -----------------
    inner_rad: float = 20 * 10**(-3) # m
    turn_count: int = 20
    nodes_per_turn: int = 133 # Set at broadly 1mm spacing
    
    depth_scalar: float = 10**(-6)
    width: float = 5.5 # mm
    width_scalar: float = 10**(-3)

    n_value: float = 30
    E_0: float = 10**(-5) # V/m

    J_e: float = 56.4 # A/mm, quoted engineering current density for per mm tape width, with below inputs
    J_e_B: float = 20
    J_e_angle: float = 0 # i.e. perpendicular
    J_e_T: float = 4.5
    YBCO_T_c0: float = 93 # K

    base_rrr_temp: float = 273 # K
    temp_limit: int = 300 # K
    normal_distance_limit: float = 10 * 10**(-3) # m, distance from ends of the tape (on inner and outer turns) that will trigger end of solver

    materials: dict = field(default_factory=default_materials)

    avg_prop_velocity_tol: float = 0.02 # i.e. percentage tolerance versus final velocity to use to set window for averaging velocity

    radial_conductivity_on: bool = True # Set to false if want to switch radial conductivity off - mainly for regression testing versus 1D model
    radial_node_lengths: bool = True # Set to false if want to use fixed node lengths - mainly for regression testing versus 1D model

    solver_steps: float = 0.001 # s
    solver_rtol: float = 1e-4

    check_tests: bool = False

    # ----------------- Variable parameters -----------------
    make_charts: bool = False

    B: float = 0 #T
    B_angle: float = 0
    J: float = 0 # A/mm^2
    T_0: float = 4.2 # K
    T_bath: float = 4.2 # K
    cooling_coeff: float = 0 # W / (m K), no cooling used initially

    find_mqe: bool = False # If False, use the below heater_pulse, otherwise find the minimum value that leads to a quench (propagating normal zone)
    heater_pulse: float = 0 # W
    heater_duration: float = 0.01 # s
    heater_turn_vs_centre: int = 0
    heater_node_vs_centre: int = 0
    post_heater_time: float = 2 # s

    contact_res_per_area: float = 0 # m^2 K / W, contact resistance per unit area
    contact_count: int = 2 # i.e. number of contact points, dependent on construction type

    properties_loc: str = BASE_DIR / "data" / "material_properties.csv"
    properties_caxis_loc: str = BASE_DIR / "data" / "material_properties_caxis.csv"
    properties_magres_loc: str = BASE_DIR / "data" / "material_properties_magres.csv"


@dataclass
class OtherParams:
    # ----------------- Variable parameters -----------------
    detection_voltage_per_node: float = 0.01 # V
    detection_delay: float = 0.005 # s
    detection_validation: float = 0.01 # s, TODO: Will in future check voltage above limit for this time
    decay_time: float = 0.0036 # s (total_current * inductance / maximum dump V)
    tap_nodes: int = 1 # number of nodes covered by taps, from heater (i.e. set to 1 if just want to use the heater node)


@dataclass
class SolverSolution:
    outcome: str
    temps: np.ndarray
    times: np.ndarray
    heater_node: int
    end_time: float
    temp_grid: np.ndarray
    norm_zone_rad: np.ndarray
    norm_zone_long: np.ndarray
    norm_zone_total: np.ndarray
    end_norm_zone_rad: float
    end_norm_zone_long: float
    end_norm_zone_total: float
    avg_prop_velocity_rad: np.ndarray
    avg_prop_velocity_long: np.ndarray
    end_prop_velocity_rad: float
    end_prop_velocity_long: float
    end_growth_rate_total: float
    voltages: np.ndarray
    U_temps: np.ndarray
    T_c: float
    total_current: float
    heater_pulse: float
    heater_duration: float
    temp_limit: int
    time_to_temp_limit_from_T_c: float
    end_peak_temp: float
    end_peak_dT_dt: float
    mqe: float
    fixed_charts: Optional[plt.Figure] = None
    evolve_charts: Optional[FuncAnimation] = None


@dataclass
class DetectionSolution:
    tap_voltages: np.ndarray
    detection_time: float
    detection_peak_temp: float
    detection_validn_peak_temp: float
    detection_validn_delay_peak_temp: float
    detection_validn_delay_decay_peak_temp: float
    target_detection_time: float


class SolverFailure(Exception):
    pass


######################################################
# ----------------- Functions -----------------
######################################################
def print_results(result):  
    print("Final normal zone = " + str(round(result['end_norm_zone_total'] * 100, 1)) + " cm")
    print("Average propagation velocity (longitudinal) = " + str(round(result['avg_prop_velocity_long'] * 100, 1)) + " cm/s")
    print("Final propagation velocity (longitudinal) = " + str(round(result['end_prop_velocity_long'] * 100, 1)) + " cm/s")
    print("Average propagation velocity (radial) = " + str(round(result['avg_prop_velocity_rad'] * 100, 1)) + " cm/s")
    print("Final propagation velocity (radial) = " + str(round(result['end_prop_velocity_rad'] * 100, 1)) + " cm/s")
    if np.isfinite(result['time_to_temp_limit_from_T_c']):
        print("Time to " + str(result['temp_limit']) + "K = " + str(int(result['time_to_temp_limit_from_T_c']*1000)) + " ms")
    else:
        print("Time to " + str(result['temp_limit']) + "K = N/A")
    print("Final peak temperature = " + str(round(result['end_peak_temp'], 1)) + " K")
    print("Final peak temperature derivative = " + str(round(result['end_peak_dT_dt'], 1)) + " K/s")
    if np.isfinite(result['mqe']):
        print("Minimum quench energy = " + str(int(result['mqe']*1000)) + " mJ (pulse = " + str(round(result['heater_pulse'], 2)) + " W, duration = " + str(int(result['heater_duration']*1000)) + " ms)")
    print("Detection time = " + str(int(result['detection_time']*1000)) + " ms")
    print("Target detection time = " + str(int(result['target_detection_time']*1000)) + " ms")
    print("Peak temperature at detection = " + str(int(result['detection_peak_temp'])) + " K")
    print("Peak temperature after dump resistor = " + str(int(result['detection_validn_delay_decay_peak_temp'])) + " K")

    display(result['fixed_charts'])
    display(HTML(result['evolve_charts'].to_jshtml()))


def interp1d_l_bounded_extrap_r_fill(x_new, x, y, left_lb, right_fill):
    result = interp1d(x, y, fill_value="extrapolate")(x_new)
    result[x_new > np.max(x)] = right_fill
    result[result < left_lb] = left_lb
    return result


def simulate_tape(p):
    T_c0 = int(p.YBCO_T_c0)
    t_ctr = p.turn_count
    n_ctr = p.nodes_per_turn
    width = p.width * p.width_scalar
    node_count = t_ctr * n_ctr
    t_heater = t_ctr // 2 + p.heater_turn_vs_centre
    n_heater = n_ctr // 2 + p.heater_node_vs_centre
    node_h = t_heater * n_ctr + n_heater
    assert t_heater in range(t_ctr), "Heater is located outside the coil"
    assert n_heater in range(n_ctr), "Heater is located outside the coil"

    # ----------------- All functions nested inside simulate_tape -----------------
    def YBCO_crit_current(T_c0, t_in, B_in, B_angle_in): # Calculate critical current given input critical temp, temperature, field and field angle
        # Critical surface fit for Fujikura REBCO tape. Parameters from Table A.3 of:
        # J. van Nugteren, "High Temperature Superconductor Accelerator Magnets",
        # PhD thesis, University of Twente, 2016.
        alpha = 0.1
        n = [1, 1.4, 4.45] # n, n_1, n_2
        g = [0.03, 0.25, 0.06, 0.058] # g_0, g_1, g_2, g_3
        v = 1.85

        B_i0_x = [250, 140] # T, in ab and c planes
        alpha_x = [68.3, 1.86] #
        p_x = [1, 0.5]
        q_x = [5, 2.5]
        gamma_x = [1.63, 2.44]

        if type(t_in) == float:
            t_in = np.array([t_in])
            take_all = False
        else:
            take_all = True

        t_div_T_c0 = t_in / T_c0
        B_i_x = np.array([
            ((1 - t_div_T_c0**n[1])**n[2] + alpha * (1 - t_div_T_c0**n[0])) * B_i0_x[0],
            (1 - t_div_T_c0**n[0]) * B_i0_x[1]
        ])
        J_c_x = np.array(
            [alpha_x[i] / B_in * (B_in/B_i_x[i])**p_x[i] * (1 - B_in/B_i_x[i])**q_x[i] *
            (B_i_x[i]/B_i0_x[i])**gamma_x[i] for i in range(len(p_x))
        ]).T
        J_c_x[B_in/B_i_x.T >= 1] = 0
        g_B_T = g[0] + g[1] * np.exp(-g[2] * B_in * np.exp(g[3] * t_in)) # Anisotropy factor
        J_c_B_T_theta = np.min(J_c_x, axis=-1) + np.maximum(0, -np.diff(J_c_x)[:,0]) / (1 + (abs(B_angle_in - math.pi/2) / g_B_T)**v)
        I_c = J_c_B_T_theta * area_crosssection_sc * 10**12 # Convert from MA/mm^2 to A, so need to back out M->1, and 10^3 in both width and depth

        if take_all:
            return I_c
        else:
            return I_c[0]


    def run_solver(pulse, n_lengths): # dT/dt solver       
        def get_outcome(sol, events):
            if sol.status == 1:
                outcome = events[np.argmax([len(e)>0 for e in sol.t_events])].label
            elif sol.success:
                outcome = "Reached time limit"
            else:
                raise SolverFailure(f"Solver failed (pulse={pulse:.3g} W, time={sol.t[-1]*1000:.1f} ms, max temp={np.max(sol.y[:node_count,-1]):.0f} K, min temp={np.min(sol.y[:node_count,-1]):.0f} K): {sol.message}")

            return outcome
        
        def dT_dt(t, n_t_j): # Define the dT/dt ODE inside run_solver that it can use pulse
            n_t = n_t_j[:node_count]
            n_t_avg_long = (n_t[:-1] + n_t[1:]) / 2
            n_t_avg_rad = (n_t[:-n_ctr] + n_t[n_ctr:]) / 2
            n_k_ab = np.interp(n_t_avg_long, temp_grid, therm_cond_ab) ## np.interp used here - runs are exited when temperature limit reached, so shouldn't go past the end of the material properties tables
            h_c = np.interp(n_t_avg_rad, temp_grid, heat_trans_coeff_c)
            n_cp = np.interp(n_t, temp_grid, spec_heat)
            n_rho_per_l = np.interp(n_t, temp_grid, resist_per_length)
            n_I_nc = np.interp(n_t, np.arange(T_c0), I_normal, right=total_current)
            n_P = (n_I_nc**2 * n_rho_per_l - p.cooling_coeff * (n_t - p.T_bath)) * n_lengths
            if t < p.heater_duration:
                n_P[node_h] += pulse

            a = -np.hstack((0, n_k_ab)) * np.hstack((0, np.diff(n_t))) / (n_lengths - d_total * math.pi / n_ctr**2)
            b = np.hstack((n_k_ab, 0)) * np.hstack((np.diff(n_t), 0)) / (n_lengths + d_total * math.pi / n_ctr**2)
            c = -np.hstack((np.zeros(n_ctr), h_c)) * np.hstack((np.zeros(n_ctr), n_t[n_ctr:] - n_t[:-n_ctr])) * (n_lengths - d_total * math.pi / n_ctr) * width
            d = np.hstack((h_c, np.zeros(n_ctr))) * np.hstack((n_t[n_ctr:] - n_t[:-n_ctr], np.zeros(n_ctr))) * (n_lengths + d_total * math.pi / n_ctr) * width

            dTdt = (a + b + c + d + n_P) / n_cp / n_lengths
            if p.check_tests:
                return np.hstack((dTdt, sum(n_P)))
            return dTdt

        # Checks functions to terminate solver if needed
        def reached_temp_limit(t, y):
            return np.max(y[:node_count]) - p.temp_limit
        reached_temp_limit.terminal = True
        reached_temp_limit.direction = 1
        reached_temp_limit.label = "Temperature has reached limit"

        left_idx = np.argmax(np.cumsum(n_lengths) > p.normal_distance_limit)
        right_idx = -(np.argmax(np.cumsum(n_lengths[::-1]) > p.normal_distance_limit)+1)

        def normal_temp_near_end(t, y):
            return max(y[left_idx], y[right_idx]) - T_norm
        normal_temp_near_end.terminal = True
        normal_temp_near_end.direction = 1
        normal_temp_near_end.label = "Normal temperature near end of the tape"

        def below_T_cs(t, y):
            return T_cs - np.max(y[:node_count])
        below_T_cs.terminal = True
        below_T_cs.direction = 1
        below_T_cs.label = "No quench - temperatures below T_cs"

        events1 = (reached_temp_limit, normal_temp_near_end)
        events2 = (reached_temp_limit, normal_temp_near_end, below_T_cs)

        jac_sparsity = sparse.diags([1,1,1,1,1], [-n_ctr,-1,0,1,n_ctr], shape=(node_count, node_count))
        y0 = np.full(node_count, p.T_0)
        if p.check_tests: # When checking tests, have to include an extra element on the end for the ohmic heating calc - then when temperatures are extracted, have to exclude this
            jac_sparsity = sparse.vstack([jac_sparsity, sparse.csr_matrix((1, node_count))])
            jac_sparsity = sparse.hstack([jac_sparsity, sparse.csr_matrix((node_count+1, 1))])      
            y0 = np.hstack((y0, 0))      

        sol1 = solve_ivp(dT_dt, (0, p.heater_duration), y0, method='Radau', events=events1, max_step=p.heater_duration/10,
                         t_eval=np.linspace(0, p.heater_duration, max(2, round(p.heater_duration/p.solver_steps))), jac_sparsity=jac_sparsity, rtol=p.solver_rtol)
        if sol1.status == 0:
            assert sol1.t[-1] == p.heater_duration, "sol1 end state does not cover full heater pulse duration"
        power_gen = None
        if sol1.status == 0 and np.max(sol1.y[:node_count,-1]) > T_cs:
            sol2 = solve_ivp(dT_dt, (p.heater_duration, max_time), sol1.y[:,-1], method='Radau', events=events2,
                             t_eval=np.linspace(p.heater_duration, max_time, max(2, round((max_time - p.heater_duration)/p.solver_steps))), jac_sparsity=jac_sparsity, rtol=p.solver_rtol)
            sol = sol2
            events = events2
            if p.check_tests:
                power_gen = np.hstack((sol1.y[-1], sol2.y[-1,1:]))
            sol.y = np.hstack((sol1.y[:node_count], sol2.y[:node_count,1:]))
            sol.t = np.hstack((sol1.t, sol2.t[1:]))
            outcome = get_outcome(sol, events)
        else:
            sol = sol1
            events = events1
            if p.check_tests:
                power_gen = sol1.y[-1]
            sol.y = sol1.y[:node_count]
            if sol1.status == 0:
                outcome = below_T_cs.label
            else:
                outcome = get_outcome(sol, events)

        if sol.status == 1:
            event_idx = np.argmax([len(e)>0 for e in sol.t_events])
            if p.check_tests:
                power_gen = np.hstack((power_gen, sol.y_events[event_idx][0][-1]))
            sol.y = np.hstack((sol.y, sol.y_events[event_idx][0][:node_count,None])) # Add the final state where the run ended, if an event happened
            sol.t = np.hstack((sol.t, sol.t_events[event_idx]))
        end_peak_dT_dt = (np.diff(np.max(sol.y, axis=0)) / np.diff(sol.t))[-1]
        time_to_temp_limit = np.nan if len(sol.t_events[0]) == 0 else sol.t_events[0][0]

        return {'outcome': outcome,
                'temps': sol.y,
                'times': sol.t,
                'end_peak_dT_dt': end_peak_dT_dt,
                'time_to_temp_limit': time_to_temp_limit,
                'power_gen': power_gen
        }


    def final_peak_dT_dt_positive(pulse, n_lengths):
        r = run_solver(pulse, n_lengths)
        return (r['end_peak_dT_dt'] > 10**(-10)) and (r['times'][-1] > p.heater_duration)


    def normal_zone_grows(pulse, n_lengths):
        r = run_solver(pulse, n_lengths)
        nzl = np.sum((r['temps'] > T_norm) * n_lengths[:,None], axis=0)
        return nzl[-1] - nzl[np.argmax(r['times'] == p.heater_duration)] > 10**(-10)
    

    # ----------------- Load materials info -----------------
    mats = copy.deepcopy(p.materials)
    for x in mats:
        mats[x]['depth'] *= p.depth_scalar
    mat_props = pd.read_csv(p.properties_loc, header=[0,1,2]).to_dict('list')
    mat_props = {((key[0].lower(), key[1].lower(), int(key[2])) if int(key[2]) != 0
        else (key[0].lower(), key[1].lower(), mats[key[1].lower()]['base_rrr'])):
        np.array(mat_props[key]) for key in mat_props if key[1].lower() in mats
    }
    mat_props_caxis = pd.read_csv(p.properties_caxis_loc, header=[0,1,2]).to_dict('list')
    mat_props_caxis = {((key[0].lower(), key[1].lower(), int(key[2])) if int(key[2]) != 0
        else (key[0].lower(), key[1].lower(), mats[key[1].lower()]['base_rrr'])):
        np.array(mat_props_caxis[key]) for key in mat_props_caxis if key[1].lower() in mats
    }
    mat_magres = pd.read_csv(p.properties_magres_loc).to_dict('list')
    mat_magres = {key.lower(): mat_magres[key] for key in mat_magres if key.lower() in mats or key.lower() == 'log_b_rrr'}

    max_temp = len(mat_props[list(mat_props)[0]])-1
    temp_grid = np.arange(max_temp+1)
    assert p.temp_limit <= temp_grid[-1], "Limit temperature is above the input data range"

    # ----------------- Setup other info -----------------
    max_time = p.heater_duration + p.post_heater_time
    d_total = sum(y['depth'] for y in mats.values())
    area_crosssection_sc = mats['ybco']['depth'] * width # m^2
    area_crosssection_current = sum(y['depth'] for y in mats.values() if y['conductor']) * width # m^2
    total_current = p.J * area_crosssection_current * 10**6 # A
    node_lengths = 2 * math.pi / n_ctr * (p.inner_rad + d_total * np.arange(node_count) * p.radial_node_lengths / n_ctr)
    node_gaps = (np.hstack((node_lengths, 0)) + np.hstack((0, node_lengths))) / 2 # Gaps between node centres, rather than node lengths

    # ----------------- Setup tables to interpolate from at each stage -----------------
    # Specific heat, resistivity, and thermal conductivity for tape cross section, for 0-300K
    for x, y in mats.items():
        base_rrr = y['base_rrr']
        if y['normal_conductor']:
            y['rho'] = mat_props[('rho', x,  base_rrr)]
        log_b_rrr = np.log10(p.B * y['rho'][p.base_rrr_temp] / y['rho']) if x in mat_magres else np.zeros(node_count)
        magres = interp1d(mat_magres['log_b_rrr'], mat_magres[x], fill_value="extrapolate")(log_b_rrr) if x in mat_magres else np.zeros_like(temp_grid)

        y['spec_heat'] = mat_props[('cp', x, base_rrr)]
        y['therm_cond_ab'] = mat_props[('k', x, base_rrr)] / (1 + magres)
        if y['normal_conductor']:
            y['rho'] *= (1 + magres)
        if p.radial_conductivity_on:
            if y['anisotropic']:
                y['therm_res_c'] = mat_props_caxis[('d_x_rho_th_slope', x, base_rrr)] + mat_props_caxis[('d_x_rho_th_intercept', x, base_rrr)] / y['depth']
            else:
                y['therm_res_c'] = 1 / y['therm_cond_ab']

    spec_heat = sum(y['spec_heat'] * y['density'] * y['depth'] for y in mats.values()) * width
    therm_cond_ab = sum(y['therm_cond_ab'] * y['depth'] for y in mats.values()) * width
    heat_trans_coeff_c = 1 / (sum(y['therm_res_c'] * y['depth'] for y in mats.values()) + p.contact_res_per_area * p.contact_count) if p.radial_conductivity_on else np.zeros_like(therm_cond_ab)
    face_resist = sum(y['depth'] * (1/y['rho']) for y in mats.values() if y['normal_conductor']) * width
    resist_per_length = 1 / face_resist

    # Calculate improvement factor using quoted engineering current density, then determine critical current given input field strength and angle for all temps, and scale using factor
    I_e = p.J_e * width * 10**3
    I_c_raw = YBCO_crit_current(p.YBCO_T_c0, p.J_e_T, p.J_e_B, p.J_e_angle)
    f_improv = I_e / I_c_raw
    I_c_B_T_theta = YBCO_crit_current(p.YBCO_T_c0, np.arange(T_c0), p.B, p.B_angle) * f_improv

    # Normal conducting current
    I_sc_bounds = np.hstack((np.zeros((T_c0,1)), np.full((T_c0,1), total_current)))
    I_sc_bounds[I_c_B_T_theta == 0] = 0
    I_sc = np.average(I_sc_bounds, axis=1)
    solve_current = total_current - p.E_0 / resist_per_length[:T_c0] * (I_sc / I_c_B_T_theta)**p.n_value - I_sc
    solve_current[I_c_B_T_theta == 0] = 0

    while np.max(np.diff(I_sc_bounds)) > 10**(-7):
        I_sc_bounds[np.arange(T_c0),(solve_current < 0)*1] = I_sc
        I_sc = np.average(I_sc_bounds, axis=1)
        solve_current = total_current - p.E_0 / resist_per_length[:T_c0] * (I_sc / I_c_B_T_theta)**p.n_value - I_sc
        solve_current[I_c_B_T_theta == 0] = 0
    I_normal = total_current - I_sc

    # ----------------- Solve for temperatures -----------------
    T_cs = np.argmax(I_c_B_T_theta < total_current)
    T_c = np.argmax(I_c_B_T_theta == 0)
    T_norm = (T_cs + T_c)/2
    heater_pulse = p.heater_pulse

    if p.find_mqe:
        # Set upper / lower bounds and check limits
        spec_heat_integral = cumulative_trapezoid(spec_heat) * node_lengths[node_h]
        therm_diff = (therm_cond_ab / spec_heat)[int(p.T_0)+1:]
        therm_diff_scalar = (int(math.sqrt(max(therm_diff * p.heater_duration)) / min(node_lengths)) + 1) # i.e. sqrt(alpha t) / node_gap, multiplied by 2 for margin, minimum node length used for conservatism
        lb = (np.interp(T_cs, temp_grid[1:], spec_heat_integral) - np.interp(p.T_0, temp_grid[1:], spec_heat_integral)) / p.heater_duration
        ub = (np.interp(p.temp_limit, temp_grid[1:], spec_heat_integral) - np.interp(p.T_0, temp_grid[1:], spec_heat_integral)) * therm_diff_scalar * 2 / p.heater_duration # Energy needed to heat one node from T_0 to temperature limit, scaled up to allow for diffusion over multiple nodes

        assert not final_peak_dT_dt_positive(0, node_lengths), "Tape quenches with no heater pulse"
        assert normal_zone_grows(ub, node_lengths), "No quench even for really high energy"
        
        # Find the window where flips between no quench -> quench
        print("================")
        try_low, try_high = lb * therm_diff_scalar, lb * therm_diff_scalar * 2
        high_quenches = False
        while try_low > lb and final_peak_dT_dt_positive(try_low, node_lengths):
            try_high, try_low = try_low, try_low/2
            high_quenches = True
            print('try_low quenches', round(try_low,4), round(try_high,4), round(try_high - try_low,4))
        if not high_quenches:
            while try_high < ub and not final_peak_dT_dt_positive(try_high, node_lengths):
                try_low, try_high = try_high, try_high*2
                print('try_high no quench', round(try_low,4), round(try_high,4), round(try_high - try_low,4))
            
        # Converge on the point where the flip happens
        print("================")
        while abs((try_high - try_low) / try_low) > 0.01:
            try_mid = (try_low + try_high)/2
            if final_peak_dT_dt_positive(try_mid, node_lengths):
                try_high = try_mid
                print('try_mid quenches', round(try_low,4), round(try_high,4), round(try_high - try_low,4))
            else: 
                try_low = try_mid
                print('try_mid no quench', round(try_low,4), round(try_high,4), round(try_high - try_low,4))
        heater_pulse = try_high
        
    r = run_solver(heater_pulse, node_lengths)
    outcome = r['outcome']
    temps = r['temps']
    times = r['times']
    end_peak_dT_dt = r['end_peak_dT_dt']
    time_to_temp_limit = r['time_to_temp_limit']


    # ----------------- Calculate values for detection (voltages and U curve) -----------------
    resistances = np.array([np.interp(np.minimum(temps_step, temp_grid[-1]), temp_grid, resist_per_length) for temps_step in temps]) * node_lengths[:,None] # Have limited temps_step value because it can overshoot by a small fraction in the solver
    voltages = np.array([np.interp(temps_step, np.arange(T_c0), I_normal, right=total_current) for temps_step in temps]) * resistances
    U_temps = cumulative_trapezoid(spec_heat / resist_per_length)
    U_temps -= np.interp(T_c, temp_grid[1:], U_temps) # i.e. set it relative to critical temperature

    # ----------------- Calculate summary values -----------------
    def interp_nzl(T_norm, t, node_gaps, node_h, centre_zone_only=False):
        if t[node_h] > T_norm:
            t_l = np.hstack((t[0], t))
            t_r = np.hstack((t, t[-1]))
            t_max = np.maximum(t_l, t_r)
            node_fracs = np.divide(t_max - T_norm, np.abs(t_r - t_l), out=(t_r > T_norm).astype(float), where=(t_r != t_l))
            node_fracs = np.clip(node_fracs, 0, 1)
            start = None
            end = None
            if centre_zone_only:
                if len(np.nonzero(node_fracs[:node_h+1] < 1)[0]) > 0:
                    start = np.nonzero(node_fracs[:node_h+1] < 1)[0][-1] # Last point transitioning to normal before heater
                if len(np.nonzero(node_fracs[node_h+1:] < 1)[0]) > 0:
                    end = np.nonzero(node_fracs[node_h+1:] < 1)[0][0] + node_h + 1 + 1 # First point transitioning from normal after heater, + additional 1 to index past the point
            return sum(node_gaps[start:end] * node_fracs[start:end])
        else:
            return 0

    norm_zone_rad = []
    norm_zone_long = []
    norm_zone_total = []
    for t in temps.T:
        norm_zone_long += [interp_nzl(T_norm, t, node_gaps, node_h, centre_zone_only=True)]
        norm_zone_total += [interp_nzl(T_norm, t, node_gaps, node_h)]
        idxs_rad = np.arange(t_ctr) * n_ctr + (node_h % n_ctr)
        node_gaps_rad = np.hstack((d_total/2, np.full(len(idxs_rad)-1, d_total), d_total/2))
        norm_zone_rad += [interp_nzl(T_norm, t[idxs_rad], node_gaps_rad, np.argwhere(idxs_rad == node_h)[0][0])]

    norm_zone_rad = np.array(norm_zone_rad)
    norm_zone_long = np.array(norm_zone_long)
    norm_zone_total = np.array(norm_zone_total)

    prop_vel_rad = np.diff(norm_zone_rad) / np.diff(times) / 2 #TODO: Develop this to only calculate before longitudinal makes full turn? (i.e. before the 2 interact). Review charts to check the progression looks sensible
    prop_vel_long = np.diff(norm_zone_long) / np.diff(times) / 2 # TODO: Develop this to only calculate before makes full turn? (i.e. before the 2 interact). Review charts to check the progression looks sensible
    growth_rate_total = np.diff(norm_zone_total) / np.diff(times)

    avg_vel_window_rad = np.abs((prop_vel_rad - prop_vel_rad[-1]) / prop_vel_rad[-1]) < p.avg_prop_velocity_tol
    avg_vel_window_long = np.abs((prop_vel_long - prop_vel_long[-1]) / prop_vel_long[-1]) < p.avg_prop_velocity_tol

    idx_rad = np.argmin(avg_vel_window_rad[::-1]) # Find the last point where there's a velocity which isn't in the window to use for averaging
    idx_long = np.argmin(avg_vel_window_long[::-1])

    avg_prop_velocity_rad = np.nan if (prop_vel_rad[-1] == 0 or not np.isfinite(prop_vel_rad[-1])) else np.average(prop_vel_rad[-idx_rad:])
    avg_prop_velocity_long = np.nan if (prop_vel_long[-1] == 0 or not np.isfinite(prop_vel_long[-1])) else np.average(prop_vel_long[-idx_long:])

    time_to_temp_limit_from_T_c = time_to_temp_limit - np.interp(T_c, np.max(temps, axis=0), times) # i.e. time to reach temperature limit after reaching critical temperature
    mqe = heater_pulse * p.heater_duration if p.find_mqe else np.nan

    solution = {
        'outcome': outcome,
        'temps': temps,
        'times': times,
        'heater_node': node_h,
        'end_time': times[-1],
        'temp_grid': temp_grid,
        'norm_zone_rad': norm_zone_rad,
        'norm_zone_long': norm_zone_long,
        'norm_zone_total': norm_zone_total,
        'end_norm_zone_rad': norm_zone_rad[-1],
        'end_norm_zone_long': norm_zone_long[-1],
        'end_norm_zone_total': norm_zone_total[-1],
        'avg_prop_velocity_rad': avg_prop_velocity_rad,
        'avg_prop_velocity_long': avg_prop_velocity_long,
        'end_prop_velocity_rad': prop_vel_rad[-1],
        'end_prop_velocity_long': prop_vel_long[-1],
        'end_growth_rate_total': growth_rate_total[-1],
        'voltages': voltages,
        'U_temps': U_temps,
        'T_c': T_c,
        'total_current': total_current,
        'heater_pulse': heater_pulse,
        'heater_duration': p.heater_duration,
        'temp_limit': p.temp_limit,
        'time_to_temp_limit_from_T_c': time_to_temp_limit_from_T_c,
        'end_peak_temp': max(temps[:,-1]),
        'end_peak_dT_dt': end_peak_dT_dt,
        'mqe': mqe
    }

    # ----------------- Setup charts -----------------
    if p.make_charts:
        # ----------------- Static charts -----------------
        #TODO: Add charts for radial / total normal zone and propagation velocities
        fixed_charts = plt.figure(figsize=(8,6))
        ax1 = fixed_charts.add_subplot(2,2,1)
        ax2 = fixed_charts.add_subplot(2,2,2,sharey=ax1)
        ax3 = fixed_charts.add_subplot(2,2,3)
        ax4 = fixed_charts.add_subplot(2,2,4,sharey=ax3)
        norm = mcolours.Normalize(vmin=0, vmax=times[-1]*10**3)
        cmap = plt.get_cmap('plasma_r')

        ax1.plot(times*1000, np.max(temps, axis=0))
        for i, t in enumerate(times):
            ax2.plot(temps[:,i], color=cmap(norm(t*10**3)))
        ax3.plot(times*1000, norm_zone_long*1000)
        ax4.plot(times[:-1]*1000, prop_vel_long*100)

        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        fixed_charts.colorbar(sm, ax=ax2, label='Time (ms)')
        ax1.set_ylabel('K')
        ax1.set_yscale('log')
        ax1.set_title('Peak temperature')
        ax2.tick_params(labelleft=False)
        ax2.set_title('Mid tape temperature profile')
        ax2.set_xlim([int(node_count)/4-1,int(node_count)*3/4-1])
        ax3.set_ylabel('mm')
        ax3.set_title('Normal zone length (longitudinal)')
        ax4.tick_params(labelleft=False)
        ax4.yaxis.set_label_position("right")
        ax4.set_ylabel('cm/s')
        ax4.set_title('Propagation velocity (longitudinal)')
        plt.subplots_adjust(wspace=0, hspace=0.3)
        plt.close(fixed_charts)

        solution['fixed_charts'] = fixed_charts

        # ----------------- Coil visualisation -----------------
        rad_in = p.inner_rad + d_total * (np.arange(node_count + n_ctr + 1) / n_ctr - 1/2) * p.radial_node_lengths
        angle_r = 2 * math.pi / n_ctr * (np.arange(node_count + n_ctr + 1) - 1/2) * p.radial_node_lengths

        coords = np.array([rad_in * np.cos(angle_r), rad_in * np.sin(angle_r)]).T
        nodes = [[coords[x] for x in [k, k+1, k+1+n_ctr, k+n_ctr]] for k in range(node_count)]

        fig = plt.figure(figsize=(12,8), constrained_layout=True)
        gs = fig.add_gridspec(1, 2, width_ratios=[4,1])
        ax_full = fig.add_subplot(gs[0,0])
        ax_zoom = fig.add_subplot(gs[0,1])
        norm = mcolours.Normalize(vmin=np.min(temps), vmax=np.max(temps))

        coils = []
        for ax, title in [(ax_full, "Full coil"), (ax_zoom, "Heater pulse region")]:
            coil = PolyCollection(
                nodes,
                array=temps[:,0],
                cmap=plt.get_cmap('plasma_r'),
                norm=norm,
                edgecolors="none",
                linewidths=0,
                antialiaseds=False
            )

            ax.add_collection(coil)
            ax.set_aspect("equal")
            ax.set_title(title)
            coils.append(coil)
        
        ax_full.autoscale_view()

        zoom_x = np.mean(np.array(nodes[node_h])[:,0])
        zoom_y = np.mean(np.array(nodes[node_h])[:,1])
        zoom_rad_x = t_ctr * d_total * 1.1 # TODO: Need to adapt this depending on where we're focusing in the coil - e.g. if it's on the top, need to flip these. And set radius to be fixed - i.e. middle of coil, but angle to be based on heater
        zoom_rad_y = (p.inner_rad + t_ctr * d_total) / 2

        ax_zoom.set_xlim(zoom_x - zoom_rad_x - p.heater_turn_vs_centre * d_total, zoom_x + zoom_rad_x + p.heater_turn_vs_centre * d_total)
        ax_zoom.set_ylim(zoom_y - zoom_rad_y, zoom_y + zoom_rad_y)
        
        fig.colorbar(coils[0], ax=[ax_full, ax_zoom], label="Temperature (K)")
        fig.suptitle("2D temperature evolution", fontsize=16, fontweight="bold")

        def update(frame):
            for coil in coils:
                coil.set_array(temps[:,frame])
            fig.suptitle(f"2D temperature evolution, {round(times[frame]*1000)} ms", fontsize=16, fontweight="bold")
            return coil, fig._suptitle

        evolve_charts = FuncAnimation(
            fig,
            update,
            frames=temps.shape[1],
            interval=100,
            blit=False
        )

        plt.close(fig)
        solution['evolve_charts'] = evolve_charts

    # ----------------- Tests to check -----------------
    if p.check_tests:
        # Diffusion test - length total and single turn
        lengths = []
        lengths_turn = []
        for time in (np.arange(20)+1)/1000:
            n = np.argmax(times>time)
            t = temps[:,n]
            threshold = min(t) + (max(t) - min(t))/2
            lengths.append(round(sum(node_lengths[t > threshold]), 4))
            lengths_turn.append(round(sum(node_lengths[node_h - n_ctr//2: node_h + n_ctr//2] * (t > threshold)[node_h - n_ctr//2: node_h + n_ctr//2]) ,4))

        print(np.array(lengths)*1000)
        print(np.array(lengths_turn)*1000)

        # Energy balance - fixed properties and no ohmic heating
        time_idx = np.argmax(times>p.heater_duration)
        heater_energy = p.heater_pulse * p.heater_duration
        heat1 = np.sum((temps - p.T_0) * spec_heat[0] * node_lengths[:,None], axis=0)[1:]
        heat_diff_pct1 = (heat1[time_idx:] - heater_energy) / heater_energy

        # Energy balance - temperature-dependent properties and no ohmic heating
        heat2 = np.array([sum((np.interp(np.minimum(t, temp_grid[-1]), temp_grid[1:], cumulative_trapezoid(spec_heat)) - 
                               np.interp(p.T_0, temp_grid[1:], cumulative_trapezoid(spec_heat))
                               ) * node_lengths) for t in temps.T])
        heat_diff_pct2 = (heat2[time_idx:] - heater_energy) / heater_energy

        # Energy balance - full properties and ohmic heating (excludes cooling)
        power_gen = r['power_gen']
        heat3 = heat2 - power_gen + heater_energy
        heat_diff_pct3 = (heat3[time_idx:] - heater_energy) / power_gen[time_idx:]

        print(np.round(heater_energy,2))
        for heat_diff_pct in (heat_diff_pct1, heat_diff_pct2, heat_diff_pct3):
            print(np.round(heat_diff_pct,4), round(min(heat_diff_pct),4), round(max(heat_diff_pct),4))
        plt.plot(heat_diff_pct3)

    return SolverSolution(**solution)


def detection_stats(p, s):
    #TODO: Replace these interp1d with np.interp if don't need the extrapolation function?
    fn_U_at_temps = interp1d(s.temp_grid[1:], s.U_temps, bounds_error=False, fill_value=(0, np.nan))
    fn_time_to_peak_temps = interp1d(np.max(s.temps, axis=0), s.times, bounds_error=False, fill_value=(0, np.nan))
    fn_peak_temps_at_time = interp1d(s.times, np.max(s.temps, axis=0), fill_value="extrapolate")

    U_limit_temp = fn_U_at_temps(s.temp_limit)[()]
    time_to_T_c = fn_time_to_peak_temps(s.T_c)[()]
    target_detection_time = U_limit_temp / s.total_current**2 - p.detection_delay - p.detection_validation - p.decay_time/2

    tap_voltages = np.sum(s.voltages[s.heater_node : s.heater_node + p.tap_nodes], axis=0) # TODO: Will need to update once insulation removed / partial, and consider both sides
    detection_time = interp1d_l_bounded_extrap_r_fill(p.detection_voltage_per_node * p.tap_nodes, tap_voltages, s.times, 0, np.nan) - time_to_T_c # Time from passing T_c to reaching voltage limit, which scales with tap spacing
    detection_peak_temp = fn_peak_temps_at_time(detection_time + time_to_T_c)[()]
    detection_validn_peak_temp = fn_peak_temps_at_time(detection_time + p.detection_validation + time_to_T_c)[()]
    detection_validn_delay_peak_temp = fn_peak_temps_at_time(detection_time + p.detection_validation + p.detection_delay + time_to_T_c)[()]
    detection_validn_delay_decay_peak_temp = fn_peak_temps_at_time(detection_time + p.detection_validation + p.detection_delay + p.decay_time/2 + time_to_T_c)[()]

    return DetectionSolution(
        tap_voltages=tap_voltages,
        detection_time=detection_time,
        detection_peak_temp=detection_peak_temp,
        detection_validn_peak_temp=detection_validn_peak_temp,
        detection_validn_delay_peak_temp=detection_validn_delay_peak_temp,
        detection_validn_delay_decay_peak_temp=detection_validn_delay_decay_peak_temp,
        target_detection_time=target_detection_time
    )


def shallow_asdict(obj):
    return {f.name: getattr(obj, f.name) for f in fields(obj)}


def run_param_grid(base_core, base_other, change_core={}, change_other={}, n_jobs=-1):
    change_core = {x: [y] if np.ndim(y) == 0 else y for x, y, in change_core.items()}
    change_other = {x: [y] if np.ndim(y) == 0 else y for x, y, in change_other.items()}
    changes_core = [dict(zip(change_core, c)) for c in product(*change_core.values())]
    changes_other = [dict(zip(change_other, c)) for c in product(*change_other.values())]

    solutions_core = Parallel(n_jobs=n_jobs)(delayed(simulate_tape)(replace(base_core, **cc)) for cc in changes_core)
    print("Core runs complete.")
    results = []
    for c_core, s_core in zip(changes_core, solutions_core):
        for c_other in changes_other:
            s_det = detection_stats(replace(base_other, **c_other), s_core)
            results.append({**c_core, **c_other, **shallow_asdict(s_core), **shallow_asdict(s_det)})

    return pd.DataFrame(results)


def tornado_plot(val, sens_results, sens_pct):
    tornado = []
    for x in sens_results:
        diffs = np.diff(sens_results[x][val])
        tornado.append({
            "parameter": x,
            "lb": -abs(diffs[0]),
            "change": abs(sum(diffs))
        })

    tornado_df = pd.DataFrame(tornado).sort_values("change").reset_index(drop=True)

    fig, ax = plt.subplots()
    y = np.arange(len(tornado_df))
    ax.barh(y, tornado_df["change"], left=tornado_df["lb"])
    ax.axvline(0, color="black", linewidth=1)
    ax.set_yticks(y)
    ax.set_yticklabels(tornado_df["parameter"])
    ax.set_xlabel(f'{val} sensitivity to {round(sens_pct*100)}% change in inputs')
    plt.tight_layout()
    plt.close(fig)

    return fig