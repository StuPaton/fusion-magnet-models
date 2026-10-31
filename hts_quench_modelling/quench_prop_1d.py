######################################################
# ----------------- Import libraries -----------------
######################################################
import copy
import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolours
from scipy.integrate import solve_ivp, cumulative_trapezoid
from scipy.interpolate import interp1d
from scipy.sparse import diags
from dataclasses import dataclass, field, replace, fields
from IPython.display import display
from itertools import product
from pathlib import Path


######################################################
# ----------------- Model parameters -----------------
######################################################
def default_materials():
    return { # density = kg m^(-3)
        'cu': {'normal_cond': True, 'density': 8960, 'base_rrr': 150, 'depth': 40},
        'ag': {'normal_cond': True, 'density': 10490, 'base_rrr': 30, 'depth': 2},
        'ha_c276': {'normal_cond': True, 'density': 8890, 'base_rrr': 0, 'depth': 50},
        'ybco': {'normal_cond': False, 'density': 5450, 'base_rrr': 0, 'depth': 2}
    }


try:
    BASE_DIR = Path(__file__).resolve().parent
except NameError:          # interactive session: use the current working directory
    BASE_DIR = Path.cwd()

@dataclass
class CoreSolverParams:
    # ----------------- Generally static parameters -----------------
    node_count: int = 401 #TODO: Use tape_length and node_count to set node_gap indirectly?
    node_gap: float = 1 # mm
    node_gap_scalar: float = 10**(-3)
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

    materials: dict = field(default_factory=default_materials)

    avg_prop_velocity_tol: float = 0.02 # i.e. percentage tolerance versus final velocity to use to set window for averaging velocity

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
    post_heater_time: float = 1 # s

    properties_loc: str = BASE_DIR / "data" / "material_properties.csv"
    properties_magres_loc: str = BASE_DIR / "data" / "material_properties_magres.csv"


@dataclass
class OtherParams:
    # ----------------- Variable parameters -----------------
    tape_limit_temp: float = 200 # K
    detection_voltage_per_node: float = 0.01 # V
    detection_delay: float = 0.005 # s
    detection_validation: float = 0.01 # s, TODO: Will in future check voltage above limit for this time
    decay_time: float = 0.0036 # s (total_current * inductance / maximum dump V)
    tap_nodes: int = 1 # number of nodes covered by taps, from centre (i.e. set to 1 if just want to use centre node)


@dataclass
class SolverSolution:
    temps: np.ndarray
    times: np.ndarray
    temp_grid: np.ndarray
    norm_zone_lengths: np.ndarray
    voltages: np.ndarray
    U_temps: np.ndarray
    T_c: float
    total_current: float
    heater_pulse: float
    heater_duration: float
    end_normal_zone: float
    avg_prop_velocity: float
    end_prop_velocity: float
    time_to_300k: float
    end_peak_temp: float
    end_peak_dT_dt: float
    mqe: float
    charts: plt.Figure


@dataclass
class DetectionSolution:
    tap_voltages: np.ndarray
    detection_time: float
    detection_peak_temp: float
    detection_validn_peak_temp: float
    detection_validn_delay_peak_temp: float
    detection_validn_delay_decay_peak_temp: float
    target_detection_time: float
    time_to_tape_limit_temp: float


######################################################
# ----------------- Functions -----------------
######################################################
def print_results(result):  
    print("Final normal zone = " + str(round(result['end_normal_zone'] * 100, 1)) + " cm")
    print("Average propagation velocity = " + str(round(result['avg_prop_velocity'] * 100, 1)) + " cm/s")
    print("Final propagation velocity = " + str(round(result['end_prop_velocity'] * 100, 1)) + " cm/s")
    if np.isfinite(result['time_to_300k']):
        print("Time to 300K = " + str(int(result['time_to_300k']*1000)) + " ms")
    else:
        print("Time to 300K = N/A")
    print("Final peak temperature = " + str(round(result['end_peak_temp'], 1)) + " K")
    print("Final peak temperature derivative = " + str(round(result['end_peak_dT_dt'], 1)) + " K/s")
    if np.isfinite(result['mqe']):
        print("Minimum quench energy = " + str(int(result['mqe']*1000)) + " mJ (pulse = " + str(round(result['heater_pulse'], 2)) + " W, duration = " + str(int(result['heater_duration']*1000)) + " ms)")
    print("Detection time = " + str(int(result['detection_time']*1000)) + " ms")
    print("Target detection time = " + str(int(result['target_detection_time']*1000)) + " ms")
    print("Time to tape limit temperature = " + str(int(result['time_to_tape_limit_temp']*1000)) + " ms")
    print("Peak temperature at detection = " + str(int(result['detection_peak_temp'])) + " K")
    print("Peak temperature after dump resistor = " + str(int(result['detection_validn_delay_decay_peak_temp'])) + " K")

    display(result['charts'])


def interp1d_l_bounded_extrap_r_fill(x_new, x, y, left_lb, right_fill):
    result = interp1d(x, y, fill_value="extrapolate")(x_new)
    result[x_new > np.max(x)] = right_fill
    result[result < left_lb] = left_lb
    return result


def simulate_tape(p):
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
            [alpha_x[n] / B_in * (B_in/B_i_x[n])**p_x[n] * (1 - B_in/B_i_x[n])**q_x[n] *
            (B_i_x[n]/B_i0_x[n])**gamma_x[n] for n in range(len(p_x))
        ]).T
        J_c_x[B_in/B_i_x.T >= 1] = 0
        g_B_T = g[0] + g[1] * np.exp(-g[2] * B_in * np.exp(g[3] * t_in)) # Anisotropy factor
        J_c_B_T_theta = np.min(J_c_x, axis=-1) + np.maximum(0, -np.diff(J_c_x)[:,0]) / (1 + (abs(B_angle_in - math.pi/2) / g_B_T)**v)
        I_c = J_c_B_T_theta * area_crosssection_sc * 10**12 # Convert from MA/mm^2 to A, so need to back out M->1, and 10^3 in both width and depth

        if take_all:
            return I_c
        else:
            return I_c[0]


    def run_solver(pulse): # dT/dt solver
        def dT_dt(t, n_t): # Define the dT/dt ODE inside run_solver that it can use pulse
            n_t_b = (n_t[:-1] + n_t[1:]) / 2
            n_k = interp1d(temp_grid, therm_cond, fill_value="extrapolate")(n_t_b)
            n_cp = interp1d(temp_grid, spec_heat, fill_value="extrapolate")(n_t)
            n_rho_per_l = interp1d(temp_grid, resist_per_length, fill_value="extrapolate")(n_t)
            n_I_nc = interp1d(np.arange(int(p.YBCO_T_c0)), I_normal, bounds_error=False, fill_value=(np.nan, total_current))(n_t)
            n_P = (n_I_nc**2 * n_rho_per_l - p.cooling_coeff * (n_t - p.T_bath)) * p.node_gap * p.node_gap_scalar
            if t < p.heater_duration:
                n_P[p.node_count // 2] += pulse
            
            a = -np.hstack((0, n_k)) * np.hstack((0, np.diff(n_t)))
            b = np.hstack((n_k, 0)) * np.hstack((np.diff(n_t), 0))
            
            return ((a + b) / p.node_gap / p.node_gap_scalar + n_P) / n_cp


        def normal_temp_near_end(t, y):
            return y[int(p.node_count*0.05)] - T_norm # Used 5% of length from end to avoid uptick due to adiabatic boundary reflection
        normal_temp_near_end.terminal = True
        normal_temp_near_end.direction = 1

        jac_sparsity = diags([1,1,1], [-1,0,1], shape=(p.node_count,p.node_count)) 
        sol1 = solve_ivp(dT_dt, (0, p.heater_duration), np.full(p.node_count, p.T_0), method='Radau', events=normal_temp_near_end, max_step=p.heater_duration/10, t_eval=np.linspace(0, p.heater_duration, int(p.heater_duration*1000)), jac_sparsity=jac_sparsity)
        assert sol1.t[-1] == p.heater_duration, "sol1 end state does not cover full heater pulse duration"
        sol2 = solve_ivp(dT_dt, (p.heater_duration, max_time), sol1.y[:,-1], method='Radau', events=normal_temp_near_end, t_eval=np.linspace(p.heater_duration, max_time, int((max_time - p.heater_duration)*1000)), jac_sparsity=jac_sparsity)
        temps = np.hstack((sol1.y, sol2.y[:,1:]))
        times = np.hstack((sol1.t, sol2.t[1:]))
        peak_dT_dt = np.diff(np.max(temps, axis=0)) / np.diff(times)

        return temps, times, peak_dT_dt


    def final_peak_dT_dt_positive(pulse):
        return run_solver(pulse)[-1][-1] > 10**(-10)


    def normal_zone_grows(pulse):
        temps, times, _ = run_solver(pulse)
        nzl = np.sum(temps > T_norm, axis=0) * p.node_gap * p.node_gap_scalar
        return nzl[-1] - nzl[np.argmax(times == p.heater_duration)] > 10**(-10)


    # ----------------- Load materials info -----------------
    mats = copy.deepcopy(p.materials)
    mat_props = pd.read_csv(p.properties_loc, header=[0,1,2]).to_dict('list')
    mat_props = {((key[0].lower(), key[1].lower(), int(key[2])) if int(key[2]) != 0
        else (key[0].lower(), key[1].lower(), mats[key[1].lower()]['base_rrr'])):
        np.array(mat_props[key]) for key in mat_props if key[1].lower() in mats
    }
    mat_magres = pd.read_csv(p.properties_magres_loc).to_dict('list')
    mat_magres = {key.lower(): mat_magres[key] for key in mat_magres if key.lower() in mats or key.lower() == 'log_b_rrr'}

    max_temp = len(mat_props[list(mat_props)[0]])-1
    temp_grid = np.arange(max_temp+1)

    # ----------------- Setup other info -----------------
    area_crosssection_sc = mats['ybco']['depth'] * p.width * p.depth_scalar * p.width_scalar # m^2
    area_crosssection_all = sum(y['depth'] for y in mats.values()) * p.width * p.depth_scalar * p.width_scalar # m^2
    area_flat = p.node_gap * p.width * p.node_gap_scalar * p.width_scalar # m^2
    total_current = p.J * area_crosssection_all * 10**6 # A
    max_time = p.heater_duration + p.post_heater_time

    # ----------------- Setup tables to interpolate from at each stage -----------------
    # Specific heat, resistivity, and thermal conductivity for tape cross section, for 0-300K
    for x, y in mats.items():
        base_rrr = y['base_rrr']
        y['rho'] = mat_props[('rho', x,  base_rrr)]
        log_b_rrr = np.log10(p.B * y['rho'][p.base_rrr_temp] / y['rho']) if x in mat_magres else np.zeros(p.node_count)
        magres = interp1d(mat_magres['log_b_rrr'], mat_magres[x], fill_value="extrapolate")(log_b_rrr) if x in mat_magres else np.zeros_like(temp_grid)

        y['spec_heat'] = mat_props[('cp', x, base_rrr)]
        y['rho'] *= (1 + magres)
        y['therm_cond'] = mat_props[('k', x, base_rrr)] / (1 + magres)

    spec_heat = sum(y['spec_heat'] * y['density'] * y['depth'] for y in mats.values()) * p.width * p.node_gap * p.depth_scalar * p.width_scalar * p.node_gap_scalar
    therm_cond = sum(y['therm_cond'] * y['depth'] for y in mats.values()) * p.width * p.depth_scalar * p.width_scalar
    face_resist = sum(y['depth'] * (1/y['rho']) for y in mats.values() if y['normal_cond']) * p.width * p.depth_scalar * p.width_scalar
    resist_per_length = 1 / face_resist

    # Calculate improvement factor using quoted engineering current density, then determine critical current given input field strength and angle for all temps, and scale using factor
    I_e = p.J_e * p.width * p.width_scalar * 10**3
    I_c_raw = YBCO_crit_current(p.YBCO_T_c0, p.J_e_T, p.J_e_B, p.J_e_angle)
    f_improv = I_e / I_c_raw
    I_c_B_T_theta = YBCO_crit_current(p.YBCO_T_c0, np.arange(int(p.YBCO_T_c0)), p.B, p.B_angle) * f_improv

    # Normal conducting current
    I_sc_bounds = np.hstack((np.zeros((int(p.YBCO_T_c0),1)), np.full((int(p.YBCO_T_c0),1), total_current)))
    I_sc_bounds[I_c_B_T_theta == 0] = 0
    I_sc = np.average(I_sc_bounds, axis=1)
    solve_current = total_current - p.E_0 / resist_per_length[:int(p.YBCO_T_c0)] * (I_sc / I_c_B_T_theta)**p.n_value - I_sc
    solve_current[I_c_B_T_theta == 0] = 0

    while np.max(np.diff(I_sc_bounds)) > 10**(-7):
        I_sc_bounds[np.arange(int(p.YBCO_T_c0)),(solve_current < 0)*1] = I_sc
        I_sc = np.average(I_sc_bounds, axis=1)
        solve_current = total_current - p.E_0 / resist_per_length[:int(p.YBCO_T_c0)] * (I_sc / I_c_B_T_theta)**p.n_value - I_sc
        solve_current[I_c_B_T_theta == 0] = 0
    I_normal = total_current - I_sc

    # ----------------- Solve for temperatures -----------------
    T_cs = np.argmax(I_c_B_T_theta < total_current)
    T_c = np.argmax(I_c_B_T_theta == 0)
    T_norm = (T_cs + T_c)/2
    heater_pulse = p.heater_pulse

    if p.find_mqe:
        # Set upper / lower bounds and check limits
        spec_heat_integral = cumulative_trapezoid(spec_heat)
        fn_cp_int = interp1d(temp_grid[1:], spec_heat_integral, fill_value="extrapolate")
        therm_diff = (therm_cond / spec_heat)[int(p.T_0)+1:]
        therm_diff_scalar = (int(math.sqrt(max(therm_diff * p.heater_duration / p.node_gap / p.node_gap_scalar))) + 1) * 2 # i.e. sqrt(alpha t) / node_gap, where alpha needs node gap applied too, multiplied by 2 for margin
        lb = (fn_cp_int(T_cs) - fn_cp_int(p.T_0)) / p.heater_duration
        ub = (spec_heat_integral[-1] - fn_cp_int(p.T_0)) * therm_diff_scalar / p.heater_duration # spread out over nodes
        assert not final_peak_dT_dt_positive(0), "Tape quenches with no heater pulse"
        assert normal_zone_grows(ub), "No quench even for really high energy"
        
        # Find the window where flips between no quench -> quench
        low, high = (lb+ub)/4, (lb+ub)/2
        while final_peak_dT_dt_positive(low) and low > lb:
            high, low = low, low/2
        while not final_peak_dT_dt_positive(high) and high < ub:
            low, high = high, high*2
            
        # Converge on the point where the flip happens
        while high - low > 10**(-4):
            mid = (low + high)/2
            if final_peak_dT_dt_positive(mid):
                high = mid
            else: 
                low = mid
        heater_pulse = high
        
    temps, times, peak_dT_dt = run_solver(heater_pulse)

    # ----------------- Calculate values for detection (normal zone lengths, voltages and U curve) -----------------
    mid = p.node_count // 2
    norm_zone_lengths = np.array([mid - interp1d(t, np.arange(mid), fill_value="extrapolate")(T_norm) if max(t) > T_norm else 0 for t in temps[:mid].T]) * p.node_gap * p.node_gap_scalar * 2
    resistances = np.array([interp1d(temp_grid, resist_per_length, fill_value="extrapolate")(temps_step) for temps_step in temps]) * p.node_gap * p.node_gap_scalar
    voltages = np.array([interp1d(np.arange(int(p.YBCO_T_c0)), I_normal, bounds_error=False, fill_value=(np.nan, total_current))(temps_step) for temps_step in temps]) * resistances
    U_temps = cumulative_trapezoid(spec_heat / resist_per_length) / p.node_gap / p.node_gap_scalar
    U_temps -= np.interp(T_c, temp_grid[1:], U_temps) # i.e. set it relative to critical temperature

    # ----------------- Calculate summary values -----------------
    prop_vel = np.diff(norm_zone_lengths) / np.diff(times) / 2
    avg_vel_window = np.abs((prop_vel - prop_vel[-1]) / prop_vel[-1]) < p.avg_prop_velocity_tol
    idx = np.argmin(avg_vel_window[::-1]) # find the last point where there's a velocity which isn't in the window to use for averaging

    end_normal_zone = norm_zone_lengths[-1]
    avg_prop_velocity = np.average(prop_vel[-idx:])
    end_prop_velocity = prop_vel[-1]
    time_to_300k = np.diff(interp1d_l_bounded_extrap_r_fill([T_c, 300], np.max(temps, axis=0), times, 0, np.nan))[0] # i.e. time to reach 300K after reaching critical temperature
    end_peak_temp = max(temps[:,-1])
    end_peak_dT_dt = peak_dT_dt[-1]
    mqe = heater_pulse * p.heater_duration if p.find_mqe else np.nan

    # ----------------- Setup charts -----------------
    if p.make_charts:
        charts = plt.figure(figsize=(8,6))
        ax1 = charts.add_subplot(2,2,1)
        ax2 = charts.add_subplot(2,2,2,sharey=ax1)
        ax3 = charts.add_subplot(2,2,3)
        ax4 = charts.add_subplot(2,2,4,sharey=ax3)
        norm = mcolours.Normalize(vmin=0, vmax=times[-1]*10**3)
        cmap = plt.get_cmap('plasma_r')

        ax1.plot(times*1000, np.max(temps, axis=0))
        for i, t in enumerate(times):
            ax2.plot(temps[:,i], color=cmap(norm(t*10**3)))
        ax3.plot(times*1000, norm_zone_lengths / p.node_gap_scalar)
        ax4.plot(times[:-1]*1000, prop_vel*100)

        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        charts.colorbar(sm, ax=ax2, label='Time (ms)')
        ax1.set_ylabel('K')
        ax1.set_yscale('log')
        ax1.set_title('Peak temperature')
        ax2.tick_params(labelleft=False)
        ax2.set_title('Mid tape temperature profile')
        ax2.set_xlim([int(p.node_count)/4-1,int(p.node_count)*3/4-1])
        ax3.set_ylabel('mm')
        ax3.set_title('Normal zone length')
        ax4.tick_params(labelleft=False)
        ax4.yaxis.set_label_position("right")
        ax4.set_ylabel('cm/s')
        ax4.set_title('Propagation velocity')
        plt.subplots_adjust(wspace=0, hspace=0.3)
    else:
        charts = plt.figure()
    plt.close(charts)

    return SolverSolution(
        temps=temps,
        times=times,
        temp_grid=temp_grid,
        norm_zone_lengths=norm_zone_lengths,
        voltages=voltages,
        U_temps=U_temps,
        T_c=T_c,
        total_current=total_current,
        heater_pulse=heater_pulse,
        heater_duration=p.heater_duration,
        end_normal_zone=end_normal_zone,
        avg_prop_velocity=avg_prop_velocity,
        end_prop_velocity=end_prop_velocity,
        time_to_300k=time_to_300k,
        end_peak_temp=end_peak_temp,
        end_peak_dT_dt=end_peak_dT_dt,
        mqe=mqe,
        charts=charts
    )


def detection_stats(p, s):
    fn_U_at_temps = interp1d(s.temp_grid[1:], s.U_temps, bounds_error=False, fill_value=(0, np.nan))
    fn_temps_at_U = interp1d(s.U_temps, s.temp_grid[1:], fill_value="extrapolate")
    fn_time_to_peak_temps = interp1d(np.max(s.temps, axis=0), s.times, bounds_error=False, fill_value=(0, np.nan))
    fn_peak_temps_at_time = interp1d(s.times, np.max(s.temps, axis=0), fill_value="extrapolate")

    U_limit_temp = fn_U_at_temps(p.tape_limit_temp)[()]
    time_to_T_c = fn_time_to_peak_temps(s.T_c)[()]
    time_to_tape_limit_temp = fn_time_to_peak_temps(p.tape_limit_temp) - time_to_T_c
    target_detection_time = U_limit_temp / s.total_current**2 - p.detection_delay - p.detection_validation - p.decay_time/2

    mid = s.voltages.shape[0] // 2
    tap_voltages = np.sum(s.voltages[mid:mid+p.tap_nodes], axis=0)
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
        target_detection_time=target_detection_time,
        time_to_tape_limit_temp=time_to_tape_limit_temp
    )


def shallow_asdict(obj):
    return {f.name: getattr(obj, f.name) for f in fields(obj)}


def run_param_grid(base_core, base_other, change_core={}, change_other={}):
    results = []
    change_core = {x: [y] if np.ndim(y) == 0 else y for x, y, in change_core.items()}
    change_other = {x: [y] if np.ndim(y) == 0 else y for x, y, in change_other.items()}
    for grid_changes_core in product(*change_core.values()):
        curr_change_core = dict(zip(change_core, grid_changes_core))
        solution_solver = simulate_tape(replace(base_core, **curr_change_core))
        print("Core run complete.")
        for grid_changes_other in product(*change_other.values()):
            curr_change_other = dict(zip(change_other, grid_changes_other))
            solution_detection = detection_stats(replace(base_other, **curr_change_other), solution_solver)
            results.append({**curr_change_core, **curr_change_other, **shallow_asdict(solution_solver), **shallow_asdict(solution_detection)})
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