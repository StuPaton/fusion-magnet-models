# ----------------- Import libraries -----------------
import math
import numpy as np
from dataclasses import replace
from pathlib import Path
from IPython.display import display
from quench_prop_1d import run_param_grid, tornado_plot, default_materials, print_results, CoreSolverParams, OtherParams

# ----------------- Setup base parameters for all runs -----------------
results = {}
defaults_core = CoreSolverParams()
defaults_other = OtherParams()

base_core_fig_2_37 = replace(defaults_core,
    B = 20,
    B_angle = 84/180 * math.pi,
    J = 970,
    T_0 = 4.2,
    heater_pulse = 5
)
base_core_fig_2_39 = replace(defaults_core,
    B = 15,
    B_angle = 84/180 * math.pi,
    J = 800,
    T_0 = 4.2,
    find_mqe = True,
    heater_duration = 0.05,
    J_e = 74,
    J_e_B = 18,
    J_e_T = 4.2
)
base_core_sect5 = replace(defaults_core,
    B = 17,
    B_angle = 86/180 * math.pi,
    T_0 = 4.5,
    heater_duration = 0.05,
    find_mqe = True
)
base_other = replace(defaults_other,
    tape_limit_temp = 200,
    detection_voltage_per_node = 0.01,
    detection_delay = 0.01,
    detection_validation = 0.005,
    decay_time = 0.0036,
    tap_nodes = 1
)

# ----------------- Setup parameters for and carry out individual runs -----------------
results = {
    "Main result, fig 2.37": run_param_grid(base_core_fig_2_37, base_other, change_core={"find_mqe": [False,True], "post_heater_time": 0.04, "make_charts": True}),
    "MQEs, fig 2.39": run_param_grid(base_core_fig_2_39, base_other, change_core={"J": np.arange(700,950,50)}),
    "Node gap x0.5": run_param_grid(base_core_fig_2_37, base_other, change_core={"node_gap": 0.5, "node_count": 801, "find_mqe": [False,True]}),
    "Node gap x2": run_param_grid(base_core_fig_2_37, base_other, change_core={"node_gap": 2, "node_count": 201, "find_mqe": [False,True]}),
    "Current and decay time sweep, fig 2.44": run_param_grid(base_core_fig_2_39, base_other, change_core={"J": np.arange(600,1050,50)}, change_other={"decay_time": np.arange(0.01,0.055,0.005), "detection_voltage_per_node": [0.01, 0.1]}),
    "Delay time sweep, fig 2.44": run_param_grid(base_core_fig_2_39, base_other, change_other={"detection_delay": np.arange(0.01,0.055,0.005), "detection_voltage_per_node": [0.01, 0.1]}),
    "Tap spacing sweep": run_param_grid(base_core_fig_2_39, base_other, change_other={"tap_nodes": np.arange(50)+1, "detection_voltage_per_node": [0.01, 0.1]}),
    "Current sweep, table 5.1": run_param_grid(base_core_sect5, base_other, change_core={"J": [97, 290, 484, 774, 1161, 1934]})
}

print_results(results['Main result, fig 2.37'].iloc[0])

sens_params = {'B': 0, 'B_angle': math.pi/2, 'J': 0, 'T_0': 0, 'heater_pulse': 0, 'J_e': 0} # Parameters to tweak, with baseline numbers to set sensitivities
sens_pct = 0.05
sens_results = {x: run_param_grid(base_core_fig_2_37, base_other, change_core={x: np.array([1-sens_pct, 1, 1+sens_pct]) * (getattr(base_core_fig_2_37, x) - sens_params[x]) + sens_params[x]}) for x in sens_params}
display(tornado_plot("avg_prop_velocity", sens_results, sens_pct))

print("Finished.")