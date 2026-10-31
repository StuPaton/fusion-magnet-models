# ----------------- Import libraries -----------------
import math
from dataclasses import replace
from pathlib import Path
from IPython.display import display
from quench_prop_2d import run_param_grid, tornado_plot, default_materials, print_results, CoreSolverParams, OtherParams

# ----------------- Setup base parameters for all runs -----------------
try:
    BASE_DIR = Path(__file__).resolve().parent
except NameError:          # interactive session: use the current working directory
    BASE_DIR = Path.cwd()

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
base_core_diff_test = replace(defaults_core,
    properties_loc = BASE_DIR / "data" / "material_properties_fixedat20k.csv",
    properties_caxis_loc = BASE_DIR / "data" / "material_properties_caxis_fixedat20k.csv",
    B = 20,
    B_angle = 84/180 * math.pi,
    J = 0,
    T_0 = 20,
    heater_pulse = 250,
    heater_duration = 0.0001,
    inner_rad = 0.05,
    nodes_per_turn = 322,
    solver_steps = 0.0001,
    post_heater_time = 0.04,
    check_tests = True
)
mats_no_kapton = default_materials()
mats_no_kapton["kapton"]["depth"] = 0
base_core_regression_1d = replace(defaults_core,
    B = 20,
    B_angle = 84/180 * math.pi,
    J = 970,
    T_0 = 4.2,
    heater_pulse = 5,
    radial_conductivity_on = False,
    radial_node_lengths = False,
    inner_rad = 0.001 * 40 / 2 / math.pi,
    turn_count = 10,
    nodes_per_turn = 40,
    materials = mats_no_kapton
)
base_other = replace(defaults_other,
    detection_voltage_per_node = 0.01,
    detection_delay = 0.01,
    detection_validation = 0.005,
    decay_time = 0.0036,
    tap_nodes = 1
)

# ----------------- Setup parameters for and carry out individual runs -----------------
results = {
    "Main result": run_param_grid(base_core_fig_2_37, base_other, change_core={"find_mqe":[False,True], "make_charts": True}),
    "1D regression": run_param_grid(base_core_regression_1d, base_other, change_core={"find_mqe":[False,True]}),
    "Diffusion test": run_param_grid(base_core_diff_test, base_other),
    "Energy balance test": run_param_grid(base_core_fig_2_37, base_other, change_core={"check_tests": True}, n_jobs=1),
    "Node spacing test": run_param_grid(base_core_fig_2_37, base_other, change_core={"heater_pulse": 6, "inner_rad": 30.69e-3, "nodes_per_turn": [100, 200, 400], "find_mqe": [False,True]})
}

print_results(results['Main result'].iloc[0])

print("Finished.")