import copy
import numpy as np
from typing import Tuple, TYPE_CHECKING

from ....utils.analysis import (
    get_unique_f_suffix,
    print_results,
    get_fprime,
    get_local_optimum,
    compute_SEC_obj,
    get_var_kwargs,
    transfer_temp_files_to_results,
    solve_single_case_with_timer,
    get_results_report,
)

if TYPE_CHECKING:
    from ...common import BaseHeatPump


#######################################################################################################
# Public API
#######################################################################################################


def analyse_RHP_case(
    HP_fn, case: dict, results_dir: str = "results"
) -> Tuple["BaseHeatPump", float]:
    """Optimise the gas-cooler offset for a single case and solve it."""
    dt_ls = []

    dT_gc_out_min_offset = (
        0.0 if "dT_gc_out_min_offset" not in case else case["dT_gc_out_min_offset"]
    )

    vars_input = {
        "dT_gc_out_min_offset": {"xi": dT_gc_out_min_offset, "bnds": (0.0, 40.0)}
    }
    case["additional_data"] = {"vars_input": vars_input}

    try:
        # Attempt to calculate the gradient at the boundary, if it fails,
        # assume that further optimisation could be possible
        fprime, dt = get_fprime(
            HP_fn,
            param=case,
            f_obj=compute_SEC_obj,
        )
        dt_ls.append(dt)
    except Exception:
        fprime = 1

    if fprime > 0.0:
        result, dt = get_local_optimum(
            HP_fn,
            param=case,
            vars_input=vars_input,
            tol=1e-4,
            f_obj=compute_SEC_obj,
            method="COBYQA",
        )
        dt_ls.append(dt)
        dT_gc_out_min_offset = result.x[0]
        is_solved = result.success
    elif fprime <= 0.0:
        dT_gc_out_min_offset = 0.0
        result = {
            "nfev": 2,
            "message": "Success",
            "success": True,
        }
        is_solved = True
    else:
        is_solved = False

    if is_solved:
        hp, dt = solve_single_case_with_timer(
            HP_fn, case | get_var_kwargs(dT_gc_out_min_offset, case)
        )
        dt_ls.append(dt)
        obj = compute_SEC_obj(hp, case)

        hp.dQdW = _get_dQdW(hp)
        get_results_report(
            hp=hp,
            param=case,
            result=result,
            obj=obj,
            dt=sum(dt_ls),
            print_to_terminal=False,
            results_dir=results_dir,
        )
        return hp, dt
    else:
        return None, dt


def analyse_unsolved_RHP_cases(
    HP_fn,
    debug: bool = True,
    results_dir: str = "results",
):
    """Retry previously unsolved cases and optionally emit debug output."""
    unsolved_cases = HP_fn.UNSOLVED_CASE_LIST_FN(results_dir=results_dir)

    if unsolved_cases is not None and unsolved_cases != []:
        for case in unsolved_cases:
            args = {
                "refrigerant": str(case["refrigerant"]),
                "dT_GC_min_tar": float(case["dT_GC_min_tar"]),
                "T_air_tar": float(case["T_air_tar"]),
                "T_air_inlet": float(case["T_air_inlet"]),
                "T_source_in": float(case["T_source_in"]),
                "n_comp": int(case["n_comp"]),
                "ihx_on": int(case["ihx_on"]),
                "dTmin_rhx": float(case["dTmin_rhx"]),
                "dT_rhx1_cold": float(case["dT_rhx1_cold"]),
                "dT_rhx2_cold": float(case["dT_rhx2_cold"]),
                "dT_rhx3_cold": float(case["dT_rhx3_cold"]),
                "comp2_power_frac": float(case["comp2_power_frac"]),
                "comp3_power_frac": float(case["comp3_power_frac"]),
                "comp_PR_uniform": bool(case["comp_PR_uniform"]),
                "comp_model_index": int(case["comp_model_index"]),
                "motor_eff": float(case["motor_eff"]),
                "sink": str(case["sink"]),
                "source": str(case["source"]),
                "Tref_evap_inlet": float(case["Tref_evap_inlet"]),
                "debug": bool(debug),
            }

            hp, _ = analyse_RHP_case(HP_fn, args, results_dir=results_dir)
            if isinstance(hp, HP_fn):
                if hp.converged and debug:
                    hp.get_report()
                    fsuffix = get_unique_f_suffix(args) + "-RESOLVED"
                    print_results(hp.output, fsuffix, results_dir=results_dir)

        transfer_temp_files_to_results("Cases-Resolved", results_dir=results_dir)


#######################################################################################################
# Local helper functions
#######################################################################################################


def _get_dQdW(hp0: "BaseHeatPump") -> float:
    """Finite-difference estimate of heating capacity to work sensitivity."""
    dtemp = 0.001
    args = copy.deepcopy(hp0.args)
    args["T_air_tar"] = hp0.T_air_tar + dtemp
    args["P_H_est"] = hp0.comp_discharge.pressure
    args["T_gc_out"] = hp0.gas_cooler1.source_out.temperature
    HP_fn = hp0.__class__
    hp1, _ = solve_single_case_with_timer(HP_fn, args)
    if hp1.feasible:
        tcop = [
            (p.T_air_tar - p.T_air_inlet) / p.COP if p.COP > 0 else np.inf
            for p in [hp0, hp1]
        ]
        dqdw = (
            dtemp / (tcop[1] - tcop[0])
            if tcop[1] - tcop[0] != 0 and tcop[0] != np.inf and tcop[1] != np.inf
            else 0.0
        )
    else:
        dqdw = 0.0
    return dqdw
