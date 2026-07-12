import copy
import math
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

_FIXED_PRESSURE_OBJECTIVE_FLAG = "fixed_discharge_pressure_objective"
_PRESSURE_SEED_FIELD = "P_H_est"
_EVAP_OUTLET_SEED_FIELD = "Tref_evap_outlet"


#######################################################################################################
# Public API
#######################################################################################################


def analyse_RHP_case(
    HP_fn,
    case: dict,
    results_dir: str = "results",
    use_fast_objective: bool = True,
    record_result: bool = True,
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
    objective_case = case
    using_fast_objective = False
    seed_hp = None
    seed_dt = None

    if use_fast_objective:
        objective_case, dt, using_fast_objective, seed_hp = _get_fast_objective_case(
            HP_fn, case, dT_gc_out_min_offset
        )
        if dt is not None:
            seed_dt = dt
            dt_ls.append(dt)

    try:
        # Attempt to calculate the gradient at the boundary, if it fails,
        # assume that further optimisation could be possible
        fprime, dt = get_fprime(
            HP_fn,
            param=objective_case,
            f_obj=compute_SEC_obj,
        )
        dt_ls.append(dt)
    except Exception:
        if using_fast_objective:
            objective_case = case
            using_fast_objective = False
            try:
                fprime, dt = get_fprime(
                    HP_fn,
                    param=objective_case,
                    f_obj=compute_SEC_obj,
                )
                dt_ls.append(dt)
            except Exception:
                fprime = 1
        else:
            fprime = 1

    if fprime > 0.0:
        result, dt = _get_local_optimum_with_exact_fallback(
            HP_fn, objective_case, case, vars_input, using_fast_objective
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
        if _can_reuse_exact_seed(seed_hp, dT_gc_out_min_offset):
            hp = seed_hp
            dt = seed_dt
        else:
            hp, dt = solve_single_case_with_timer(
                HP_fn, _get_exact_solve_args(case, dT_gc_out_min_offset, seed_hp)
            )
            dt_ls.append(dt)
        obj = compute_SEC_obj(hp, case)

        hp.dQdW = _get_dQdW(hp)
        if record_result:
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


def _get_fast_objective_case(
    HP_fn, case: dict, dT_gc_out_min_offset: float
) -> tuple[dict, float | None, bool, "BaseHeatPump | None"]:
    """Return a fixed-pressure objective case seeded by an exact solve."""
    try:
        hp, dt = solve_single_case_with_timer(
            HP_fn, _get_exact_solve_args(case, dT_gc_out_min_offset)
        )
    except Exception:
        return case, None, False, None

    seed = _get_exact_state_seed_from_hp(hp)
    if seed is None:
        return case, dt, False, hp

    return (
        case
        | seed
        | {_FIXED_PRESSURE_OBJECTIVE_FLAG: True},
        dt,
        True,
        hp,
    )


def _get_local_optimum_with_exact_fallback(
    HP_fn,
    objective_case: dict,
    exact_case: dict,
    vars_input: dict,
    using_fast_objective: bool,
):
    """Optimise using the objective case, falling back to exact if needed."""
    try:
        result, dt = get_local_optimum(
            HP_fn,
            param=objective_case,
            vars_input=vars_input,
            tol=1e-4,
            f_obj=compute_SEC_obj,
            method="COBYQA",
        )
    except Exception:
        if not using_fast_objective:
            raise
        return get_local_optimum(
            HP_fn,
            param=exact_case,
            vars_input=vars_input,
            tol=1e-4,
            f_obj=compute_SEC_obj,
            method="COBYQA",
        )

    if using_fast_objective and not result.success:
        fallback_result, fallback_dt = get_local_optimum(
            HP_fn,
            param=exact_case,
            vars_input=vars_input,
            tol=1e-4,
            f_obj=compute_SEC_obj,
            method="COBYQA",
        )
        return fallback_result, dt + fallback_dt

    return result, dt


def _get_exact_solve_args(
    case: dict,
    dT_gc_out_min_offset: float,
    seed_hp: "BaseHeatPump | None" = None,
) -> dict:
    """Return final-solve arguments with any objective-only flags removed."""
    args = case | get_var_kwargs(dT_gc_out_min_offset, case)
    args.pop(_FIXED_PRESSURE_OBJECTIVE_FLAG, None)
    seed = (
        _get_exact_state_seed_from_hp(seed_hp)
        if _can_reuse_exact_seed(seed_hp, dT_gc_out_min_offset)
        else None
    )
    return _apply_exact_state_seed(args, seed)


def _get_exact_state_seed_from_hp(hp: "BaseHeatPump | None") -> dict | None:
    """Return exact pressure and evaporator-outlet seeds from a solved model."""
    if hp is None or not _is_successfully_solved(hp):
        return None

    pressure = getattr(getattr(hp, "comp_discharge", None), "pressure", None)
    evap_outlet = getattr(hp, _EVAP_OUTLET_SEED_FIELD, None)
    if not _is_finite_positive(pressure) or not _is_finite_number(evap_outlet):
        return None

    return {
        _PRESSURE_SEED_FIELD: float(pressure),
        _EVAP_OUTLET_SEED_FIELD: float(evap_outlet),
    }


def _apply_exact_state_seed(args: dict, seed: dict | None) -> dict:
    """Return args updated with an exact-state seed, without mutating inputs."""
    return dict(args) if seed is None else dict(args) | seed


def _is_finite_number(value: object) -> bool:
    """Return True for values that can be represented as a finite float."""
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def _is_finite_positive(value: object) -> bool:
    """Return True for finite positive values."""
    return _is_finite_number(value) and float(value) > 0.0


def _can_reuse_exact_seed(
    seed_hp: "BaseHeatPump | None", dT_gc_out_min_offset: float
) -> bool:
    """Return True when the exact seed already solves the selected offset."""
    if seed_hp is None or not _is_successfully_solved(seed_hp):
        return False

    try:
        return (
            abs(
                float(seed_hp.dT_gc_out_min_offset)
                - float(dT_gc_out_min_offset)
            )
            <= 1e-9
        )
    except (AttributeError, TypeError, ValueError):
        return False


def _is_successfully_solved(hp: "BaseHeatPump") -> bool:
    """Return True when a seed solve produced usable exact results."""
    return bool(getattr(hp, "feasible", False) and getattr(hp, "converged", False))


def _get_dQdW(hp0: "BaseHeatPump") -> float:
    """Finite-difference estimate of heating capacity to work sensitivity."""
    dtemp = 0.001
    args = copy.deepcopy(hp0.args)
    args["T_air_tar"] = hp0.T_air_tar + dtemp
    args["P_H_est"] = hp0.comp_discharge.pressure
    args["T_gc_out"] = hp0.gas_cooler1.source_out.temperature
    args["Tref_evap_outlet"] = hp0.Tref_evap_outlet
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
