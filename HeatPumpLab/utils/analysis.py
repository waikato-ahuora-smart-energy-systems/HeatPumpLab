"""Analysis utilities for running and post-processing heat pump cases."""

import math
import numpy as np
import os

from typing import TYPE_CHECKING
from scipy.optimize import minimize as loc_min, approx_fprime

from .properties_package import PropertyPackage
from .parallel_solving import solve_in_parallel
from .readWriteOutput import (
    transfer_temp_files_to_results,
    print_result_to_csv,
    print_result_to_terminal,
)
from .decorator_functions import timer_func
from .solve_profile import profile_section

if TYPE_CHECKING:
    from ..cycles.common.heatPump import BaseHeatPump

__all__ = ["get_parameter_sweep"]

_FAILED_OBJECTIVE_PENALTY = 1e6
_LOCAL_OPTIMISATION_MAX_EVALUATIONS = 80
_LOCAL_OPTIMISATION_MAX_ITERATIONS = 160


#######################################################################################################
# Public API
#######################################################################################################


def get_parameter_sweep(
    HP_fn,
    filename: str = "HP_results",
    param_kwargs: dict = None,
    n_cpu: int = os.cpu_count() - 1,
    results_dir: str = "results",
) -> None:
    """Solve a batch of cases concurrently and persist intermediate results."""
    file_name = filename
    if param_kwargs is None:
        param_kwargs = {}
    transfer_temp_files_to_results(file_name, results_dir=results_dir)
    cases = HP_fn.CASE_LIST_FN(
        file_name=file_name,
        results_dir=results_dir,
        **param_kwargs,
    )
    solve_in_parallel(
        HP_fn=HP_fn,
        cases=cases,
        f_solve=HP_fn.ANAYLSE_CASE_FN,
        n_cpu=n_cpu,
        solve_kwargs={"results_dir": results_dir},
    )
    transfer_temp_files_to_results(file_name, results_dir=results_dir)


@timer_func
def solve_single_case_with_timer(HP_fn, args: dict) -> "BaseHeatPump":
    """Instantiate and solve a single HeatPump case while timing execution."""
    hp = HP_fn(**args)
    hp.solve()
    return hp


def get_results_report(
    hp: "BaseHeatPump",
    param: dict,
    result: dict,
    obj: float,
    dt: float,
    print_to_terminal: bool = False,
    results_dir: str = "results",
) -> None:
    """Format, persist, and optionally print a solved case summary."""
    with profile_section("reporting"):
        hp.get_report()
        if "weight" in param["additional_data"]:
            hp.output.update(
                {
                    "Weight": param["additional_data"]["weight"],
                    "SEC_0": param["additional_data"]["SEC_0"],
                    "P_0": param["additional_data"]["P_0"],
                }
            )
        hp.output.update(
            {
                "Obj_function": obj,
                "Function evaluations": result["nfev"],
                "Message": result["message"],
                "Success": result["success"],
                "Time_to_optimise": dt,
            }
        )
        fsuffix = get_unique_f_suffix(param)
        print_results(
            results=dict(hp.output),
            fsuffix=fsuffix,
            print_to_terminal=print_to_terminal,
            results_dir=results_dir,
        )


#######################################################################################################
# Local helper functions
#######################################################################################################


def _get_initial_val(vars_input: dict) -> list:
    """Return the initial optimisation variables in argument order."""
    return [vars_input[key]["xi"] for key in vars_input]


def _get_bounds(vars_input: dict) -> list:
    """Return bounds corresponding to the optimisation variables."""
    return [vars_input[key]["bnds"] for key in vars_input]


def _get_obj(HP_fn, vars: np.ndarray, param: dict, f_obj) -> float:
    """Evaluate the objective function for a candidate variable vector."""
    with profile_section("objective_eval"):
        hp, _ = solve_single_case_with_timer(
            HP_fn, param | get_var_kwargs(vars.tolist(), param)
        )
        obj = float(f_obj(hp, param))
        if not math.isfinite(obj):
            obj = _FAILED_OBJECTIVE_PENALTY
        return obj * 1000


#######################################################################################################
# External helper functions
#######################################################################################################


@timer_func
def get_local_optimum(
    HP_fn, param: dict, vars_input: dict, tol: float, f_obj, method: str
) -> dict:
    """Run a local optimisation over the specified variable set."""
    return loc_min(
        fun=lambda x: _get_obj(HP_fn, x, param, f_obj),
        x0=_get_initial_val(vars_input),
        bounds=_get_bounds(vars_input),
        method=method,
        tol=tol,
        options={
            "maxfev": _LOCAL_OPTIMISATION_MAX_EVALUATIONS,
            "maxiter": _LOCAL_OPTIMISATION_MAX_ITERATIONS,
        },
    )


@timer_func
def get_fprime(HP_fn, param: dict, f_obj) -> dict:
    """Estimate the derivative of the objective at the lower bound."""
    vars_input = param.get("additional_data", {}).get("vars_input", {})
    x0 = _get_initial_val(vars_input) if vars_input else 0.0
    return (
        approx_fprime(xk=x0, f=lambda x: _get_obj(HP_fn, x, param, f_obj), epsilon=0.1)
        * -1
    )


def compute_SEC_obj(hp: "BaseHeatPump", param: dict) -> float:
    """Objective function: specific electricity consumption."""
    return hp.SEC


def get_var_kwargs(vars: list[float], param: dict) -> dict[str, float]:
    """Map optimisation vector back to keyword arguments."""
    d_var = dict()
    if not isinstance(vars, list | np.ndarray):
        vars = [vars]
    for i, key in enumerate(param["additional_data"]["vars_input"]):
        d_var[key] = vars[i]
    return d_var


def print_results(
    results: dict,
    fsuffix: str,
    print_to_terminal: bool = True,
    results_dir: str = "results",
) -> None:
    """Persist results to CSV and optionally print to stdout."""
    if print_to_terminal:
        print_result_to_terminal(results)
    print_result_to_csv(
        output=results,
        suffix=fsuffix,
        file_name="temp/R",
        results_dir=results_dir,
    )


def get_unique_f_suffix(param, fld_only=False) -> str:
    """Build a unique filename suffix encoding refrigerant and key inputs."""
    suffix: str = "-"

    r = PropertyPackage(param["refrigerant"])
    fld_comp = r.fld_name.replace(".FLD", "").split("*")
    if len(fld_comp) == 1:
        suffix = "-" + param["refrigerant"]
    else:
        for i in range(len(fld_comp)):
            suffix += str(fld_comp[i]) + str(r.fld_composition[i]).replace(".", "_")
            if i + 1 < len(fld_comp):
                suffix += "*"
    if fld_only:
        for ele in param:
            if ele != "additional_data" and ele != "refrigerant":
                suffix = suffix + "-" + "{:.0f}".format(param[ele])
    return suffix
