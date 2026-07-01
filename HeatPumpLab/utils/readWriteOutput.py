"""Helpers for printing and persisting heat pump case results."""

import os
from glob import glob
import csv
from pathlib import Path
from typing import Any

__all__ = [
    "print_result_to_terminal",
    "print_result_to_csv",
    "transfer_temp_files_to_results",
    "read_csv_to_results",
]

#######################################################################################################
# Public API
#######################################################################################################


def print_result_to_terminal(output: dict[str, Any]):
    """Pretty-print a result dictionary to the terminal."""
    key: str
    value: Any

    if isinstance(output, dict):
        for key, value in output.items():
            print(f"{key} {value}")
    else:
        print("Failed to print due to wrong output format.")


def print_result_to_csv(
    output: dict[str, Any] | None = None,
    file_name: str = "R",
    file_ext: str = ".csv",
    suffix: str = "",
    results_dir: str | Path = "results",
):
    """Append a result dictionary to a CSV file, creating it if needed."""
    if not isinstance(output, dict) or not output:
        print("Failed to append result to CSV file due to wrong output format.")
        return

    try:
        file_name = _get_results_path(file_name, file_ext, suffix, results_dir)
        file_name.parent.mkdir(parents=True, exist_ok=True)
        file_exists = os.path.isfile(file_name)
        if file_exists:
            with open(file_name, "a", newline="") as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=output.keys())
                writer.writerow(output)

        else:
            with open(file_name, "w", newline="") as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=output.keys())
                writer.writeheader()
                writer.writerow(output)

    except (OSError, csv.Error) as exc:
        print(f"Failed to append result to CSV file: {exc}")


def transfer_temp_files_to_results(
    file_name: str = "untitled", results_dir: str | Path = "results"
) -> None:
    """Move temporary CSV files into the main results file."""
    results_path = Path(results_dir)
    output_path = _get_results_path(file_name, ".csv", "", results_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    csv_file_list = glob(str(results_path / "temp" / "*.csv"))
    header_saved = False
    with open(output_path, "a") as fout:
        for f in reversed(csv_file_list):
            with open(f) as fin:
                header = next(fin)
                if not header_saved:
                    fout.write(header)
                    header_saved = True
                for line in fin:
                    fout.write(line)
            os.remove(f)


def read_csv_to_results(
    file_name: str = "untitled",
    file_ext: str = ".csv",
    results_dir: str | Path = "results",
) -> list[dict[str, Any]]:
    """Reads a CSV file from the results directory into a list of dictionaries."""
    file_path = _get_results_path(file_name, file_ext, "", results_dir)

    if not os.path.exists(file_path):
        print(f"File not found: {file_path}")
        return []

    try:
        with open(file_path, mode="r", newline="", encoding="utf-8") as csvfile:
            reader = csv.DictReader(csvfile)
            return list(reader)
    except (OSError, csv.Error) as e:
        print(f"Failed to read CSV file '{file_path}': {e}")
        return []


def _get_results_path(
    file_name: str,
    file_ext: str,
    suffix: str,
    results_dir: str | Path,
) -> Path:
    """Return a result path while preserving legacy results-relative file names."""
    results_path = Path(results_dir)
    relative_path = Path(file_name + suffix + file_ext)
    if relative_path.parts and relative_path.parts[0] == "results":
        relative_path = Path(*relative_path.parts[1:])
    return results_path / relative_path
