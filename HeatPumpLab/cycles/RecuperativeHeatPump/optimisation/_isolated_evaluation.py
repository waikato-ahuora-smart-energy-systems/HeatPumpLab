"""Process isolation for exact thermodynamic objective evaluations."""

from __future__ import annotations

from multiprocessing import get_context
from multiprocessing.connection import Connection
from typing import Any


class ExactEvaluationTimeout(TimeoutError):
    """Raised when an isolated exact solve exceeds its wall-clock limit."""


class ExactEvaluationWorkerError(RuntimeError):
    """Raised when an isolated exact-solve worker exits unexpectedly."""


def _exact_evaluation_worker(connection: Connection) -> None:
    """Evaluate requests until the parent asks this disposable worker to stop."""
    try:
        while True:
            request = connection.recv()
            if request is None:
                return
            try:
                # Imported here so the subprocess, rather than its parent, owns all
                # thermodynamic solve state used by exact objective evaluations.
                from ._objective import _OuterDesignObjective

                objective = _OuterDesignObjective(
                    max_evals=1,
                    base_case=request["base_case"],
                    refrigerant=request["refrigerant"],
                    bounds=request["bounds"],
                    results_dir=request["results_dir"],
                    use_fast_objective=request["use_fast_objective"],
                    penalty=request["penalty"],
                    cache_decimals=request["cache_decimals"],
                    molar_masses=request["molar_masses"],
                    exact_evaluation_timeout_seconds=None,
                )
                connection.send(("ok", objective.evaluate(request["x"])))
            except BaseException as exc:  # Keep parent alive after native-worker errors.
                connection.send(("error", f"{type(exc).__name__}: {exc}"))
    except (EOFError, BrokenPipeError, OSError):
        return
    finally:
        connection.close()


class IsolatedExactEvaluator:
    """Run exact solves in a recyclable subprocess with a hard timeout."""

    def __init__(
        self,
        *,
        timeout_seconds: float,
        max_tasks: int,
        context_name: str = "spawn",
    ) -> None:
        self.timeout_seconds = float(timeout_seconds)
        self.max_tasks = int(max_tasks)
        self.context_name = context_name
        self._process = None
        self._connection: Connection | None = None
        self._completed_tasks = 0

    def evaluate(self, request: dict[str, Any]):
        """Return one exact evaluation or terminate a non-returning worker."""
        self._ensure_started()
        assert self._connection is not None
        try:
            self._connection.send(request)
        except (BrokenPipeError, EOFError, OSError) as exc:
            self._stop(force=True)
            raise ExactEvaluationWorkerError(
                f"exact-evaluation worker could not accept request: {exc}"
            ) from exc

        if not self._connection.poll(self.timeout_seconds):
            self._stop(force=True)
            raise ExactEvaluationTimeout(
                "exact thermodynamic evaluation timed out after "
                f"{self.timeout_seconds:g} seconds"
            )

        try:
            status, payload = self._connection.recv()
        except (BrokenPipeError, EOFError, OSError) as exc:
            self._stop(force=True)
            raise ExactEvaluationWorkerError(
                "exact-evaluation worker exited without returning a result"
            ) from exc

        self._completed_tasks += 1
        if self._completed_tasks >= self.max_tasks:
            self._stop(force=False)
        if status != "ok":
            self._stop(force=True)
            raise ExactEvaluationWorkerError(str(payload))
        return payload

    def close(self) -> None:
        """Stop the current worker, if one exists."""
        self._stop(force=False)

    def _ensure_started(self) -> None:
        if self._process is not None and self._process.is_alive():
            return
        self._stop(force=True)
        context = get_context(self.context_name)
        parent_connection, child_connection = context.Pipe(duplex=True)
        process = context.Process(
            target=_exact_evaluation_worker,
            args=(child_connection,),
            daemon=True,
        )
        process.start()
        child_connection.close()
        self._process = process
        self._connection = parent_connection
        self._completed_tasks = 0

    def _stop(self, *, force: bool) -> None:
        process = self._process
        connection = self._connection
        self._process = None
        self._connection = None
        self._completed_tasks = 0
        if process is None:
            if connection is not None:
                connection.close()
            return

        if not force and process.is_alive() and connection is not None:
            try:
                connection.send(None)
            except (BrokenPipeError, EOFError, OSError):
                force = True
        if connection is not None:
            connection.close()

        process.join(timeout=1.0)
        if process.is_alive():
            process.terminate()
            process.join(timeout=2.0)
        if process.is_alive():
            process.kill()
            process.join(timeout=2.0)
        process.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, traceback) -> None:
        self.close()

    def __del__(self) -> None:
        try:
            self._stop(force=True)
        except Exception:
            pass
