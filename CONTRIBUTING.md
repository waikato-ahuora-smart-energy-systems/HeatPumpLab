# Contributing to HeatPumpLab

Thanks for your interest. Bug reports, questions, new cycle models and fixes are all welcome.

## Reporting a problem or asking for a feature

Open an [issue](https://github.com/waikato-ahuora-smart-energy-systems/HeatPumpLab/issues/new/choose). For a
wrong or failed solve, include the cycle, the inputs (or a case link from the web apps), the property backend
(CoolProp or REFPROP) and what you expected.

## Setting up

```sh
git clone https://github.com/waikato-ahuora-smart-energy-systems/HeatPumpLab.git
cd HeatPumpLab
uv sync --group dev                  # add --extra refprop if you have a REFPROP licence
uv run pytest                        # the full suite; REFPROP tests skip without REFPROP
uv run ruff check --select F HeatPumpLab tests
```

The JavaScript tests need [Node.js](https://nodejs.org) 20 or later on your PATH; they skip without it.

## Making a change

1. Branch from `develop`. `main` only receives pull requests from `develop`, and only when the tests pass.
2. Keep the change focused, and add or update tests in `tests/` for any behaviour you change. Cycle results
   are compared against reference values, so a change that moves a COP should say why in the pull request.
3. Results should not depend on the property backend beyond the documented differences: run the affected tests
   with CoolProp, and with `HEATPUMPLAB_PROPERTY_BACKEND=refprop` if you can.
4. Open a pull request into `develop`. CI lints, rebuilds the web apps and runs the suite on CoolProp.

### Web apps

The pages in `apps/steam` and `apps/dryer` are built from `apps/src`. Edit the sources, rebuild with
`uv run python apps/src/steam/build.py` or `uv run python apps/src/dryer/build.py`, and commit the sources and
the rebuilt page together; CI fails if they disagree. See [apps/src/README.md](apps/src/README.md). If you change
a cycle model in Python, change its JavaScript port to match: the tests check the two against each other.

### Adding a cycle

New cycles go in `HeatPumpLab/cycles/<Name>/` with a config dataclass, a model class with `solve()`, a result,
tests in `tests/`, a row in the README's cycle table and, where there is one, the reference it follows.

## Licence

HeatPumpLab is licensed under the [Apache License 2.0](LICENSE). By contributing, you agree that your
contribution is licensed under the same terms.
