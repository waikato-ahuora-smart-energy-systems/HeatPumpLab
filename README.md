# Recuperative Heat Pump Model

High temperature recuperative heat pump model
By Tim Walmsley and Keri Adamson, University of Waikato, 2021-2025 (tim.walmsley@waikato.ac.nz)


## Software requirements

- [Anaconda](https://www.anaconda.com/products/individual).
- [Refprop](https://www.nist.gov/srd/refprop)
## Create conda environment

Create a `Ahuora_env` conda environment with the python packages specified in [environment.yml](environment.yml).
```bash
conda env create -f environment.yml
```
This step is only required once. If it becomes necessary to remove the `Ahuora_env` environment (e.g. to start over), this can be done using:
```bash
conda env remove -n Ahuora_env
```

## Installation of this package into `Ahuora_env` env

```bash
conda activate Ahuora_env
pip install .
```

## Update of `Ahuora_env` env

This step should only be needed when a new python package update is available or a new package is needed to be installed

```bash
conda env update --prune
```


## Additional tools

- autopep8 helps to format your code in conformance with the pep8 style guide. Usage example:
```
autopep8 -i classes\basicM_v5.py
```
