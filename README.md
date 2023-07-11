# Transcritical-recuperative-heat-hump-model

High temperature transcritical-recuperative heat pump model
By Tim Walmsley and Keri Adamson, University of Waikato, 2021-2023 (tim.walmsley@waikato.ac.nz)


## Software requirements

- [Anaconda](https://www.anaconda.com/products/individual).
- [Refprop](https://www.nist.gov/srd/refprop)
## Create conda environment

Create a `HTHP_env` conda environment with the python packages specified in [environment.yml](environment.yml).
```bash
conda env create -f environment.yml
```
This step is only required once. If it becomes necessary to remove the `HTHP_env` environment (e.g. to start over), this can be done using:
```bash
conda env remove -n HTHP_env
```

## Installation of this package into `HTHP_env` env

```bash
conda activate HTHP_env
pip install .
```

## Update of `HTHP_env` env

This step should only be needed when a new python package update is available or a new package is needed to be installed

```bash
conda env update --prune
```


## Additional tools

- autopep8 helps to format your code in conformance with the pep8 style guide. Usage example:
```
autopep8 -i classes\basicM_v5.py
```