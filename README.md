# Grading the BLAS (Quick Start)

REQUIRED:
1. a path to a static or dynamic library that exposes the symbol `dgemm_` using the standard BLAS API.
2. Python 3 or higher (for matplotlib dependency)

__Note: all paths in this README are relative to the root of this repository.__

## Provide a mk file describing your DGEMM
See examples in `./config_files/dgemm_configs`.

## Run the test
```
make run DGEMM_CONFIG_MK=path/to/your_dgemm_config.mk TEST_CONFIG_MK=./config_files/test_configs/2c_sweep.mk
```

## Plot the results

First, a one-time installation of the matplotlib dependency in a python virtual environment
```
./scripts/install_plot_dependencies.sh
```

Make sure the environment is active
```
source ./scripts/plot_env.sh
```

Plot
```
python3 scripts/plot_max_values_across_dims_and_scales.py results/raw_data/test2c/*/*/*/f_of_n.fp64
```

Find output png written to `./results`