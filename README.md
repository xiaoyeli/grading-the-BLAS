# Grading the BLAS

Measures the floating-point accuracy of a DGEMM implementation by sweeping over matrix dimensions and input scaling factors, computing componentwise errors against a reference, and producing a figure summarizing the results.

## Requirements

1. A shared or static library exposing `dgemm_` via the standard BLAS Fortran API
2. A C++17 compiler
3. Python 3 (for the virtual environment and matplotlib)

_All paths below are relative to the repository root._

---

## One-time setup

Install the Python plotting dependencies into the local virtual environment:

```bash
./scripts/install_plot_dependencies.sh
```

---

## Provide a config file for your DGEMM

Create a `.mk` file describing the library under test. See `config_files/dgemm_configs/` for examples (OpenBLAS, BLIS, cuBLAS).

Key variables:

| Variable | Purpose |
|---|---|
| `DGEMM_IMPL_NAME` | Label used in output paths |
| `LIB_PATH` / `LDFLAGS` | Path to the `.so` and linker flags |
| `CXX` | Compiler (`g++` or `nvcc` for cuBLAS) |
| `USE_CUDA` | Set to `1` when targeting cuBLAS |
| `OUTER_SCALING` / `INNER_SCALING` | Optional diagonal pre/post-scaling around the DGEMM call |

---

## Run

```bash
make run DGEMM_CONFIG_MK=path/to/your_dgemm_config.mk \
         TEST_CONFIG_MK=config_files/test_configs/2c_sweep.mk
```

This will:
1. Compile the test harness
2. Create a timestamped experiment directory `results/YYYYMMDD_HHMMSS/`
3. Copy both config files into that directory to document the parameters used to generate the data
4. Run the sweep and write raw binary results under `results/YYYYMMDD_HHMMSS/raw_data/`
5. Automatically generate `results/YYYYMMDD_HHMMSS/figure.png`

---

## Output structure

```
results/
└── YYYYMMDD_HHMMSS/
    ├── your_dgemm_config.mk       <- config snapshot
    ├── 2c_sweep.mk                <- config snapshot
    ├── figure.png                 <- generated automatically
    └── raw_data/
        └── <INPUT_TYPE>/
            └── scl<s>/
                └── m<m>_n<n>_k<k>/
                    └── <DGEMM_IMPL_NAME>/
                        ├── f_of_n.fp64
                        └── rel_errors.fp64
```

`f_of_n.fp64` contains the normalized error metric `f(n) = |C - C_ref| / (|A|·|B| · u)` where `u = 2^-53`. `rel_errors.fp64` contains componentwise relative errors. Both are column-major binary `double` arrays.

---

## Changelog

### 2026-06

**New behavior:**

- `make run` now runs the plot script automatically on completion.
- Config files are saved alongside results to document the parameters used to generate the data.

**Breaking changes for existing users:**

- **New output layout.** Results are no longer written to `results/raw_data/`. Each `make run` now creates a unique timestamped directory `results/YYYYMMDD_HHMMSS/` containing `raw_data/` and a copy of both config files. Existing result directories under `results/raw_data/` are unaffected but will not be picked up by the updated plot script.

- **Plot script interface changed.** `plot_max_values_across_dims_and_scales.py` renamed to `plot_sweep_results.py` and it now takes a single experiment directory as its argument instead of a glob of file paths:
  ```bash
  # old
  python3 scripts/plot_max_values_across_dims_and_scales.py results/raw_data/test2c/*/*/*/f_of_n.fp64
  # new
  python3 scripts/plot_sweep_results.py results/YYYYMMDD_HHMMSS
  ```
  The script searches for `f_of_n.fp64` files recursively under `<dir>/raw_data/` and writes `figure.png` as a sibling to `raw_data/`.