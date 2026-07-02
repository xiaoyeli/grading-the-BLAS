# Grading the BLAS

A floating-point accuracy benchmarking harness for DGEMM implementations. The harness sweeps over matrix dimensions and input scaling factors, computes error metrics against a high-precision reference, and produces figures summarizing the results.

## Requirements

1. A shared or static library exposing `dgemm_` via the standard BLAS Fortran API
2. A C++17 compiler
3. Python 3 (for data visualization)

## One-time setup

Install the Python plotting dependencies into the local virtual environment:

```bash
./scripts/install_plot_dependencies.sh
```

---

## Specifying a Test Configuration

### DGEMM config (`config_files/dgemm_configs/*.mk`)

Describes the implementation under test.

| Variable | Purpose |
|---|---|
| `DGEMM_IMPL_NAME` | Label used in output paths |
| `LIB_PATH` | Path to the library exposing the `dgemm_` symbol under test |
| `CXX` | Compiler (`g++`, or `nvcc` when `USE_CUDA=1`) |
| `CXX_FLAGS` / `LDFLAGS` | Compiler flags / linker flags |
| `USE_CUDA` | Set to `1` to wrap `cublasDgemm` instead of `dgemm_` |
| `RUN_PREFIX` | optional; string to be prepended to CL invocation of test binary, i.e., `${RUN_PREFIX} ./test1` |

See `config_files/dgemm_configs/` for ready-made examples (OpenBLAS, BLIS, cuBLAS, Strassen variants, Ozaki variants).

To add a new implementation: copy any existing file, set `DGEMM_IMPL_NAME`, and point `LIB_PATH`/`LDFLAGS` at your `.so`.

### Sweep config (`config_files/common.mk`)

Controls the parameter sweep shared across all four tests.

| Variable | Purpose |
|---|---|
| `LOG2_MIN_DIM` / `LOG2_MAX_DIM` | Parameterizes the sweep of dimensions according to pseudocode below |
| `LOG2_MIN_SCALING_FACTOR` / `LOG2_MAX_SCALING_FACTOR` / `LOG2_INC_SCALING_FACTOR` | Parameterizes the sweep of scaling factors according to pseudocode below |
| `REFERENCE_GEMM` | `1` = double-precision reference, `2` = long-double reference |
| `DUMP_A_B` | `1` = also write `A` and `B` to disk |
| `VERBOSE` | `1` = print `A`, `B`, and `C` to stdout |

```
for ( int d = LOG2_MIN_DIM; d <= LOG2_MAX_DIM; ++d ){
    // A and B have shape (2^d) x (2^d)
    ...
    for ( int s = LOG2_MIN_SCALING_FACTOR; s <= LOG2_MAX_SCALING_FACTOR; s += LOG2_INC_SCALING_FACTOR ){
        // entries in A and B span the orders of 2^-s to 2^s (note that LOG2_INC_SCALING_FACTOR must be >0)
        ...
        // execute DGEMM under test
        ...
}
```

---

## Running One or More Tests

```bash
# Run all tests
make grade DGEMM_CONFIG_MK=...

# Run individual tests
make run_test1 DGEMM_CONFIG_MK=...
make run_test2 DGEMM_CONFIG_MK=...
make run_test3 DGEMM_CONFIG_MK=...
make run_test4 DGEMM_CONFIG_MK=...
```

Figures are automatically written to `results/figures/<input_type>_<impl>.png`.

- **test1**: single heatmap — percentage of target values recovered (0–100%, `Reds_r` colormap)
- **test2 / test3 / test4**: line plot (one line per scale factor, log-log axes) + heatmap of the scalar max metric

Rerunning any target overwrites its output directory; no timestamped directories are created.

---

## Tests

The harness has four test binaries, each targeting a different class of difficult input.

> TODO: decribe each

---

## Output structure

```
results/raw_data/
├── test1c/
│   └── <DGEMM_IMPL_NAME>/scl<s>/m<m>_n<n>_k<k>/
│       └── values.fp64     # full C matrix (column-major float64)
│
├── test2c/
│   ├── C_ref/scl<s>/m<m>_n<n>_k<k>/
|   |   └── values.fp128            # full C reference matrix (column-major float128)
|   |
│   └── <DGEMM_IMPL_NAME>/scl<s>/m<m>_n<n>_k<k>/
|       ├── values.fp64             # full C matrix (column-major float64)
│       └── max_f_of_n.fp64         # scalar float64
│
├── test3/
│   ├── C_ref/scl<s>/m<m>_n<n>_k<k>/
|   |   └── values.fp128            # full C reference matrix (column-major float128)
|   |
│   └── <DGEMM_IMPL_NAME>/scl<s>/m<m>_n<n>_k<k>/
|       ├── values.fp64             # full C matrix (column-major float64)
│       └── max_f_of_n.fp64         # scalar float64
│
└── test2c_w_inner_permute/
    └── <DGEMM_IMPL_NAME>/scl<s>/m<m>_n<n>_k<k>/
        ├── values.fp64             # full C matrix (column-major float64)
        └── max_rel_error.fp64      # scalar float64
```

Shared reference matrices like `C_ref` are cached on first run and reused automatically.