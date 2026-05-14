import re
import numpy as np
from pathlib import Path

def load_matrix(path: str) -> np.ndarray:
    path = Path(path)
    match = re.search(r'm(\d+)_n(\d+)', str(path))
    if not match:
        raise ValueError(f"Could not extract m and n from path: {path}")
    m, n = int(match.group(1)), int(match.group(2))
    if path.suffix == ".fp64":
        data = np.fromfile(path, dtype=np.float64)
    elif path.suffix == ".fp128":
        data = np.fromfile(path, dtype=np.float128)
    else:
        raise ValueError(f"file extension is not one of .fp64 or .fp128: {path}")
    if data.size != m * n:
        raise ValueError(f"Expected {m*n} elements, got {data.size}")
    return data.reshape((m, n), order='F')

def parse_path_metadata(path: str) -> dict:
    """Extract scale factor s, dimension d (assuming m == n == k), impl name, and value descriptor from path."""
    path_str = str(path)
    scl_match = re.search(r'scl(\d+)', path_str)
    dim_match = re.search(r'm(\d+)_n(\d+)_k(\d+)', path_str)
    if not (scl_match and dim_match):
        raise ValueError(f"Could not parse scl/m/n/k from path: {path}")
    s = int(scl_match.group(1))
    m, n, k = (int(dim_match.group(i)) for i in (1, 2, 3))
    if not (m == n == k):
        raise ValueError(f"Expected m == n == k, got m={m}, n={n}, k={k} in {path}")
    impl_name = Path(path).parent.name
    value_descriptor = Path(path).stem
    return {"s": s, "d": m, "impl": impl_name, "value_descriptor": value_descriptor}