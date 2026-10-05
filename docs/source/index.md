# z7py

z7py implements the Z7 indexing system for the IGEO7 (ISEA7H) aperture 7
hexagonal Discrete Global Grid System (DGGS) in Python. Z7 is based on the
Generalized Balanced Ternary (GBT) numeral system with an alternating
clockwise/counter-clockwise rotation pattern, as described in
[Kmoch et al. (2025)](https://doi.org/10.5194/agile-giss-6-32-2025).

Z7 indices are bit-packed into a single `uint64` (4 bits for the base cell,
20 x 3 bits for the resolution digits), which makes neighbour traversal,
parent/child navigation and monotonic range arithmetic pure integer
operations. The hot paths are JIT-compiled with Numba.

The package is deliberately thin: it depends only on NumPy and Numba. It needs
neither the DGGRID binary nor a geospatial stack, so it can be imported inside
Numba `nopython` kernels and on bare compute workers.

```python
from z7py import z7

raw = z7.encode_z7int(3, [1] * 10)        # base cell 3, 10 resolution digits
z7.index_to_z7string(raw)                 # '031111111111'
z7.get_neighbours(raw)                    # 6 neighbours as uint64
```

:::{note}
z7py is work in progress and currently an alpha release. The function names
and signatures can still change between minor versions.
:::

Find the source repository on [GitHub](https://github.com/allixender/z7py).
For running the DGGRID software itself from Python (grid generation, cell
geometries, coordinate transformations), see
[dggrid4py](https://dggrid4py.readthedocs.io/).

```{toctree}
:maxdepth: 2
:caption: User guide

installation
usage
resolutions
```

```{toctree}
:maxdepth: 2
:caption: Background

concepts
range_index
```

```{toctree}
:maxdepth: 2
:caption: Reference

api
changelog
development
references
```
