# z7py

<img src="https://raw.githubusercontent.com/allixender/z7py/main/images/igeo7_logo.svg" width="100" alt="IGEO7 logo">

[![PyPI](https://img.shields.io/pypi/v/z7py?label=PyPI)](https://pypi.org/project/z7py/)
[![prefix.dev](https://img.shields.io/badge/dynamic/json?url=https%3A%2F%2Frepo.prefix.dev%2Fallixender%2Fgeo%2Fnoarch%2Frepodata.json&query=%24%5B%27packages.conda%27%5D%5B*%5D.version&prefix=v&label=prefix.dev&color=F7CC49)](https://prefix.dev/channels/@allixender/geo/packages/z7py)
[![Documentation Status](https://readthedocs.org/projects/z7py/badge/?version=latest)](https://z7py.readthedocs.io/en/latest/)

z7py implements the Z7 indexing system for the IGEO7 (ISEA7H) aperture 7
hexagonal Discrete Global Grid System (DGGS) in Python. Z7 is based on the
Generalized Balanced Ternary (GBT) numeral system with an alternating
clockwise/counter-clockwise rotation pattern, as described in
[Kmoch et al. (2025)](https://doi.org/10.5194/agile-giss-6-32-2025).

Z7 indices are bit-packed into a single `uint64`, which makes neighbour
traversal, parent/child navigation and monotonic range arithmetic pure integer
operations. z7py is deliberately thin: it depends only on **NumPy** and
**Numba**. It needs neither the DGGRID binary nor a geospatial stack, so it can
be imported inside Numba `nopython` kernels and on bare compute workers.

This is work in progress; the function names and signatures can still change.

## Installation

```bash
pip install z7py
```

With pixi, or with conda/mamba, from the
[prefix.dev channel](https://prefix.dev/channels/@allixender/geo):

```bash
pixi project channel add https://prefix.dev/allixender/geo
pixi add z7py

conda install -c https://prefix.dev/allixender/geo z7py
```

## Quick start

```python
import numpy as np
from z7py import z7

raw = z7.encode_z7int(3, [1] * 10)        # base cell 3, 10 resolution digits
print(z7.get_resolution(raw))             # 10
print(z7.index_to_z7string(raw))          # '031111111111'

# neighbour traversal is pure integer arithmetic
print(z7.get_neighbours(raw))            # 6 neighbours, UINT64_MAX where invalid

# monotonic mapping preserves sort order -> range compression
mono = z7.z7_to_monotonic_int(raw, 10)
assert z7.monotonic_int_to_z7(mono, 10) == raw
```

## Documentation

The documentation is hosted at <https://z7py.readthedocs.io/>:

- [Usage](https://z7py.readthedocs.io/en/latest/usage.html): index representations, parents, neighbours, monotonic mapping, resolution lookups;
- [Resolutions](https://z7py.readthedocs.io/en/latest/resolutions.html): number of cells, cell area and characteristic length scale per resolution;
- [Z7 indexing concepts](https://z7py.readthedocs.io/en/latest/concepts.html): GBT, rotation pattern, bit layout and neighbour traversal;
- [Sorted range index for Z7](https://z7py.readthedocs.io/en/latest/range_index.html): design notes for an Xarray/Zarr range index on top of z7py;
- [API reference](https://z7py.readthedocs.io/en/latest/api.html);
- [Development](https://z7py.readthedocs.io/en/latest/development.html): Pixi setup, tests, building the docs and releasing.

## Development

```bash
pixi run test
```

See the [development page](https://z7py.readthedocs.io/en/latest/development.html)
for the details.

## Acknowledgements

- DGGRID, 2019-2026 Kevin Sahr <sahrk@sou.edu>, https://github.com/sahrk/DGGRID/
- Z7, 2025-2026 Weston James Renoud, https://github.com/wrenoud/Z7/
- 2025-2026 Javier Jimenez Shaw, https://github.com/jjimenezshaw

## License

Apache License 2.0, see [LICENSE](https://github.com/allixender/z7py/blob/main/LICENSE) and [NOTICE](https://github.com/allixender/z7py/blob/main/NOTICE).

## Citing

If you use z7py in academic work, please cite the IGEO7 paper and the software
itself; see [CITATION.cff](https://github.com/allixender/z7py/blob/main/CITATION.cff).

Kmoch, A., Sahr, K., Chan, W. T., & Uuemaa, E. (2025). IGEO7: A new
hierarchically indexed hexagonal equal-area discrete global grid system.
_AGILE: GIScience Series, 6_(32). https://doi.org/10.5194/agile-giss-6-32-2025

The full reference list is given in the
[documentation](https://z7py.readthedocs.io/en/latest/references.html).
