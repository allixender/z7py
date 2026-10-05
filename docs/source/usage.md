# Usage

All functions live in the `z7py.z7` module and are also re-exported from the
top-level `z7py` package. The examples below use one cell throughout: base
cell 3 with ten resolution digits of value 1.

```python
import numpy as np
from z7py import z7
```

## Three representations of a Z7 index

A Z7 index identifies one cell by its base cell (0 to 11) and up to 20
resolution digits (0 to 6). z7py handles the three representations of an index
which are in use:

- the packed integer (`np.uint64`), which all Numba kernels operate on;
- the hex string (16 characters), which is the DGGRID `Z7` output format;
- the digit string (two characters for the base cell, then one character per
  resolution digit), which DGGRID calls `Z7_STRING`.

```python
raw = z7.encode_z7int(3, [1] * 10)
raw                                   # np.uint64(3623467586827583487)

z7.z7int_to_z7hex(raw)                # '324924927fffffff'
z7.index_to_z7string(raw)             # '031111111111'
```

The conversions work in every direction:

```python
z7.z7hex_to_z7int("324924927fffffff") == raw        # True
z7.z7hex_to_z7string("324924927fffffff")            # '031111111111'
z7.z7string_to_index("031111111111") == raw         # True
z7.encode_z7hex_index(3, [1] * 10)                  # '324924927fffffff'
```

`decode_z7int` and `decode_z7hex_index` return the base cell and all 20 digit
slots. Unused slots carry the padding value 7 (see [](concepts.md)):

```python
z7.decode_z7hex_index("324924927fffffff")
# (3, [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7])
```

## Inspecting an index

```python
z7.get_resolution(raw)                # 10
z7.get_base_cell(raw)                 # 3
z7.get_digit(raw, 1)                  # 1  (positions are 1-based)
z7.get_digits(raw)                    # uint8 array of all 20 digit slots
```

For the string representations the resolution and the position of a cell
inside its parent are available without decoding to an integer first. The
last element of the returned tuple states whether the cell is the centre child
(digit 0) of its parent:

```python
z7.get_z7string_resolution("031111111111")   # 10
z7.get_z7string_local_pos("031111111111")    # ('03111111111', '1', False)
z7.get_z7hex_resolution("324924927fffffff")  # 10
```

## Parents

`get_parent` moves one level up by default, or to any coarser resolution:

```python
z7.index_to_z7string(z7.get_parent(raw))      # '03111111111'  (resolution 9)
z7.index_to_z7string(z7.get_parent(raw, 5))   # '0311111'      (resolution 5)
```

## Neighbours

`get_neighbours` returns the six neighbours of a cell as a `uint64` array.
The function resolves carries across parent cells and across base cells:

```python
[z7.index_to_z7string(n) for n in z7.get_neighbours(raw)]
# ['024242424242', '031111111113', '024242424264',
#  '031111111115', '024242424246', '031111111110']
```

At every resolution, the twelve cells centred on the icosahedron vertices are
pentagons and have only five neighbours. The missing neighbour is returned as the maximum
`uint64` value:

```python
INVALID = np.uint64(np.iinfo(np.uint64).max)
pentagon = z7.encode_z7int(0, [0, 0, 0])
[z7.index_to_z7string(n) if n != INVALID else None
 for n in z7.get_neighbours(pentagon)]
# ['00001', None, '00003', '00004', '00005', '00006']
```

The neighbouring base cells of a base cell are available directly:

```python
z7.get_base_cell_neighbours(3)        # [2, 0, 8, 7, 4]
```

`get_neighbour` is the low-level building block that `get_neighbours` uses
internally. It applies the GBT addition for one direction inside a single base
cell and returns the result together with the remaining carry. A carry other
than 0 means that the neighbour lies in another base cell, which
`get_neighbour` does not resolve; use `get_neighbours` for complete results.

## Monotonic integer mapping

`z7_to_monotonic_int` maps an index at a given resolution onto a dense number
line, on which the seven children of a parent are consecutive integers. This
mapping is the basis for the range compression described in
[](range_index.md).

```python
mono = z7.z7_to_monotonic_int(raw, 10)           # 894504955
z7.monotonic_int_to_z7(mono, 10) == raw          # True

siblings = [z7.encode_z7int(3, [1] * 9 + [d]) for d in range(7)]
[int(z7.z7_to_monotonic_int(s, 10)) for s in siblings]
# [894504954, 894504955, 894504956, 894504957,
#  894504958, 894504959, 894504960]
```

## Resolution statistics

The number of cells, the cell area and the characteristic length scale (CLS)
per resolution are available as lookups (see [](resolutions.md) for the full
table):

```python
z7.get_num_cells(10)                  # 2824752492
z7.get_cell_area_km2(10)              # 0.1805700228709557
z7.get_cell_area_m2(20)               # 0.0006392419283112329
z7.get_cls_m(10)                      # 479.4882
z7.get_resolution_stats(10)
# {'num_cells': 2824752492, 'area_km2': 0.1805700228709557,
#  'area_m2': 180570.02287095567, 'cls_km': 0.4794882, 'cls_m': 479.4882}
```

The `find_resolution_by_*` functions select the resolution for a target value.
With `prefer="closest"` (the default) the resolution with the smallest
absolute difference is returned. `prefer="larger"` and `prefer="smaller"`
restrict the search to resolutions whose value is not below, or not above, the
target:

```python
z7.find_resolution_by_cls_m(500.0)                     # 10 (CLS 479 m)
z7.find_resolution_by_cls_m(500.0, prefer="larger")    # 9  (CLS 1269 m)
z7.find_resolution_by_area_m2(1e6)                     # 9
z7.find_resolution_by_num_cells(1_000_000)             # 6
```

## Authalic and geodetic latitudes

`z7py.latitudes` converts between geodetic latitudes on the WGS84 ellipsoid
and authalic latitudes (the latitudes on a sphere with the same surface area
as the ellipsoid), in degrees by default:

```python
from z7py import geodetic_to_authalic, authalic_to_geodetic

xi = geodetic_to_authalic(58.38)      # 58.2653451741872
authalic_to_geodetic(xi)              # 58.38000000000001
```

## Use inside Numba kernels

The index operations are compiled with `@njit(cache=True)` and can be called
from other `nopython` functions, for example to convert a whole chunk of
indices without leaving compiled code:

```python
import numba as nb
from z7py import z7_to_monotonic_int

@nb.njit
def to_monotonic(chunk, resolution):
    out = np.empty(len(chunk), dtype=np.uint64)
    for i in range(len(chunk)):
        out[i] = z7_to_monotonic_int(chunk[i], resolution)
    return out

to_monotonic(z7.get_neighbours(raw), 10)
# array([741497528, 894504957, 741497544, 894504959, 741497532, 894504954],
#       dtype=uint64)
```
