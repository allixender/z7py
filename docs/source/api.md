# API reference

All functions of `z7py.z7` listed below are also importable from the top-level
`z7py` package, as are `geodetic_to_authalic`, `authalic_to_geodetic`,
`WGS84_A` and `WGS84_F` from `z7py.latitudes`.

Functions marked as compiled in the source (`@nb.njit(cache=True)`) expect
`np.uint64` indices and can be called from other Numba `nopython` functions.

## z7py.z7

```{eval-rst}
.. currentmodule:: z7py.z7
```

### Encoding and conversion

```{eval-rst}
.. autofunction:: encode_z7int
.. autofunction:: decode_z7int
.. autofunction:: encode_z7hex_index
.. autofunction:: decode_z7hex_index
.. autofunction:: z7int_to_z7hex
.. autofunction:: z7hex_to_z7int
.. autofunction:: z7hex_to_z7string
.. autofunction:: z7string_to_index
.. autofunction:: index_to_z7string
```

### Index inspection

```{eval-rst}
.. autofunction:: get_base_cell
.. autofunction:: get_digit
.. autofunction:: get_digits
.. autofunction:: get_resolution
.. autofunction:: first_non_zero
.. autofunction:: get_z7hex_resolution
.. autofunction:: get_z7hex_local_pos
.. autofunction:: get_z7string_resolution
.. autofunction:: get_z7string_local_pos
```

### Hierarchy

```{eval-rst}
.. autofunction:: get_parent
.. autofunction:: get_parent_at
```

### Neighbour traversal

```{eval-rst}
.. autofunction:: get_neighbours
.. autofunction:: get_neighbour
.. autofunction:: get_base_cell_neighbours
.. autofunction:: get_base_cell_neighbour
.. autofunction:: neighbour_addition_cw
.. autofunction:: neighbour_addition_ccw
.. autofunction:: neighbour_addition_ccw_mod
```

### Monotonic integer mapping

```{eval-rst}
.. autofunction:: z7_to_monotonic_int
.. autofunction:: monotonic_int_to_z7
```

### Resolution statistics

```{eval-rst}
.. py:data:: RESOLUTION_STATS
   :type: dict[int, dict]

   Dictionary of the grid statistics per resolution (0 to 20) with the keys
   ``num_cells``, ``area_km2``, ``area_m2``, ``cls_km`` and ``cls_m``.

.. autofunction:: get_resolution_stats
.. autofunction:: get_num_cells
.. autofunction:: get_cell_area_m2
.. autofunction:: get_cell_area_km2
.. autofunction:: get_cls_m
.. autofunction:: get_cls_km
.. autofunction:: find_resolution_by_value
.. autofunction:: find_resolution_by_cls_m
.. autofunction:: find_resolution_by_area_m2
.. autofunction:: find_resolution_by_num_cells
```

## z7py.latitudes

```{eval-rst}
.. automodule:: z7py.latitudes
   :members: geodetic_to_authalic, authalic_to_geodetic, geodetic_to_authalic_custom, authalic_to_geodetic_custom, compute_fourier_coeffs, third_flattening, normalized_meridian_arc_unit, horner, fourier_sin

.. py:data:: z7py.latitudes.WGS84_A
   :value: 6378137.0

   Semi-major axis of the WGS84 ellipsoid (in metres).

.. py:data:: z7py.latitudes.WGS84_F
   :value: 1.0 / 298.257223563

   Flattening of the WGS84 ellipsoid.
```
