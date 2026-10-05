# Changelog

## 0.1.1

- Cell areas (`get_cell_area_m2`, `get_cell_area_km2`, `RESOLUTION_STATS`) are
  derived from the resolution 0 area and the aperture (area / 7^resolution)
  and are returned in full floating point precision. Before, the areas were
  stored as printed by DGGRID and were 0.0 from resolution 18 onwards. The new
  values agree with the previous ones at the printed precision (7 decimals
  for km², 1 decimal for m²), but they are not bit-identical: for example
  `get_cell_area_m2(14)` returns 75.206... instead of 75.2.
- The internal lookup tables are no longer keyword arguments of
  `get_neighbours`, `get_neighbour`, `neighbour_addition_cw`,
  `neighbour_addition_ccw` and `neighbour_addition_ccw_mod` (`bcn`, `excl`,
  `rotations_arr`, `pole0`, `cw0`, `cw1`, `ccw0`, `ccw1`, `mod7`). The
  functions read the tables as module constants. The results are unchanged,
  and `get_neighbours` runs about twice as fast.
- `get_parent_at` is exported from the top-level `z7py` package.
- All public functions have docstrings.
- Documentation on [Read the Docs](https://z7py.readthedocs.io/), and a
  shorter README.
- `pixi.lock` is tracked in the repository.

## 0.1.0

- First release as a standalone package.
