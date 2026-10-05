# Resolutions

IGEO7 has 12 base cells at resolution 0, and every refinement step multiplies
the number of cells by approximately 7 (aperture 7). The number of cells at
resolution *r* is 10 x 7^r + 2. A Z7 index holds up to 20 resolution digits,
so z7py covers the resolutions 0 to 20.

The table lists the values that are stored in `z7py.z7.RESOLUTION_STATS` and
returned by `get_resolution_stats`, `get_num_cells`, `get_cell_area_km2`,
`get_cell_area_m2`, `get_cls_km` and `get_cls_m`. The characteristic length
scale (CLS) is the diameter of a spherical cap of the same area as a cell, as
defined in DGGRID.

| Resolution | Number of cells | Cell area (km²) | Cell area (m²) | CLS (km) | CLS (m) |
|---:|---:|---:|---:|---:|---:|
| 0 | 12 | 51,006,562 | 51,006,562,172,409 | 8,199.5003701 | 8,199,500.3701 |
| 1 | 72 | 7,286,652 | 7,286,651,738,916 | 3,053.2232428 | 3,053,223.2428 |
| 2 | 492 | 1,040,950 | 1,040,950,248,417 | 1,151.6430095 | 1,151,643.0095 |
| 3 | 3,432 | 148,707 | 148,707,178,345 | 435.1531492 | 435,153.1492 |
| 4 | 24,012 | 21,244 | 21,243,882,621 | 164.4655799 | 164,465.5799 |
| 5 | 168,072 | 3,035 | 3,034,840,374 | 62.1617764 | 62,161.7764 |
| 6 | 1,176,492 | 433.549 | 433,548,625 | 23.4949231 | 23,494.9231 |
| 7 | 8,235,432 | 61.936 | 61,935,518 | 8.8802451 | 8,880.2451 |
| 8 | 57,648,012 | 8.848 | 8,847,931 | 3.3564171 | 3,356.4171 |
| 9 | 403,536,072 | 1.264 | 1,263,990 | 1.2686064 | 1,268.6064 |
| 10 | 2,824,752,492 | 0.181 | 180,570 | 0.4794882 | 479.4882 |
| 11 | 19,773,267,432 | 0.0258 | 25,796 | 0.1812295 | 181.2295 |
| 12 | 138,412,872,012 | 0.00369 | 3,685 | 0.0684983 | 68.4983 |
| 13 | 968,890,104,072 | 5.264e-04 | 526.443 | 0.0258899 | 25.8899 |
| 14 | 6,782,230,728,492 | 7.521e-05 | 75.206 | 0.0097855 | 9.7855 |
| 15 | 47,475,615,099,432 | 1.074e-05 | 10.744 | 0.0036986 | 3.6986 |
| 16 | 332,329,305,696,012 | 1.535e-06 | 1.535 | 0.0013979 | 1.3979 |
| 17 | 2,326,305,139,872,072 | 2.193e-07 | 0.219 | 0.0005284 | 0.5284 |
| 18 | 16,284,135,979,104,492 | 3.132e-08 | 0.0313 | 0.0001997 | 0.1997 |
| 19 | 113,988,951,853,731,432 | 4.475e-09 | 0.00447 | 0.0000755 | 0.0755 |
| 20 | 797,922,662,976,120,012 | 6.392e-10 | 6.392e-04 | 0.0000285 | 0.0285 |

The number of cells and the CLS are the values reported by DGGRID. The cell
area is the area of a hexagonal cell; the twelve pentagonal cells have 5/6 of
this area. Because every refinement step divides the cell area by exactly 7,
z7py derives the areas from the resolution 0 value (51,006,562.1724089 km², a
tenth of the surface of the WGS84 authalic sphere) and returns them in full
floating point precision. The areas are rounded in the table above for
display only.

:::{note}
Up to version 0.1.0 the areas were stored as printed by DGGRID (7 decimals for
km², 1 decimal for m²), so that they were returned as 0.0 from resolution 18
onwards. The values of version 0.1.1 agree with the previous ones at that
printed precision.
:::
