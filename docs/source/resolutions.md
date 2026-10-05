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
| 0 | 12 | 51,006,562.1724089 | 51,006,562,172,408.9 | 8,199.5003701 | 8,199,500.3701 |
| 1 | 72 | 7,286,651.7389156 | 7,286,651,738,915.6 | 3,053.2232428 | 3,053,223.2428 |
| 2 | 492 | 1,040,950.2484165 | 1,040,950,248,416.5 | 1,151.6430095 | 1,151,643.0095 |
| 3 | 3,432 | 148,707.1783452 | 148,707,178,345.2 | 435.1531492 | 435,153.1492 |
| 4 | 24,012 | 21,243.8826207 | 21,243,882,620.7 | 164.4655799 | 164,465.5799 |
| 5 | 168,072 | 3,034.8403744 | 3,034,840,374.4 | 62.1617764 | 62,161.7764 |
| 6 | 1,176,492 | 433.5486249 | 433,548,624.9 | 23.4949231 | 23,494.9231 |
| 7 | 8,235,432 | 61.9355178 | 61,935,517.8 | 8.8802451 | 8,880.2451 |
| 8 | 57,648,012 | 8.8479311 | 8,847,931.1 | 3.3564171 | 3,356.4171 |
| 9 | 403,536,072 | 1.2639902 | 1,263,990.2 | 1.2686064 | 1,268.6064 |
| 10 | 2,824,752,492 | 0.1805700 | 180,570.0 | 0.4794882 | 479.4882 |
| 11 | 19,773,267,432 | 0.0257957 | 25,795.7 | 0.1812295 | 181.2295 |
| 12 | 138,412,872,012 | 0.0036851 | 3,685.1 | 0.0684983 | 68.4983 |
| 13 | 968,890,104,072 | 0.0005264 | 526.4 | 0.0258899 | 25.8899 |
| 14 | 6,782,230,728,492 | 0.0000752 | 75.2 | 0.0097855 | 9.7855 |
| 15 | 47,475,615,099,432 | 0.0000107 | 10.7 | 0.0036986 | 3.6986 |
| 16 | 332,329,305,696,012 | 0.0000015 | 1.5 | 0.0013979 | 1.3979 |
| 17 | 2,326,305,139,872,072 | 0.0000002 | 0.2 | 0.0005284 | 0.5284 |
| 18 | 16,284,135,979,104,492 | 0.0000000 | 0.0 | 0.0001997 | 0.1997 |
| 19 | 113,988,951,853,731,432 | 0.0000000 | 0.0 | 0.0000755 | 0.0755 |
| 20 | 797,922,662,976,120,012 | 0.0000000 | 0.0 | 0.0000285 | 0.0285 |

:::{note}
The cell areas are stored with a fixed number of decimals (7 for km², 1 for
m²). For this reason the areas from resolution 18 onwards are returned as 0.0,
although the cells have a non-zero area. Use the CLS values for these
resolutions.
:::
