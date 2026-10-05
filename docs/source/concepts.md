# Z7 indexing concepts

Z7 is an indexing system for aperture 7 hexagonal discrete global grid
systems, as coined in [Kmoch et al. (2025)][kmoch2025]. It is based on the
Generalized Balanced Ternary (GBT) numeral system described in
[Lucas, Gibson (1982)][lucas1982], [van Roessel (1988)][vanRoessel1988], and
[Wikipedia (2025)][wikipedia2025]. z7py explores neighbour traversal directly
in the Z7 domain, without generating cell geometries.

## Rotation pattern

Lucas and Gibson (1982) gave examples of GBT with exclusively counter-clockwise (CCW) rotation. The disadvantage of this
approach is that the orientation of the hexagons at the different hierarchy levels quickly diverge from each other.
[White et al. (1992)][white1992] proposed an alternating rotation direction (CW, CCW, CW, ...) to maintain an alignment
of hexagon
orientations across hierarchy levels. Kmoch et al. (2025) adopted this alternating rotation pattern for Z7 assigning CW
to odd resolutions and CCW to even resolutions.

## GBT addition

The single-digit GBT addition tables for both rotation directions are
illustrated in the README of the [Z7 C++ library](https://github.com/wrenoud/Z7/)
by Weston Renoud. Figures for this documentation are still to be added.

## Bit layout

The `z7py` implementation uses a bit-packed 64-bit representation to enable extremely fast neighbor traversal and range arithmetic.

Indices are stored as `uint64` with the following structure:
- **Bits 63–60:** Base Cell ID (0–11).
- **Bits 59–0:** 20 refinement digits, each occupying 3 bits.
- **Value 7:** A special padding digit indicating the end of the identifier.

Digit 1 occupies bits 59–57, digit 2 occupies bits 56–54, and so on down to
digit 20 in bits 2–0. This layout allows for efficient bit-shifting to extract
digits and maintains lexicographical sorting consistency.

As an example, base cell 3 with ten resolution digits of value 1 is packed as
follows (the hex string is the same value written in base 16):

```text
base cell   digits 1-10 (3 bits each)        digits 11-20 (padding, value 7)
0011        001 001 001 001 001 ... 001      111 111 111 ... 111

hex    0x324924927fffffff
uint64 3623467586827583487
```

## Neighbor traversal and cascading carries

Neighbor finding is implemented using GBT arithmetic tables (`_GBT_CW_*` and `_GBT_CCW_*`) compiled with Numba.
1. **Alternating Rotations:** The implementation automatically switches between Clockwise (CW) and Counter-Clockwise (CCW) addition tables based on whether the resolution level is odd or even.
2. **Cascading Carries:** Moving across a cell boundary often triggers a "carry" to the parent level. The algorithm recursively ripples these carries up the hierarchy (similar to carrying a 1 in decimal addition) until the carry is 0 or a base-cell boundary is reached.
3. **Bitwise Efficiency:** All digit extractions and updates use raw bit-shifts (`>>`) and masks (`&`), ensuring O(1) performance for single-level operations and avoiding expensive string or list allocations.

## Numba implementation choices

In `z7py`, several choices maximize performance:
- **`@njit(cache=True)`:** Compiles the indexing logic to machine code, bypassing the Python Global Interpreter Lock (GIL) and allowing Dask to execute the range calculation across multiple threads efficiently.
- **Lookup Tables in Memory:** Constants like `_BASE_CELL_NEIGHBOURS_RAW` are stored as contiguous numpy arrays, allowing the Numba compiler to optimize memory access patterns during index traversal.
- **Static Typing:** By forcing `np.uint64` and `np.uint8`, we avoid expensive type-checking and boxing/unboxing overhead in the hot path of the RangeIndex creation.

The full reference list is given in [](references.md).

[wikipedia2025]: https://en.wikipedia.org/wiki/Generalized_balanced_ternary

[kmoch2025]: https://agile-giss.copernicus.org/articles/6/32/2025/agile-giss-6-32-2025.pdf

[vanRoessel1988]: https://web.archive.org/web/20250523045258/https://www.asprs.org/wp-content/uploads/pers/1988journal/nov/1988_nov_1565-1570.pdf

[white1992]: https://doi.org/10.1559/152304092783786636

[lucas1982]: https://apps.dtic.mil/sti/tr/pdf/ADA125710.pdf
