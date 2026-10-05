# Installation

z7py is published on [PyPI](https://pypi.org/project/z7py/) and on
[prefix.dev](https://prefix.dev/channels/@allixender/geo). It requires
Python 3.10 or newer, NumPy (1.24 or newer) and Numba (0.60 or newer).

## pip

```bash
pip install z7py
```

## pixi

```bash
pixi project channel add https://prefix.dev/allixender/geo
pixi add z7py
```

## conda or mamba

```bash
conda install -c https://prefix.dev/allixender/geo z7py
```

## Not yet on conda-forge

`conda install -c conda-forge z7py` is the intended end state, but the package
is not published there yet. Until then use PyPI or the channel above.

## From source

```bash
git clone https://github.com/allixender/z7py.git
cd z7py
pip install -e ".[test]"
pytest -q
```

See [](development.md) for the Pixi-based development setup.
