# Configuration file for the Sphinx documentation builder.

import os

# Numba dispatchers are opaque to autodoc. With the JIT disabled, nb.njit
# returns the plain Python function, so signatures and docstrings are found.
os.environ["NUMBA_DISABLE_JIT"] = "1"

import z7py  # noqa: E402

# -- Project information

project = "z7py"
copyright = "2025-2026, Alexander Kmoch and z7py contributors"
author = "Alexander Kmoch"

release = z7py.__version__
version = release

# -- General configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "myst_parser",
]

myst_enable_extensions = ["dollarmath", "amsmath", "colon_fence"]
myst_heading_anchors = 3

autodoc_member_order = "bysource"

intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "numba": ("https://numba.readthedocs.io/en/stable/", None),
    "dggrid4py": ("https://dggrid4py.readthedocs.io/en/latest/", None),
}
intersphinx_disabled_domains = ["std"]

# -- Options for HTML output

html_theme = "sphinx_rtd_theme"
html_logo = "../../images/igeo7_logo.svg"
html_title = f"z7py {release}"
