# Development

## Development with Pixi

This project uses [Pixi](https://pixi.sh), a modern, high-performance package manager built on the Conda ecosystem. Pixi is particularly well-suited for scientific projects because it provides strict reproducibility through lockfiles and handles complex C/C++ dependencies (like those often found in GIS and numerical libraries) better than standard `pip`.

### Getting started

1.  **Install Pixi:**
    ```bash
    curl -fsSL https://pixi.sh/install.sh | bash
    ```
    *(For other platforms, see the [installation guide](https://pixi.sh/latest/#installation)).*

2.  **Setup and Run Tests:**
    Pixi automatically manages the environment for you. You don't need to `pip install` anything.
    ```bash
    pixi run test
    ```

3.  **Enter the Development Environment:**
    If you want to run a script or start a REPL within the project environment:
    ```bash
    pixi shell
    python
    ```

### Jupyter integration

For interactive science and experimentation, you can easily use this Pixi environment with Jupyter:

1.  **Add Jupyter to the project:**
    ```bash
    pixi add jupyterlab ipykernel
    ```

2.  **Run Jupyter Lab directly:**
    ```bash
    pixi run jupyter lab
    ```

3.  **Using with VS Code or external Jupyter:**
    If you prefer using VS Code, simply select the Python interpreter located in `.pixi/envs/default/bin/python`. VS Code will automatically recognize the environment and its packages.

### Why Pixi?

-   **Reproducibility:** The `pixi.lock` file ensures that every collaborator uses the exact same versions of all dependencies, including Python itself and system-level libraries.
-   **No Activation Needed:** Unlike `conda` or `venv`, you don't need to manually activate environments to run tasks. `pixi run <task>` handles it instantly.
-   **Unified Conda & Pip:** Pixi can install and manage dependencies from both Conda channels (like `conda-forge`) and PyPI (`pip`) in the same environment. This solves the "missing package" problem that often plagues tools like `poetry` or `micromamba`.
-   **Project-Local Environments:** Unlike Conda/Micromamba which use global named environments, Pixi stores the environment inside the project directory (`.pixi/`). This makes it much easier to run isolated experiments with different package versions without polluting your system or forgetting which "env" was for which project.
-   **Multi-language:** It handles Python, R, C++, and more, making it ideal for projects that bridge high-level analysis and low-level kernels.

### The lock file

`pixi.lock` is tracked in the repository. It pins the `default` and the `docs`
environment for the four platforms listed in `pixi.toml` (`linux-64`,
`osx-64`, `osx-arm64` and `win-64`). The lock file is solved for all four
platforms at once, independent of the machine it is generated on. After a
change of the dependencies in `pixi.toml` or `pyproject.toml`, update and
commit it:

```bash
pixi lock
```

A version bump of z7py itself does not change the lock file.

## Development without Pixi

A plain virtual environment works as well:

```bash
pip install -e ".[test]"
pytest -q
```

## Building the documentation

The documentation is written in MyST Markdown, built with Sphinx and hosted on
[Read the Docs](https://z7py.readthedocs.io/). To build it locally:

```bash
pixi run -e docs docs
```

or, without Pixi:

```bash
pip install -e ".[docs]"
sphinx-build -b html docs/source docs/build/html
```

The result is written to `docs/build/html/index.html`. The API reference is
generated from the docstrings in `z7py/z7.py` and `z7py/latitudes.py`. The
build sets `NUMBA_DISABLE_JIT=1`, because Sphinx cannot read signatures and
docstrings from compiled Numba dispatchers.

## Releasing

The version number is kept in four files, which have to be updated together:

- `z7py/__init__.py` (`__version__`, the source for the PyPI build);
- `pixi.toml` (`version`);
- `recipe.yaml` (`context.version`, the source for the conda build);
- `CITATION.cff` (`version`).

Publishing a GitHub release with a tag that matches the version (for example
`v0.1.1`) runs `.github/workflows/publish.yml`. This workflow builds and
uploads the conda package to the prefix.dev channel (this requires the
`PREFIX_API_KEY` repository secret) and the sdist and wheel to PyPI (via
Trusted Publishing). The workflow refuses a release whose tag does not match
`__version__`. In a manual run of the workflow the targets can be selected
individually, including a dry run against TestPyPI:

```bash
gh workflow run publish.yml -f pypi=testpypi
gh workflow run publish.yml -f pypi=pypi -f conda=true
```
