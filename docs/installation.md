# Installation

`kronrod` can be installed from PyPI, GitHub, or built from source.

---

## Prerequisites

- **Python 3.10+**
- **NumPy** (installed automatically as a dependency)

For source builds you additionally need:

- A C compiler (`gcc` or `clang`)
- `meson` and `meson-python` build system
- `numpy` (for `f2py` compilation)

## [PyPI](https://pypi.org/project/kronrod)

For using the PyPI package in your project, add the following to your
configuration file:

=== "pyproject.toml"

    ```toml
    [project]
    dependencies = [
        "kronrod"
    ]
    ```

=== "requirements.txt"

    ```text
    kronrod
    ```

### pip

```bash
pip install --upgrade kronrod
```

### uv

```bash
# Add to a uv project
uv add kronrod

# Or install into the current environment
uv pip install kronrod
```

### pipenv

```bash
pipenv install kronrod
```

### poetry

```bash
poetry add kronrod
```

### pdm

```bash
pdm add kronrod
```

### hatch

```bash
hatch add kronrod
```

## [git](https://github.com/eggzec/kronrod)

Install the latest development version directly from the repository:

```bash
pip install --upgrade "git+https://github.com/eggzec/kronrod.git#egg=kronrod"
```

### Building locally

Clone and build from source to modify the C code or test local changes:

```bash
git clone https://github.com/eggzec/kronrod.git
cd kronrod
uv pip install .
```

This invokes the `meson` build system to compile the C sources via `f2py`
and install the resulting extension module.

!!! warning "C compiler required"
    Source builds require a working C compiler. On most Linux distributions,
    install `gcc`:

    === "Debian / Ubuntu"
        ```bash
        sudo apt install gcc
        ```
    === "Fedora"
        ```bash
        sudo dnf install gcc
        ```
    === "macOS (Homebrew)"
        ```bash
        brew install gcc
        ```
    === "Windows"
        Install [MinGW-w64](https://www.mingw-w64.org/) with gcc or use MSYS2.

## Verifying the installation

After installation, verify that the package loads correctly:

```python
import kronrod

x, w1, w2 = kronrod.kronrod(3)
print(x)  # [0.96049127 0.77459667 0.43424375 0.        ]
print(w1)  # [0.10465623 0.26848809 0.40139741 0.45091654]
print("kronrod is working!")
```

## Dependencies

- Python >=3.10
- [numpy](https://pypi.org/project/numpy)
