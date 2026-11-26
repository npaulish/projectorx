# Projector Extensions

![PyPI version](https://img.shields.io/pypi/v/projectorx.svg)

Tools to generate required projectors for Wannierization based on pseudo-atomic orbitals (PAOs). In case a pseudopotential does not contain enough orbitals, PAOs are extended with additional orbitals from OpenMX library: https://www.openmx-square.org.

For more details refer to Y. Jiang et al. "Robust Wannierization including magnetization and spin-orbit coupling via projectability disentanglement", arXiv:2507.06840, https://doi.org/10.48550/arXiv.2507.06840


* PyPI package: https://pypi.org/project/projectorx/
* Free software: MIT License
* Documentation: [docs/](docs/index.md)

## Features

* Extend a UPF pseudopotential with additional projectors: `projectorx extend-upf INPUT_FILE [--output-dir DIR] [--orbital 4p,5s]`.
* Extend every UPF pseudopotential in a directory: `projectorx extend-library INPUT_DIR [--output-dir DIR]`.
* Missing orbitals are built as hydrogenic functions, either orthogonal to the inner shell with the same l or fitted to an OpenMX pseudo-atomic orbital (downloaded on first use).
* Spin-orbit pseudopotentials are supported.

See [Usage](docs/usage.md) for the options, how orbitals are chosen, and the output format.


## Credits

This package was created with [Cookiecutter](https://github.com/audreyfeldroy/cookiecutter) and the [audreyfeldroy/cookiecutter-pypackage](https://github.com/audreyfeldroy/cookiecutter-pypackage) project template.
