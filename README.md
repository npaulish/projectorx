# Projector Extensions

![PyPI version](https://img.shields.io/pypi/v/projectorx.svg)

Tools to generate required projectors for Wannierization based on pseudo-atomic orbitals (PAOs). In case a pseudopotential does not contain enough orbitals, PAOs are extended with additional orbitals from OpenMX library: https://www.openmx-square.org.

For more details the paper (see citation below).

* PyPI package: https://pypi.org/project/projectorx/
* Free software: MIT License
* Documentation: [docs/](docs/index.md)

## Features

* Extend a UPF pseudopotential with additional projectors: `projectorx extend-upf INPUT_FILE [--output-dir DIR] [--orbital 4p,5s]`.
* Extend every UPF pseudopotential in a directory: `projectorx extend-library INPUT_DIR [--output-dir DIR]`.
* Missing orbitals are built as hydrogenic functions, either orthogonal to the inner shell with the same l or fitted to an OpenMX pseudo-atomic orbital (downloaded on first use).
* Spin-orbit pseudopotentials are supported.

See [Usage](docs/usage.md) for the options, how orbitals are chosen, and the output format.


## Acknowledgements

If you use this code in your work, we would be grateful if you could cite:

[Y. Jiang et al., *Robust Wannierization including magnetization and spin-orbit coupling via projectability disentanglement*, **npj Comput. Mater. 11**, 353 (2025)](https://doi.org/10.1038/s41524-025-01835-8)



