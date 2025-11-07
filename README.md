# Projector Extensions

![PyPI version](https://img.shields.io/pypi/v/projectorx.svg)
[![Documentation Status](https://readthedocs.org/projects/projectorx/badge/?version=latest)](https://projectorx.readthedocs.io/en/latest/?version=latest)

Tools to generate required projectors for Wannierization based on pseudo-atomic orbitals (PAOs). In case a pseudopotential does not contain enough orbitals, PAOs are extended with additional orbitals from OpenMX library: https://www.openmx-square.org.

For more details refer to Y. Jiang et al. "Robust Wannierization including magnetization and spin-orbit coupling via projectability disentanglement", arXiv:2507.06840, https://doi.org/10.48550/arXiv.2507.06840


* PyPI package: https://pypi.org/project/projectorx/
* Free software: MIT License
* Documentation: https://projectorx.readthedocs.io.

## Features

* Extend a UPF pseudopotential with additional projectors via `projectorx extend-upf INPUT_FILE OUTPUT_DIR`.
* Extend an existing AiiDA pseudo family with additional atomic projectors via `projectorx extend-family PSEUDO_FAMILY_LABEL`.


## Credits

This package was created with [Cookiecutter](https://github.com/audreyfeldroy/cookiecutter) and the [audreyfeldroy/cookiecutter-pypackage](https://github.com/audreyfeldroy/cookiecutter-pypackage) project template.
