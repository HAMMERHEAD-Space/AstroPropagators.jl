# AstroPropagators

[![CI](https://github.com/HAMMERHEAD-Space/AstroPropagators.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/HAMMERHEAD-Space/AstroPropagators.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![codecov](https://codecov.io/gh/HAMMERHEAD-Space/AstroPropagators.jl/branch/main/graph/badge.svg?token=47G4OLV6PD)](https://codecov.io/gh/HAMMERHEAD-Space/AstroPropagators.jl)
[![](https://img.shields.io/badge/docs-dev-blue.svg)][docs-dev-url]
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![DOI](https://zenodo.org/badge/672594990.svg)](https://doi.org/10.5281/zenodo.16954455)

## Description

This project implements several propagation methods for the **AstroPropagators.jl** ecosystem using equations of motion for different orbital element sets. The force models can either be user-supplied or built up with the **AstroForceModels.jl** package. This package was inspired by the [THALASSA library](https://github.com/woodywu-arizona/thalassa). A full list of implemented propagators can be found below:

- [x] Cowell
- [x] Gauss Variational Equations
- [ ] EDromo
- [ ] Kustaanheimo-Stiefel
- [x] Milankovich
- [ ] Stiefel-Scheifel
- [X] Unified State Model
    - [x] USM7
    - [x] USM6
    - [x] USMEM


NOTE: Although the code for some of these are present, only the checked ones currently supported!

## Installation

```julia
julia> using Pkg
julia> Pkg.add("AstroPropagators")
```

## Documentation

For more information, see the [documentation][docs-stable-url].

## Citing

If you use `AstroPropagators.jl` in your work, please consider citing it and `AstroForceModels.jl`.

```bibtex
@software{jordan_murphy_2025_16954386,
  author       = {Jordan Murphy},
  title        = {HAMMERHEAD-Space/AstroForceModels.jl: v0.3.13},
  month        = aug,
  year         = 2025,
  publisher    = {Zenodo},
  version      = {v0.3.13},
  doi          = {10.5281/zenodo.16954386},
  url          = {https://doi.org/10.5281/zenodo.16954386},
}
```

```bibtex
@software{jordan_murphy_2025_16954456,
  author       = {Jordan Murphy},
  title        = {HAMMERHEAD-Space/AstroPropagators.jl: v0.1.4},
  month        = aug,
  year         = 2025,
  publisher    = {Zenodo},
  version      = {v0.1.4},
  doi          = {10.5281/zenodo.16954456},
  url          = {https://doi.org/10.5281/zenodo.16954456},
}
```

[docs-dev-url]: https://HAMMERHEAD-Space.github.io/AstroPropagators.jl/dev/
[docs-stable-url]: https://HAMMERHEAD-Space.github.io/AstroPropagators.jl/dev/
