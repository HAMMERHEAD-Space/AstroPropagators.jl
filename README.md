# AstroPropagators

## IMPORTANT NOTE: While the existing functions will work, for the most performance and differentiation capability, use one of the PR and instructions detailed there.

[![CI](https://github.com/jmurphy6895/AstroPropagators.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/jmurphy6895/AstroPropagators.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![codecov](https://codecov.io/gh/jmurphy6895/AstroPropagators.jl/branch/main/graph/badge.svg?token=47G4OLV6PD)](https://codecov.io/gh/jmurphy6895/AstroPropagators.jl)
[![](https://img.shields.io/badge/docs-stable-blue.svg)][docs-stable-url]
[![](https://img.shields.io/badge/docs-dev-blue.svg)][docs-dev-url]
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

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

[docs-dev-url]: https://jmurphy6895.github.io/AstroPropagators.jl/dev/
[docs-stable-url]: https://jmurphy6895.github.io/AstroPropagators.jl/dev/
