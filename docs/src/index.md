AstroPropagators.jl
================================

AstroPropagators.jl is a comprehensive Julia package for high-performance orbital trajectory propagation using multiple coordinate systems and integration methods. Designed as part of the **SatelliteToolbox.jl** ecosystem, it provides efficient, accurate, and differentiable implementations of classical and modern orbital propagation techniques.

## Features

This package provides robust implementations of the following propagation methods:
- [x] Cowell
- [x] Gauss Variational Equations
- [x] EDromo
- [x] Kustaanheimo-Stiefel
- [x] Milankovich
- [x] Stiefel-Scheifel
- [X] Unified State Model
    - [x] USM7
    - [x] USM6
    - [x] USMEM
- [ ] GEqOE

## Design Philosophy

AstroPropagators.jl is built with the following principles:

### Performance-Oriented
- **Zero-allocation** propagation loops for long-term integration
- **Optimized coordinate transformations** using StaticArrays
- **Efficient force model evaluation** through AstroForceModels.jl integration

### Multiple Coordinate Systems
- **Coordinate system agnostic**: Choose the best method for your problem
- **Automatic conversions** between different element sets
- **Regularized methods** for challenging orbit regimes

### Scientific Accuracy
- **Validated implementations** against reference solutions and literature
- **Numerical stability** through appropriate regularization techniques
- **Conservation properties** maintained where applicable

### Flexible Integration
- **User-defined force models** or pre-built models from AstroForceModels.jl
- **Multiple integrators** from DifferentialEquations.jl ecosystem
- **Event handling** for maneuvers and mission events

## Installation

Install AstroPropagators.jl using Julia's package manager:

```julia
julia> using Pkg
julia> Pkg.add("AstroPropagators")
```

## Applications

AstroPropagators.jl is suitable for a wide range of astrodynamics applications:

### Mission Design and Analysis
- **Trajectory optimization**: Multi-revolution orbit transfers
- **Station-keeping analysis**: Long-term orbit maintenance
- **Constellation design**: Multi-satellite system propagation
- **Collision avoidance**: High-accuracy orbit prediction

### Scientific Research
- **Orbit determination**: Precise trajectory reconstruction
- **Celestial mechanics**: N-body problem investigations  
- **Astrodynamics research**: Testing new propagation methods
- **Space debris studies**: Long-term orbital evolution

### Operational Applications
- **Satellite operations**: Real-time orbit propagation
- **Ground station tracking**: Visibility predictions
- **Launch analysis**: Trajectory validation
- **Formation flying**: Relative motion studies

### Educational Use
- **Orbital mechanics teaching**: Demonstrate different methods
- **Method comparison**: Understand propagation trade-offs
- **Algorithm development**: Test new coordinate systems
- **Numerical analysis**: Study integration techniques

## Integration with SatelliteToolbox.jl

AstroPropagators.jl works seamlessly with the broader SatelliteToolbox.jl ecosystem:

- **[AstroCoords.jl](https://github.com/jmurphy6895/AstroCoords.jl)**: Coordinate system transformations
- **[AstroForceModels.jl](https://github.com/jmurphy6895/AstroForceModels.jl)**: Force model implementations
- **[SatelliteToolbox.jl](https://github.com/JuliaSpace/SatelliteToolbox.jl)**: Astrodynamics functionality

## Documentation Structure

This documentation is organized as follows:

- **[Usage Guide](man/usage.md)**: Comprehensive examples and tutorials
- **[API Reference](man/api.md)**: Complete function and type documentation  
- **[Propagators](propagators/index.md)**: Detailed descriptions of each propagation method
- **[Library Reference](lib/library.md)**: Auto-generated API documentation

*Performance varies with orbit characteristics and force model complexity*

## Contributing

We welcome contributions to AstroPropagators.jl! Whether you're:

- **Reporting bugs** or requesting features
- **Improving documentation** or adding examples
- **Implementing new propagators** or optimizations
- **Adding test cases** or benchmarks

Please see our [contribution guidelines](https://github.com/jmurphy6895/AstroPropagators.jl/blob/master/CONTRIBUTING.md) and feel free to open issues or pull requests on [GitHub](https://github.com/jmurphy6895/AstroPropagators.jl).

## Citation

If you use AstroPropagators.jl in your research, please consider citing:

```bibtex
@software{astropropagators_jl,
  title = {AstroPropagators.jl: High-Performance Orbital Propagation Methods},
  author = {Jordan Murphy},
  url = {https://github.com/jmurphy6895/AstroPropagators.jl},
  year = {2024}
}
```

## License

AstroPropagators.jl is released under the MIT License. See [LICENSE](https://github.com/jmurphy6895/AstroPropagators.jl/blob/master/LICENSE) for details.

## Acknowledgments

This package builds upon decades of research in orbital mechanics and astrodynamics, with particular inspiration from:

- **SatelliteToolbox**: Tons of astrodynamics models and functions to build off of
- **THALASSA library**: Reference implementation for several methods
- **DifferentialEquations.jl**: Robust ODE integration capabilities
- **AstroCoords.jl**: Coordinate system transformation foundation
- **Classical Literature**: Battin, Vallado, Montenbruck, and others
