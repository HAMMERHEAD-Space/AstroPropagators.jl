# API Reference

This page provides a comprehensive reference for the AstroPropagators.jl API. All propagation methods use the unified `propagate` interface while supporting different coordinate systems and integration options.

## Core Interface

### propagate Function

The primary interface for all propagation methods:

```@docs
propagate
```

## Propagator Types

AstroPropagators defines different propagation methods through type dispatch:

```@docs
AbstractPropType
CowellPropagator
GaussVEPropagator
EDromoPropagator
KSPropagator
MilankovichPropagator
StiefelScheiffePropagator
USM7Propagator
USM6Propagator
USMEMPropagator
```

## Equations of Motion

Each propagator implements its specific equations of motion:

### Cowell Method

```@docs
Cowell_EOM
Cowell_EOM!
```

### Gauss Variational Equations

```@docs
GaussVE_EOM
GaussVE_EOM!
```

### EDromo Method

```@docs
EDromo_EOM
EDromo_EOM!
EDromo_time_condition
end_EDromo_integration
```

### Kustaanheimo-Stiefel Method

```@docs
KS_EOM
KS_EOM!
```

### Unified State Model

```@docs
USM7_EOM
USM7_EOM!
USM6_EOM
USM6_EOM!
USMEM_EOM
USMEM_EOM!
```

### Milankovich Method

```@docs
Milankovich_EOM
Milankovich_EOM!
```

## Configuration

### RegularizedCoordinateConfig

For regularized methods (EDromo, K-S, etc.):

```@docs
RegularizedCoordinateConfig
```

## Coordinate Transformations

### EDromo Transformations

```@docs
cartesian_to_edromo
edromo_to_cartesian
```

### Kustaanheimo-Stiefel Transformations

```@docs
cartesian_to_ks
ks_to_cartesian
```

## Event Handling

### Impulsive Maneuvers

```@docs
EDromo_burn
impulsive_burn_edromo!
```

## Utility Functions

### General Utilities

```@docs
build_dynamics_model
RTN_frame
```

## Type Definitions

## Integration with DifferentialEquations.jl

AstroPropagators leverages the DifferentialEquations.jl ecosystem:

### Solver Selection

```julia
using OrdinaryDiffEqVerner, OrdinaryDiffEqTsit5, OrdinaryDiffEqRadau

# High-precision integrators
sol = propagate(u0, p, models, tspan;
                prop_type = CowellPropagator(),
                ODE_solver = Vern8())  # 8th-order Verner

# Balanced performance
sol = propagate(u0, p, models, tspan;
                prop_type = CowellPropagator(),
                ODE_solver = Tsit5())  # 5th-order Runge-Kutta

# Stiff problems
sol = propagate(u0, p, models, tspan;
                prop_type = CowellPropagator(),
                ODE_solver = Radau5())  # Implicit Radau
```

### Tolerance Settings

```julia
# High precision (recommended for precise orbits)
sol = propagate(u0, p, models, tspan;
                prop_type = CowellPropagator(),
                abstol = 1e-14, reltol = 1e-14)

# Balanced accuracy and speed
sol = propagate(u0, p, models, tspan;
                prop_type = CowellPropagator(),
                abstol = 1e-10, reltol = 1e-10)

# Fast integration (for preliminary analysis)
sol = propagate(u0, p, models, tspan;
                prop_type = CowellPropagator(),
                abstol = 1e-8, reltol = 1e-8)
```

### Callback Integration

```julia
using DifferentialEquations

# Event detection callback
function altitude_condition(u, t, integrator)
    altitude = norm(u[1:3]) - 6378.137  # km
    altitude - 200.0  # Trigger at 200 km altitude
end

cb = ContinuousCallback(altitude_condition, terminate!)

sol = propagate(u0, p, models, tspan;
                prop_type = CowellPropagator(),
                callback = cb)
```

## Examples and Tutorials

For practical examples and step-by-step tutorials, see:

- **[Usage Guide](usage.md)**: Comprehensive usage examples
- **[Propagator Documentation](../propagators/)**: Detailed method descriptions
- **Test suite**: Working examples in the `/test` directory
- **Benchmarks**: Performance comparisons in the `/benchmark` directory

## Support and Development

### Contributing

Contributions are welcome! Please see the development guidelines:

1. **Fork the repository** and create a feature branch
2. **Write tests** for new functionality
3. **Follow coding style** conventions
4. **Document new features** with docstrings
5. **Submit a pull request** with a clear description

### Reporting Issues

Please report bugs and feature requests on the GitHub issue tracker:

- **Provide a minimal working example** that demonstrates the issue
- **Include version information** for Julia and all packages
- **Describe expected vs. actual behavior**
- **Include error messages** and stack traces if applicable

### Getting Help

For questions and support:

- **Check the documentation** first
- **Search existing issues** on GitHub
- **Ask questions** in the Julia Discourse forum with the "astro" tag
- **Review the source code** for implementation details

## Version History

See [CHANGELOG.md](../../../CHANGELOG.md) for detailed version history and breaking changes.