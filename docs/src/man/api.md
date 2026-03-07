# API Reference

This page provides a comprehensive reference for the AstroPropagators.jl API. All propagation methods use the unified `propagate` / `propagate!` interface while supporting different coordinate systems and integration options.

## Type Hierarchy

```
AbstractPropagator
├── AbstractStandardPropagator        # No config needed
│   ├── CowellPropagator             # Cartesian EOM
│   ├── GaussVEPropagator            # Gauss VE (Keplerian elements)
│   ├── MilankovichPropagator        # Milankovich elements
│   ├── USM7Propagator               # Unified State Model (7-element)
│   ├── USM6Propagator               # Unified State Model (6-element, MRP)
│   └── USMEMPropagator              # Unified State Model (exponential map)
└── AbstractRegularizedPropagator     # Requires RegularizedCoordinateConfig
    ├── EDromoPropagator             # EDromo elements
    ├── KSPropagator                 # Kustaanheimo-Stiefel
    ├── StiSchePropagator            # Stiefel-Scheifele
    └── GEqOEPropagator              # Generalized Equinoctial Elements
```

## `propagate` / `propagate!` Interface

Two variants are available following Julia naming conventions:

- **`propagate`** — builds an out-of-place `ODEProblem{false}` using `eom`. Best for `StaticArrays` and small state vectors.
- **`propagate!`** — builds an in-place `ODEProblem{true}` using `eom!`. Preferred for mutable arrays and larger systems.

Both accept the same arguments and return an `ODESolution`.

### Standard Propagators

```julia
sol = propagate(CowellPropagator(), u0, p, models, tspan)
sol = propagate!(CowellPropagator(), u0, p, models, tspan)
sol = propagate(USM7Propagator(), u0, p, models, tspan; solver=Vern8())
```

### Regularized Propagators

Regularized propagators require a `RegularizedCoordinateConfig` as a positional argument:

```julia
sol = propagate(EDromoPropagator(), u0, p, models, tspan, config)
sol = propagate!(EDromoPropagator(), u0, p, models, tspan, config)
sol = propagate(GEqOEPropagator(), u0, p, models, tspan, config)
```

## `eom` / `eom!` Dispatch

The EOM dispatch layer provides a unified interface to all equations of motion:

- **`eom(prop, u, p, t, models[, config])`** — out-of-place, returns derivative vector
- **`eom!(prop, du, u, p, t, models[, config])`** — in-place, writes into `du`

```julia
# Out-of-place
du = eom(CowellPropagator(), u, p, t, models)
du = eom(EDromoPropagator(), u, p, ϕ, models, config)

# In-place
eom!(CowellPropagator(), du, u, p, t, models)
eom!(EDromoPropagator(), du, u, p, ϕ, models, config)
```

## Integration with DifferentialEquations.jl

AstroPropagators leverages the DifferentialEquations.jl ecosystem:

### Solver Selection

```julia
using OrdinaryDiffEqVerner, OrdinaryDiffEqTsit5, OrdinaryDiffEqRadau

# High-precision integrators
sol = propagate(CowellPropagator(), u0, p, models, tspan;
                solver = Vern8())  # 8th-order Verner

# Balanced performance
sol = propagate(CowellPropagator(), u0, p, models, tspan;
                solver = Tsit5())  # 5th-order Runge-Kutta

# Stiff problems
sol = propagate(CowellPropagator(), u0, p, models, tspan;
                solver = Radau5())  # Implicit Radau
```

### Tolerance Settings

```julia
# High precision (recommended for precise orbits)
sol = propagate(CowellPropagator(), u0, p, models, tspan;
                abstol = 1e-14, reltol = 1e-14)

# Balanced accuracy and speed
sol = propagate(CowellPropagator(), u0, p, models, tspan;
                abstol = 1e-10, reltol = 1e-10)

# Fast integration (for preliminary analysis)
sol = propagate(CowellPropagator(), u0, p, models, tspan;
                abstol = 1e-8, reltol = 1e-8)
```

### Callback Integration

Additional keyword arguments are forwarded directly to `OrdinaryDiffEq.solve`:

```julia
using DifferentialEquations

function altitude_condition(u, t, integrator)
    altitude = norm(u[1:3]) - 6378.137  # km
    altitude - 200.0  # Trigger at 200 km altitude
end

cb = ContinuousCallback(altitude_condition, terminate!)

sol = propagate(CowellPropagator(), u0, p, models, tspan;
                callback = cb)
```

## Extending with New Propagators

Define a new propagator by subtyping and implementing both `eom` and `eom!`:

```julia
struct MyPropagator <: AbstractStandardPropagator end

@inline function AstroPropagators.eom(::MyPropagator, u, p, t, models)
    # compute and return derivative vector (out-of-place)
end

@inline function AstroPropagators.eom!(::MyPropagator, du, u, p, t, models)
    # compute derivatives and write into du (in-place)
end
```

For regularized propagators:

```julia
struct MyRegularizedPropagator <: AbstractRegularizedPropagator end

@inline function AstroPropagators.eom(::MyRegularizedPropagator, u, p, t, models, config)
    # compute and return derivative vector (out-of-place)
end

@inline function AstroPropagators.eom!(
    ::MyRegularizedPropagator, du, u, p, t, models, config
)
    # compute derivatives and write into du (in-place)
end
```

## Examples and Tutorials

For practical examples and step-by-step tutorials, see:

- **[Usage Guide](usage.md)**: Comprehensive usage examples
- **[Propagator Documentation](../propagators/index.md)**: Detailed method descriptions
- **Test suite**: Working examples in the `/test` directory
- **Benchmarks**: Performance comparisons in the `/benchmark` directory
