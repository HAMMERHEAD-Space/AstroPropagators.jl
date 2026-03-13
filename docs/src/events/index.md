# Event System

AstroPropagators.jl provides a comprehensive event detection and handling system built on DifferentialEquations.jl's `ContinuousCallback` infrastructure. The system is designed around two core principles:

1. **Coordinate-agnostic**: All detectors and maneuver builders dispatch on `Type{<:AstroCoord}` (e.g., `Cartesian`, `EDromo`), not on propagator types. This means you can use them with `propagate()` *or* directly with `ODEProblem`.
2. **Decoupled triggers and effects**: Maneuvers separate *when* something happens (triggers) from *what* happens (effects), enabling flexible composition.

## Architecture

The event system is organized into four layers:

```
┌─────────────────────────────────────────────────────┐
│                   Builder API                        │
│  build_maneuver_callback, build_event_callback, ...  │
├─────────────────────────────────────────────────────┤
│              Triggers & Effects                       │
│  TimeTrigger, EventTrigger, FixedDeltaV, ...         │
├─────────────────────────────────────────────────────┤
│            Event Condition Functions                  │
│  apside_condition, eclipse_condition, ...            │
├─────────────────────────────────────────────────────┤
│              State Adapter Layer                      │
│  get_cartesian, get_keplerian, get_physical_time     │
└─────────────────────────────────────────────────────┘
```

### State Adapter Layer

The bottom layer converts raw integrator state vectors to `AstroCoord` types. All higher layers go through this, so coordinate-specific logic is centralized here.

Key functions:
- `get_cartesian(u, t, μ, ::Type{C})` — convert state to `Cartesian`
- `get_keplerian(u, t, μ, ::Type{C})` — convert state to `Keplerian`
- `get_physical_time(u, t, ::Type{C})` — recover physical time (important for regularized sets)
- `coord_type(::AbstractPropagator)` — bridge from propagator types to coord types

### Event Condition Functions

Condition functions return a closure `(u, t, integrator) -> Float64` where zero-crossings trigger events. They are organized into categories:

- **[Orbital Detectors](orbital_detectors.md)**: Apsides, node crossings, anomaly crossings, RAAN
- **[Eclipse Detectors](eclipse_detectors.md)**: Shadow entry/exit using AstroForceModels shadow models
- **[Geometric Detectors](geometric_detectors.md)**: Altitude, latitude, longitude, beta angle
- **[Utility Detectors](utility_detectors.md)**: Date/time crossing, boolean composition, time shifting

### Maneuver System

The maneuver system decouples *when* a maneuver fires from *what* happens:

- **[Maneuvers](maneuvers.md)**: Triggers (`TimeTrigger`, `EventTrigger`), effects (`FixedDeltaV`, `ComputedDeltaV`), actions, and the builder API

## Quick Start

### With `propagate()`

```julia
using AstroPropagators, AstroForceModels, AstroCoords, ComponentArrays
using SatelliteToolboxTransformations

u0_cart = [-1076.225324679696, -6765.896364327722, -332.3087833503755,
            9.356857417032581, -3.3123476319597557, -1.1880157328553503]

JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
grav_model = KeplerianGravityAstroModel()
μ = grav_model.μ
p = ComponentVector(; JD=JD, μ=μ)
models = CentralBodyDynamicsModel(grav_model)
tspan = (0.0, 86400.0)

# Detect all apsides
cond = apside_condition(Cartesian)
log = []
cb = ContinuousCallback(cond,
    (integrator) -> push!(log, (t=integrator.t, type=:periapsis)),
    (integrator) -> push!(log, (t=integrator.t, type=:apoapsis)),
)

sol = propagate(CowellPropagator(), u0_cart, p, models, tspan; callback=cb)
```

### With `ODEProblem` directly

```julia
using OrdinaryDiffEqVerner, SciMLBase

f(u, p, t) = Cowell_EOM(u, p, t, models)
prob = ODEProblem{false}(f, u0_cart, tspan, p)

# Same condition function — no change needed
cond = apside_condition(Cartesian)
cb = ContinuousCallback(cond,
    (integrator) -> push!(log, (t=integrator.t, type=:periapsis)),
    (integrator) -> push!(log, (t=integrator.t, type=:apoapsis)),
)

sol = solve(prob, Vern9(); abstol=1e-13, reltol=1e-13, callback=cb)
```

### With regularized coordinates (EDromo)

Regularized formulations require a `RegularizedCoordinateConfig` as an additional argument:

```julia
config = RegularizedCoordinateConfig(u0_cart, μ; W=W, t₀=0.0, flag_time=PhysicalTime())
ϕ₀ = compute_initial_phi(u0_cart, μ, config)
u0 = Array(EDromo(Cartesian(u0_cart), μ, ϕ₀, config))

# Same detector, just pass the coord type and config
cond = apside_condition(EDromo, config)
end_cb = build_termination_callback(86400.0, EDromo, config)

cb = CallbackSet(
    ContinuousCallback(cond, affect!, affect!),
    end_cb,
)

sol = propagate(EDromoPropagator(), u0, p, models, (ϕ₀, ϕ₀ + 6π), config; callback=cb)
```

## Supported Coordinate Types

All detectors work with any `AstroCoord` type. The `coord_type` bridge maps propagators to their native coordinate:

| Propagator | AstroCoord Type | Regularized | Config Required |
|---|---|---|---|
| `CowellPropagator()` | `Cartesian` | No | No |
| `GaussVEPropagator()` | `Keplerian` | No | No |
| `MilankovichPropagator()` | `Milankovich` | No | No |
| `USM7Propagator()` | `USM7` | No | No |
| `USM6Propagator()` | `USM6` | No | No |
| `USMEMPropagator()` | `USMEM` | No | No |
| `EDromoPropagator()` | `EDromo` | Yes | Yes |
| `KSPropagator()` | `KustaanheimoStiefel` | Yes | Yes |
| `StiSchePropagator()` | `StiefelScheifele` | Yes | Yes |
| `GEqOEPropagator()` | `GEqOE` | Yes | Yes |

## Integrator Parameter Requirements

Most detectors only require `integrator.p.μ`. Some additionally require:

| Parameter | Required By |
|---|---|
| `integrator.p.μ` | All detectors and maneuver builders |
| `integrator.p.JD` | Eclipse, longitude, beta angle detectors |
