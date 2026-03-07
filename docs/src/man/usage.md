# Usage

AstroPropagators.jl provides a unified `propagate` interface for orbital trajectory propagation across multiple coordinate formulations, integrated with AstroForceModels.jl for force modeling and AstroCoords.jl for coordinate transformations.

## Quick Start

### Keplerian Propagation (Two-Body)

```julia
using AstroPropagators
using AstroForceModels, AstroCoords
using ComponentArrays
using SatelliteToolboxTransformations

# Initial Cartesian state [km, km/s]
u0_cart = [-1076.225324679696, -6765.896364327722, -332.3087833503755,
            9.356857417032581, -3.3123476319597557, -1.1880157328553503]

# Parameters
JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
grav_model = KeplerianGravityAstroModel()
μ = grav_model.μ
p = ComponentVector(; JD=JD, μ=μ)

# Force model (two-body only)
models = CentralBodyDynamicsModel(grav_model)
tspan = (0.0, 86400.0)  # 1 day

# Propagate with Cowell (Cartesian)
sol = propagate(CowellPropagator(), u0_cart, p, models, tspan)
```

### High-Fidelity Propagation

```julia
using AstroPropagators
using AstroForceModels, AstroCoords
using ComponentArrays
using SatelliteToolboxGravityModels, SatelliteToolboxTransformations
using SpaceIndices

SpaceIndices.init()
eop_data = fetch_iers_eop()
grav_coeffs = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96))

# Initial Cartesian state [km, km/s]
u0_cart = [-1076.225324679696, -6765.896364327722, -332.3087833503755,
            9.356857417032581, -3.3123476319597557, -1.1880157328553503]

JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)

# Gravity harmonics (36×36)
grav_model = GravityHarmonicsAstroModel(;
    gravity_model=grav_coeffs, eop_data=eop_data, order=36, degree=36,
)
μ = GravityModels.gravity_constant(grav_model.gravity_model) / 1E9
p = ComponentVector(; JD=JD, μ=μ)

# Third-body, SRP, and drag perturbations
sun_model = ThirdBodyModel(; body=SunBody(), eop_data=eop_data)
moon_model = ThirdBodyModel(; body=MoonBody(), eop_data=eop_data)
srp_model = SRPAstroModel(;
    satellite_srp_model=CannonballFixedSRP(0.2),
    sun_data=sun_model, eop_data=eop_data, shadow_model=Conical(),
)
drag_model = DragAstroModel(;
    satellite_drag_model=CannonballFixedDrag(0.2),
    atmosphere_model=JB2008(), eop_data=eop_data,
)

models = CentralBodyDynamicsModel(
    grav_model, (sun_model, moon_model, srp_model, drag_model),
)
tspan = (0.0, 86400.0)

# Propagate with Cowell
sol = propagate(CowellPropagator(), u0_cart, p, models, tspan)
```

## The `propagate` API

The propagator type is always the first argument. Standard and regularized propagators have different signatures:

```julia
# Standard propagators — no config needed
sol = propagate(CowellPropagator(), u0, p, models, tspan)

# Regularized propagators — config is a required positional argument
sol = propagate(EDromoPropagator(), u0, p, models, tspan, config)
```

### Keyword Arguments

All keyword arguments are forwarded to `OrdinaryDiffEq.solve`:

```julia
sol = propagate(CowellPropagator(), u0, p, models, tspan;
    solver=Vern9(),    # ODE solver (default: Vern9)
    abstol=1e-13,      # absolute tolerance
    reltol=1e-13,      # relative tolerance
    callback=cb,       # DiffEq callbacks
    saveat=0.0:60.0:86400.0,  # save every 60 seconds
)
```

## Choosing a Propagator

### Standard Propagators

| Propagator | Coordinates | Best For |
|---|---|---|
| `CowellPropagator()` | Cartesian `[x,y,z,vx,vy,vz]` | General purpose, reference solutions |
| `GaussVEPropagator()` | Keplerian `[a,e,i,Ω,ω,f]` | Small perturbations, element insight |
| `MilankovichPropagator()` | Milankovich `[hx,hy,hz,ex,ey,ez,L]` | Non-singular, perturbation analysis |
| `USM7Propagator()` | USM7 (quaternion) `[C,Rf1,Rf2,ϵ1,ϵ2,ϵ3,η]` | Non-singular, all orbit types |
| `USM6Propagator()` | USM6 (MRP) `[C,Rf1,Rf2,σ1,σ2,σ3]` | Minimal state, non-singular |
| `USMEMPropagator()` | USMEM (exp. map) `[C,Rf1,Rf2,a1,a2,a3]` | Minimal state, non-singular |

### Regularized Propagators

| Propagator | Coordinates | Best For |
|---|---|---|
| `EDromoPropagator()` | EDromo elements | High eccentricity, long-term propagation |
| `KSPropagator()` | Kustaanheimo-Stiefel | Collision orbits, extreme eccentricity |
| `StiSchePropagator()` | Stiefel-Scheifele | Hamiltonian preservation, all orbit types |
| `GEqOEPropagator()` | Generalized Equinoctial | Conservative perturbations, high accuracy |

## Coordinate Conversion

Use AstroCoords.jl to convert between the initial Cartesian state and each propagator's coordinate system:

```julia
# Standard propagators
u0_kep  = Array(Keplerian(Cartesian(u0_cart), μ))
u0_mil  = Array(Milankovich(Cartesian(u0_cart), μ))
u0_usm7 = Array(USM7(Cartesian(u0_cart), μ))

# Regularized propagators need a config
config = RegularizedCoordinateConfig(u0_cart, μ; W=W, t₀=0.0, flag_time=PhysicalTime())
ϕ = AstroCoords.compute_initial_phi(u0_cart, μ, config)

u0_edromo  = Array(EDromo(Cartesian(u0_cart), μ, ϕ, config))
u0_ks      = Array(KustaanheimoStiefel(Cartesian(u0_cart), μ, config))
u0_stische = Array(StiefelScheifele(Cartesian(u0_cart), μ, ϕ, config))

# GEqOE uses a simpler config (only W needed)
config_geqoe = RegularizedCoordinateConfig(; W=W)
u0_geqoe = Array(GEqOE(Cartesian(u0_cart), μ, config_geqoe))
```

## Converting Back to Cartesian

After propagation, convert the final state back to Cartesian:

```julia
# Standard propagators
cart_final = Array(Cartesian(Keplerian(sol.u[end]), μ))

# Regularized propagators
cart_final = Array(Cartesian(GEqOE(sol.u[end]), μ, config))
```

## Using `ODEProblem` Directly

The `propagate` function is a thin wrapper around `ODEProblem` + `solve`. For full
control — custom callbacks, event handling, or composing the EOM closure into a larger
system — you can build the problem yourself.

### Standard Propagators

```julia
using OrdinaryDiffEqVerner, SciMLBase

f!(du, u, p, t) = Cowell_EOM!(du, u, p, t, models)

prob = ODEProblem(f!, u0_cart, tspan, p)
sol = solve(prob, Vern9(); abstol=1e-13, reltol=1e-13)
```

The `eom!` dispatch layer works identically:

```julia
f!(du, u, p, t) = eom!(CowellPropagator(), du, u, p, t, models)
```

### Regularized Propagators

Regularized EOMs require the `config` argument in the closure. EDromo, K-S, and
Stiefel-Scheifele also integrate in fictitious time and need an end-of-integration
callback:

```julia
config = RegularizedCoordinateConfig(u0_cart, μ; W=W, t₀=0.0, flag_time=PhysicalTime())

ϕ₀ = compute_initial_phi(u0_cart, μ, config)
u0 = Array(EDromo(Cartesian(u0_cart), μ, ϕ₀, config))

f!(du, u, p, t) = EDromo_EOM!(du, u, p, t, models, config)

prob = ODEProblem(f!, u0, (ϕ₀, ϕ₀ + 6π), p)
sol = solve(prob, Vern9(); abstol=1e-13, reltol=1e-13,
    callback=end_EDromo_integration(86400.0, config),
)
```

GEqOE integrates in physical time, so no callback is needed:

```julia
config_geqoe = RegularizedCoordinateConfig(; W=W)
u0 = Array(GEqOE(Cartesian(u0_cart), μ, config_geqoe))

f!(du, u, p, t) = GEqOE_EOM!(du, u, p, t, models, config_geqoe)

prob = ODEProblem(f!, u0, tspan, p)
sol = solve(prob, Vern9(); abstol=1e-13, reltol=1e-13)
```

## Extending with New Propagators

Define a new propagator by subtyping and implementing `eom!`:

```julia
# Standard propagator
struct MyPropagator <: AbstractStandardPropagator end

@inline function AstroPropagators.eom!(::MyPropagator, du, u, p, t, models)
    # compute derivatives, write into du
end

# Regularized propagator
struct MyRegPropagator <: AbstractRegularizedPropagator end

@inline function AstroPropagators.eom!(::MyRegPropagator, du, u, p, t, models, config)
    # compute derivatives, write into du
end
```

## Troubleshooting

1. **Unit consistency**: All quantities must use km, km/s, km³/s² units
2. **Coordinate matching**: The initial state must be in the propagator's coordinate system
3. **Config requirement**: Regularized propagators require a `RegularizedCoordinateConfig` — the compiler will error if you forget it
4. **SpaceIndices**: Call `SpaceIndices.init()` before using atmospheric drag models
5. **EOP data**: Earth orientation parameters are required for gravity harmonics, SRP, drag, and third-body models
