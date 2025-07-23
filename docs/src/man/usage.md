# Usage

AstroPropagators.jl provides a comprehensive framework for orbital trajectory propagation using multiple coordinate systems and integration methods. The package integrates with AstroForceModels.jl and AstroCoords.jl ecosystems while providing access to differential equation solvers through OrdinaryDiffEq.jl.

## Basic Concepts

### Equations of Motion (EOM) Functions

Each propagator implements its specific equations of motion that can be used with ODE solvers:

- `Cowell_EOM!(du, u, p, t, models)`: Cartesian coordinates [km, km/s]
- `GaussVE_EOM!(du, u, p, t, models)`: Keplerian elements [km, rad, rad, rad, rad, rad]
- `EDromo_EOM!(du, u, p, t, models, config)`: Regularized elements
- `KS_EOM!(du, u, p, t, models, config)`: Kustaanheimo-Stiefel coordinates
- `USM7_EOM!(du, u, p, t, models)`: Unified State Model elements

### DynamicsModel Structure

Force models are wrapped in a `DynamicsModel` struct

## Complete Usage Example

Here's a simple example showing different propagation methods:

```julia
using AstroPropagators
using AstroForceModels, AstroCoords
using ComponentArrays
using OrdinaryDiffEq

# Initial conditions (ISS-like orbit in km, km/s)
r0 = [6378.137 + 408, 0.0, 0.0]  # Position [km]
v0 = [0.0, 7.660, 0.0]           # Velocity [km/s]  
u0_cart = [r0; v0]

# Parameters
μ = 3.986004415e5  # km³/s²
JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
p = ComponentVector(JD=JD, μ=μ)
tspan = (0.0, 86400.0)  # 24 hours

# Force models
gravity_model = KeplerianGravityAstroModel(μ=μ)
models = DynamicsModel(gravity_model)

# 1. Cowell Method (Cartesian)
function cowell_dynamics!(du, u, p, t)
    Cowell_EOM!(du, u, p, t, models)
end

prob_cowell = ODEProblem(cowell_dynamics!, u0_cart, tspan, p)
sol_cowell = solve(prob_cowell, Vern8(), abstol=1e-12, reltol=1e-12)

# 2. Gauss VE Method (Keplerian)
# Convert to Keplerian elements
cart_state = Cartesian(u0_cart, μ=μ)
kep_state = Keplerian(cart_state)
u0_kep = [kep_state.a, kep_state.e, kep_state.i, kep_state.Ω, kep_state.ω, kep_state.f]

function gaussve_dynamics!(du, u, p, t)
    GaussVE_EOM!(du, u, p, t, models)
end

prob_gauss = ODEProblem(gaussve_dynamics!, u0_kep, tspan, p)
sol_gauss = solve(prob_gauss, Vern8(), abstol=1e-12, reltol=1e-12)

# 3. EDromo Method (Regularized)
# Create configuration and convert to EDromo elements  
config = RegularizedCoordinateConfig(u0_cart, μ; flag_time=PhysicalTime())
u0_edromo = cart2EDromo(u0_cart, μ, 0.0, config)

function edromo_dynamics!(du, u, p, t)
    EDromo_EOM!(du, u, p, t, models, config)
end

prob_edromo = ODEProblem(edromo_dynamics!, u0_edromo, tspan, p)
sol_edromo = solve(prob_edromo, Vern8(), abstol=1e-12, reltol=1e-12)

println("All propagation methods completed successfully!")
println("Cowell final position: $(sol_cowell.u[end][1:3]) km")
println("Gauss VE final elements: $(sol_gauss.u[end])")
println("EDromo final state: $(sol_edromo.u[end])")
```

## Troubleshooting

1. **Unit consistency**: Ensure all quantities use km-based units
2. **Coordinate matching**: Match initial conditions to propagator requirements  
3. **Force model units**: Verify force models use correct unit system
4. **RegularizedCoordinateConfig**: Use `RegularizedCoordinateConfig(state, μ)` constructor
5. **Integration tolerances**: Balance accuracy with computational cost