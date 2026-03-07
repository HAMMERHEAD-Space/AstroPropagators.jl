# Stiefel-Scheifele Method

The Stiefel-Scheifele method is a canonical regularization technique that transforms the orbital motion problem using universal variables and canonical transformations. This method provides excellent numerical stability for eccentric orbits while preserving the Hamiltonian structure of the dynamics.

## Physical Description

The Stiefel-Scheifele method employs:

- **Universal variables**: Elements that work for all orbit types (elliptic, parabolic, hyperbolic)
- **Canonical transformation**: Preserves the symplectic structure of Hamiltonian mechanics
- **Regularized time**: Uses a transformed time coordinate for uniform integration steps
- **Energy preservation**: Maintains conservation properties better than non-canonical methods



## Usage Examples

!!! note "Fictitious time"
    Stiefel-Scheifele integrates in a fictitious time variable `ϕ`, not physical
    time `t`. The `tspan` argument is therefore the range of `ϕ`, and a
    `ContinuousCallback` (via `end_StiSche_integration`) terminates the integration
    when the desired physical time is reached.

### Keplerian Propagation

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

W = (
    potential(Cartesian(u0_cart), p, 0.0, grav_model) -
    potential(Cartesian(u0_cart), p, 0.0, KeplerianGravityAstroModel(; μ=μ))
)
config = RegularizedCoordinateConfig(u0_cart, μ; W=W, t₀=0.0, flag_time=PhysicalTime())

ϕ₀ = compute_initial_phi(u0_cart, μ, config)
u0 = Array(StiefelScheifele(Cartesian(u0_cart), μ, ϕ₀, config))
tspan = (ϕ₀, ϕ₀ + 6π)

sol = propagate(
    StiSchePropagator(), u0, p, models, tspan, config;
    callback=end_StiSche_integration(86400.0, config),
)
```

### High-Fidelity Propagation

```julia
using AstroPropagators, AstroForceModels, AstroCoords, ComponentArrays
using SatelliteToolboxGravityModels, SatelliteToolboxTransformations, SpaceIndices

SpaceIndices.init()
eop_data = fetch_iers_eop()
grav_coeffs = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96))

u0_cart = [-1076.225324679696, -6765.896364327722, -332.3087833503755,
            9.356857417032581, -3.3123476319597557, -1.1880157328553503]

JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
grav_model = GravityHarmonicsAstroModel(;
    gravity_model=grav_coeffs, eop_data=eop_data, order=36, degree=36,
)
μ = GravityModels.gravity_constant(grav_model.gravity_model) / 1E9
p = ComponentVector(; JD=JD, μ=μ)

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

W = (
    potential(Cartesian(u0_cart), p, 0.0, grav_model) -
    potential(Cartesian(u0_cart), p, 0.0, KeplerianGravityAstroModel(; μ=μ))
)
config = RegularizedCoordinateConfig(u0_cart, μ; W=W, t₀=0.0, flag_time=PhysicalTime())

ϕ₀ = compute_initial_phi(u0_cart, μ, config)
u0 = Array(StiefelScheifele(Cartesian(u0_cart), μ, ϕ₀, config))
tspan = (ϕ₀, ϕ₀ + 6π)

sol = propagate(
    StiSchePropagator(), u0, p, models, tspan, config;
    callback=end_StiSche_integration(86400.0, config),
)
```

### Using `ODEProblem` Directly

```julia
using OrdinaryDiffEqVerner, SciMLBase

# Using the same force model and config setup from above...
ϕ₀ = compute_initial_phi(u0_cart, μ, config)
u0 = Array(StiefelScheifele(Cartesian(u0_cart), μ, ϕ₀, config))
tspan = (ϕ₀, ϕ₀ + 6π)

# In-place (what propagate! uses internally)
f!(du, u, p, t) = StiSche_EOM!(du, u, p, t, models, config)
prob = ODEProblem(f!, u0, tspan, p)
sol = solve(
    prob, Vern9();
    abstol=1e-13, reltol=1e-13,
    callback=end_StiSche_integration(86400.0, config),
)

# Out-of-place (what propagate uses internally)
f(u, p, t) = StiSche_EOM(u, p, t, models, config)
prob = ODEProblem{false}(f, u0, tspan, p)
sol = solve(
    prob, Vern9();
    abstol=1e-13, reltol=1e-13,
    callback=end_StiSche_integration(86400.0, config),
)
```
## Applications

The Stiefel-Scheifele method is ideal for:

### High-Precision Applications
- **Long-term integration**: Multi-year propagation with conservation
- **Precision orbit determination**: Scientific mission requirements
- **Fundamental research**: Studies requiring mathematical rigor

### Variable Orbit Types
- **Interplanetary missions**: Hyperbolic departure and arrival
- **Comet studies**: Highly eccentric natural orbits
- **Transfer trajectories**: Orbits that change type during mission

### Theoretical Studies
- **Hamiltonian mechanics**: Canonical formulation benefits
- **Conservation analysis**: Energy and momentum preservation
- **Method validation**: Gold standard for regularized methods

## References

[1]: Stiefel, E. L., and G. Scheifele. "Linear and Regular Celestial Mechanics." Springer-Verlag, 1971.
[2]: Scheifele, G. "On numerical integration of perturbed linear oscillating systems." Zeitschrift für angewandte Mathematik und Physik, 22(2), pp. 186-210, 1971.
[3]: Battin, Richard H. "An Introduction to the Mathematics and Methods of Astrodynamics." AIAA, 1999.
[4]: Bond, V. R., and M. C. Allman. "Modern Astrodynamics: Fundamentals and Perturbation Methods." Princeton University Press, 1996. 