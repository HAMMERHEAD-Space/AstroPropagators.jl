# Kustaanheimo-Stiefel Method

The Kustaanheimo-Stiefel (K-S) method is a regularization technique that transforms the 3-dimensional two-body problem into a 4-dimensional harmonic oscillator. This method completely regularizes the equations of motion, making it ideal for numerical stability and long-term propagation.

## Physical Description

The K-S transformation uses four variables `[u₁, u₂, u₃, u₄]` to represent the three-dimensional position through the bilinear transformation:

```
x = u₁² - u₂² - u₃² + u₄²
y = 2(u₁u₂ - u₃u₄)  
z = 2(u₁u₃ + u₂u₄)
```

The method transforms the independent variable from physical time `t` to a fictitious time `s`, providing complete regularization.

## Key Features

- **Complete Regularization**: Handles collision orbits and arbitrary eccentricity
- **4D Harmonic Oscillator**: Transforms Kepler problem to oscillator equation
- **Quaternion Representation**: Uses spinor-like 4D variables
- **Linear Time Element**: Provides uniform integration properties

## Usage Examples

!!! note "Fictitious time"
    K-S integrates in a fictitious time variable `s`, not physical time `t`.
    The `tspan` argument is therefore the range of `s`, and a `ContinuousCallback`
    (via `end_KS_integration`) terminates the integration when the desired
    physical time is reached.

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

u0 = Array(KustaanheimoStiefel(Cartesian(u0_cart), μ, config))
tspan = (0.0, 9π)

sol = propagate(
    KSPropagator(), u0, p, models, tspan, config;
    callback=end_KS_integration(86400.0, config),
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

u0 = Array(KustaanheimoStiefel(Cartesian(u0_cart), μ, config))
tspan = (0.0, 9π)

sol = propagate(
    KSPropagator(), u0, p, models, tspan, config;
    callback=end_KS_integration(86400.0, config),
)
```

### Using `ODEProblem` Directly

```julia
using OrdinaryDiffEqVerner, SciMLBase

# Using the same force model and config setup from above...
u0 = Array(KustaanheimoStiefel(Cartesian(u0_cart), μ, config))
tspan = (0.0, 9π)

f!(du, u, p, t) = KS_EOM!(du, u, p, t, models, config)

prob = ODEProblem(f!, u0, tspan, p)
sol = solve(
    prob, Vern9();
    abstol=1e-13, reltol=1e-13,
    callback=end_KS_integration(86400.0, config),
)
```

## Applications

The K-S method excels for:

- **Collision orbits**: Trajectories passing through the central body
- **Extreme eccentricity**: Orbits with e ≥ 1 (parabolic/hyperbolic)
- **Close encounters**: Very low perigee passages
- **Numerical stability**: When other methods fail due to singularities

## References

[1]: Kustaanheimo, P., and E. Stiefel. "Perturbation theory of Kepler motion based on spinor regularization." Journal für die reine und angewandte Mathematik, 218, pp. 204-219, 1965.
[2]: Stiefel, E. L., and G. Scheifele. "Linear and Regular Celestial Mechanics." Springer-Verlag, 1971.
[3]: Battin, Richard H. "An Introduction to the Mathematics and Methods of Astrodynamics." AIAA, 1999. 