# EDromo Method

The EDromo (Element-based Dromo) method is a modern regularized orbital propagation technique that uses a set of non-singular orbital elements and a fictitious time coordinate. Developed by Baù, Bombardelli, Peláez, and Lorenzini, this method is particularly effective for numerical stability and long-term propagation.

## Physical Description

EDromo uses eight elements `[ζ₁, ζ₂, ζ₃, ζ₄, ζ₅, ζ₆, ζ₇, ζ₈]` that avoid the singularities of classical Keplerian elements:

- **ζ₁, ζ₂**: Eccentricity components related to `e*cos(ω+Ω)` and `e*sin(ω+Ω)`
- **ζ₃**: Reciprocal of semi-latus rectum (1/p)
- **ζ₄, ζ₅**: Inclination components related to `tan(i/2)*cos(Ω)` and `tan(i/2)*sin(Ω)`
- **ζ₆**: Angular variable related to argument of latitude
- **ζ₇**: Time element for time regularization
- **ζ₈**: Time element (depends on `flag_time` setting)

The time element (ζ₈) in EDromo depends on the `flag_time` parameter in the `RegularizedCoordinateConfig`:

- **`PhysicalTime()`**: The time element directly represents the physical time
  ```julia
  ζ₈ = t₀ / TU
  ```

- **`ConstantTime()`**: The time element includes a constant offset plus orbital motion
  ```julia
  ζ₈ = t₀ / TU + ζ₃^(1.5) * (ζ₁*sin(ϕ) - ζ₂*cos(ϕ) - ϕ)
  ```

- **`LinearTime()`**: The time element includes a linear term plus orbital motion
  ```julia
  ζ₈ = t₀ / TU + ζ₃^(1.5) * (ζ₁*sin(ϕ) - ζ₂*cos(ϕ))
  ```

The choice of time formulation affects numerical stability and computational efficiency for different orbit types.

## Key Features

- **Non-singular elements**: Works for all orbit types including circular and equatorial orbits
- **Fictitious time regularization**: Eliminates numerical difficulties near collision
- **High accuracy**: Maintains precision for long-term propagation
- **Multiple time formulations**: Configurable time element for different applications

## Advantages and Disadvantages

### Advantages
- **Regularization**: Excellent performance for high-eccentricity and collision orbits
- **Numerical stability**: Uniform step sizes in fictitious time
- **Conservation properties**: Better long-term conservation of energy and angular momentum
- **Versatility**: Handles all orbit types with single formulation

### Disadvantages
- **Complexity**: More complex mathematics and implementation
- **Setup overhead**: Requires initial coordinate transformations and configuration
- **Computational cost**: Additional overhead from regularization transformations
- **Physical interpretation**: Less intuitive physical meaning of state variables
- **Time coordinate**: Requires tracking both fictitious and physical time

## Implementation Details

The EDromo method in AstroPropagators.jl requires:

1. **Configuration**: RegularizedCoordinateConfig specifying units and options
2. **Initial transformation**: Convert from Cartesian to EDromo elements
3. **Integration**: Propagate in fictitious time φ
4. **Final transformation**: Convert back to desired coordinate system

## Usage Examples

!!! note "Fictitious time"
    EDromo integrates in a fictitious time variable `ϕ`, not physical time `t`.
    The `tspan` argument is therefore the range of `ϕ`, and a `ContinuousCallback`
    (via `end_EDromo_integration`) terminates the integration when the desired
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

ϕ₀ = compute_initial_phi(u0_cart, μ, config)
u0 = Array(EDromo(Cartesian(u0_cart), μ, ϕ₀, config))
tspan = (ϕ₀, ϕ₀ + 6π)

sol = propagate(
    EDromoPropagator(), u0, p, models, tspan, config;
    callback=end_EDromo_integration(86400.0, config),
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
u0 = Array(EDromo(Cartesian(u0_cart), μ, ϕ₀, config))
tspan = (ϕ₀, ϕ₀ + 6π)

sol = propagate(
    EDromoPropagator(), u0, p, models, tspan, config;
    callback=end_EDromo_integration(86400.0, config),
)
```

### Using `ODEProblem` Directly

```julia
using OrdinaryDiffEqAdamsBashforthMoulton, SciMLBase

# Using the same force model and config setup from above...
ϕ₀ = compute_initial_phi(u0_cart, μ, config)
u0 = Array(EDromo(Cartesian(u0_cart), μ, ϕ₀, config))
tspan = (ϕ₀, ϕ₀ + 6π)

f!(du, u, p, t) = EDromo_EOM!(du, u, p, t, models, config)

prob = ODEProblem(f!, u0, tspan, p)
sol = solve(
    prob, VCABM();
    abstol=1e-13, reltol=1e-13,
    callback=end_EDromo_integration(86400.0, config),
)
```

## Applications

EDromo is particularly valuable for:

### Asteroid and Comet Studies
- **Close approach prediction**: Accurate near-Earth object trajectories
- **Impact assessment**: Collision probability calculations
- **Long-term evolution**: Multi-orbital period propagation

### Spacecraft Mission Analysis
- **Highly elliptical orbits**: Molniya, GTO, and interplanetary trajectories
- **Gravity assist maneuvers**: Close planetary encounters
- **Orbit determination**: High-precision trajectory reconstruction

### Research Applications
- **Method validation**: Testing numerical integration techniques
- **Perturbation studies**: Effects on high-eccentricity orbits
- **Dynamical astronomy**: Long-term orbital evolution studies

## References

[1]: Baù, G., Bombardelli, C., Peláez, J., and Lorenzini, E. "Nonsingular orbital elements for special perturbations in the two-body problem." MNRAS 454(3), pp. 2890-2908, 2015.
[2]: Amato, D., Bombardelli, C., Baù, G., Morand, V., and Rosengren, A. J. "Non-averaged regularized formulations as an alternative to semi-analytical theories for highly perturbed orbits." Celestial Mechanics and Dynamical Astronomy, 131(5), pp. 1-36, 2019.
[3]: Roa, J., and Peláez, J. "Dynamics of a particle under the gravitational potential of a massive annulus." Celestial Mechanics and Dynamical Astronomy, 123(4), pp. 435-461, 2015.
[4]: Amato, D. "THALASSA: Orbit propagator for near-Earth and cislunar space." Astrophysics Source Code Library, ascl-1905, 2019. 