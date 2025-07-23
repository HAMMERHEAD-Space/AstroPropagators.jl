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

### High-Eccentricity Orbit Propagation

```julia
using AstroPropagators
using AstroForceModels, AstroCoords
using ComponentArrays
using OrdinaryDiffEq

# High-eccentricity orbit (comet-like)
r0 = [20000.0, 0.0, 0.0]     # Aphelion position [km]
v0 = [0.0, 2.5, 0.0]         # Low velocity [km/s]
u0_cart = [r0; v0]

# Create RegularizedCoordinateConfig with different time options
μ = 3.986004415e5  # km³/s²

# Option 1: Physical time (simplest)
config_physical = RegularizedCoordinateConfig(u0_cart, μ; flag_time=PhysicalTime())

# Convert to EDromo elements (using linear time for this example)
u0_edromo = cart2EDromo(u0_cart, μ, 0.0, config_physical)

# Parameters
JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
p = ComponentVector(JD=JD, μ=μ)
tspan = (0.0, 86400.0 * 30)  # 30 days

# Define force models
gravity_model = KeplerianGravityAstroModel(μ=μ)
models = DynamicsModel(gravity_model)

# Create ODE problem using EDromo EOM
function dynamics!(du, u, p, t)
    EDromo_EOM!(du, u, p, t, models, config_linear)
end

prob = ODEProblem(dynamics!, u0_edromo, tspan, p)
sol = solve(prob, Vern8(), abstol=1e-12, reltol=1e-12)

println("EDromo propagation completed")
println("Final EDromo state: $(sol.u[end])")
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