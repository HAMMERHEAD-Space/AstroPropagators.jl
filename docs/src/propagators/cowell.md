# Cowell Method

The Cowell method is the most direct approach to orbital propagation, integrating the equations of motion in Cartesian coordinates. Named after Philip H. Cowell, this method solves the second-order differential equation of motion by directly integrating position and velocity vectors in an inertial reference frame.

## Physical Description

The Cowell method is based on Newton's second law applied to orbital motion:

```
d²r/dt² = a_total = Σ F_i / m
```

Where:
- `r` is the position vector
- `a_total` is the total acceleration from all forces
- `F_i` represents individual force contributions
- `m` is the spacecraft mass

The method converts the second-order differential equation into a system of first-order equations:

```
dr/dt = v
dv/dt = a_total(r, v, t)
```

This results in a 6-dimensional state vector `[rx, ry, rz, vx, vy, vz]` that is integrated over time.

## Advantages and Disadvantages

### Advantages
- **Conceptual Simplicity**: Direct physical interpretation of position and velocity
- **Force Model Flexibility**: Can handle any combination of force models without modification
- **High Precision**: No approximations in the coordinate system representation
- **Numerical Stability**: Well-conditioned for most orbit types
- **Universal Applicability**: Works for all orbit geometries and eccentricities

### Disadvantages
- **Computational Cost**: Requires frequent force evaluations
- **Step Size Sensitivity**: Small integration steps needed for high-eccentricity orbits
- **No Inherent Regularization**: May struggle with near-collision trajectories

## Usage Examples

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
tspan = (0.0, 86400.0)

sol = propagate(CowellPropagator(), u0_cart, p, models, tspan)
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
tspan = (0.0, 86400.0)

sol = propagate(CowellPropagator(), u0_cart, p, models, tspan)
```

### Using `ODEProblem` Directly

For full control over the integration — custom callbacks, event handling, or composing
the EOM closure into a larger system — you can bypass `propagate` and build the
`ODEProblem` yourself.

```julia
using OrdinaryDiffEqVerner, SciMLBase

# Using the same force model setup from above...

# Build the EOM closure — captures models (and config for regularized propagators)
f!(du, u, p, t) = Cowell_EOM!(du, u, p, t, models)

# Or equivalently, through the dispatch layer:
# f!(du, u, p, t) = eom!(CowellPropagator(), du, u, p, t, models)

prob = ODEProblem(f!, u0_cart, tspan, p)
sol = solve(prob, Vern9(); abstol=1e-13, reltol=1e-13)
```

This is exactly what `propagate` does internally, so the two approaches produce
identical results.

## Optimal Use Cases
- **High-precision applications**: Precise orbit determination, collision analysis
- **Complex force environments**: Multiple simultaneous perturbations
- **Irregular force models**: Non-conservative forces, thrust maneuvers
- **Validation studies**: Reference solutions for other propagation methods

## Applications

The Cowell method is particularly well-suited for:

### Precise Orbit Determination
- **Satellite tracking**: High-accuracy orbit reconstruction from observations
- **Laser ranging analysis**: Centimeter-level precision requirements
- **Formation flying**: Relative motion with tight accuracy constraints

### Mission Operations
- **Collision avoidance**: Accurate prediction of close approaches
- **Maneuver planning**: Precise delta-V calculations and execution
- **Ground station scheduling**: Accurate visibility predictions

### Scientific Applications
- **Gravity field determination**: Recovery of gravitational coefficients
- **Atmospheric research**: Precise drag coefficient estimation
- **Relativity tests**: High-precision tests of general relativity

### Validation and Benchmarking
- **Reference solutions**: Gold standard for validating other methods
- **Algorithm development**: Testing new propagation techniques
- **Monte Carlo studies**: Statistical analysis requiring many propagations

## References

[1]: Cowell, P. H., and A. C. D. Crommelin. "Investigation of the motion of Halley's comet from 1759 to 1910." Greenwich Observations, 1910.
[2]: Vallado, David A. "Fundamentals of Astrodynamics and Applications." 4th ed., 2013.
[3]: Battin, Richard H. "An Introduction to the Mathematics and Methods of Astrodynamics." AIAA, 1999.
[4]: Montenbruck, Oliver, and Eberhard Gill. "Satellite Orbits: Models, Methods and Applications." Springer, 2000. 