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

## Usage Example

```julia
using AstroPropagators
using AstroForceModels, AstroCoords
using ComponentArrays
using OrdinaryDiffEq

# Near-collision initial conditions
u0_cart = [6378.137 + 100, 0.0, 0.0, 0.0, 11.0, 0.0]  #km, km/s

# Create RegularizedCoordinateConfig from state  
μ = 3.986004415e5  # km³/s²
config = RegularizedCoordinateConfig(u0_cart, μ; flag_time=PhysicalTime())

# Convert to K-S coordinates
u0_ks = cart2KS(u0_cart, μ, config)

# Parameters
JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
p = ComponentVector(JD=JD, μ=μ)

# Force models
gravity_model = KeplerianGravityAstroModel(μ=μ)
models = DynamicsModel(gravity_model)

# Create ODE problem using K-S EOM
function dynamics!(du, u, p, t)
    KS_EOM!(du, u, p, t, models, config)
end

prob = ODEProblem(dynamics!, u0_ks, (0.0, 3600.0), p)
sol = solve(prob, Vern8(), abstol=1e-14, reltol=1e-14)

println("K-S propagation completed successfully")
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