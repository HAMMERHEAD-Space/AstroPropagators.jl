# Stiefel-Scheifele Method

The Stiefel-Scheifele method is a canonical regularization technique that transforms the orbital motion problem using universal variables and canonical transformations. This method provides excellent numerical stability for eccentric orbits while preserving the Hamiltonian structure of the dynamics.

## Physical Description

The Stiefel-Scheifele method employs:

- **Universal variables**: Elements that work for all orbit types (elliptic, parabolic, hyperbolic)
- **Canonical transformation**: Preserves the symplectic structure of Hamiltonian mechanics
- **Regularized time**: Uses a transformed time coordinate for uniform integration steps
- **Energy preservation**: Maintains conservation properties better than non-canonical methods



## Usage Example

```julia
using AstroPropagators
using AstroForceModels, AstroCoords
using ComponentArrays
using OrdinaryDiffEq

# Initial conditions (works for any eccentricity)
r0 = [6378.137 + 500, 0.0, 0.0]  # Position [km]
v0 = [0.0, 9.5, 0.0]             # Velocity [km/s] (hyperbolic)
u0_cart = [r0; v0]

# Create RegularizedCoordinateConfig from state
μ = 3.986004415e5  # km³/s²
config = RegularizedCoordinateConfig(u0_cart, μ; flag_time=PhysicalTime())

# Convert to Stiefel-Scheifele elements
u0_ss = cart2StiefelScheifele(u0_cart, μ, 0.0, config)

# Parameters
JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
p = ComponentVector(JD=JD, μ=μ)

# Force models
gravity_model = KeplerianGravityAstroModel(μ=μ)
models = DynamicsModel(gravity_model)

# Create ODE problem using Stiefel-Scheifele EOM
function dynamics!(du, u, p, t)
    StiefelScheifele_EOM!(du, u, p, t, models, config)
end

prob = ODEProblem(dynamics!, u0_ss, (0.0, 86400.0), p)
sol = solve(prob, Vern8(), abstol=1e-12, reltol=1e-12)

println("Stiefel-Scheifele propagation completed")
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