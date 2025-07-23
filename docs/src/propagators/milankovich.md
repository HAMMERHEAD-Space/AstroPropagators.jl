# Milankovich Elements

The Milankovich elements are an alternative non-singular orbital element set that provides a robust parameterization of orbital motion. Named after Milutin Milanković, these elements are particularly useful for perturbed orbital dynamics and avoid the singularities present in classical Keplerian elements.

## Physical Description

Milankovich elements use a 7-dimensional state vector that combines angular momentum and eccentricity vectors with a true longitude:

### Element Definitions

The Milankovich state vector consists of 7 elements: `[hx, hy, hz, ex, ey, ez, L]`

#### **Angular Momentum Vector Components (hx, hy, hz)**
- **hx**: X-component of the specific angular momentum vector
- **hy**: Y-component of the specific angular momentum vector  
- **hz**: Z-component of the specific angular momentum vector

The angular momentum vector **H** = **r** × **v** is computed as:
```julia
H = cross(r, v)
```

**Physical meaning**: The angular momentum vector defines the orbital plane and its magnitude represents the orbital energy distribution.

#### **Eccentricity Vector Components (ex, ey, ez)**
- **ex**: X-component of the eccentricity vector
- **ey**: Y-component of the eccentricity vector
- **ez**: Z-component of the eccentricity vector

The eccentricity vector **e** is computed as:
```julia
e = cross(v, H) / μ - r / norm(r)
```

**Physical meaning**: The eccentricity vector points from the focus toward periapsis and its magnitude equals the orbital eccentricity.

#### **True Longitude (L)**
- **L**: True longitude = Ω + ω + f

**Physical meaning**: The sum of right ascension of ascending node, argument of periapsis, and true anomaly. This provides a continuous angular measure of position.

#### **Conservation Properties**

In the unperturbed two-body problem:
- **H** (angular momentum vector) is conserved
- **e** (eccentricity vector) is conserved  
- **L** (true longitude) varies linearly with time

## Applications

### **Perturbation Analysis**
- **Secular evolution**: Long-term orbital changes due to perturbations
- **Resonance studies**: Analysis of orbital resonances with other bodies
- **Stability analysis**: Investigation of orbital stability under perturbations

### **Spacecraft Mission Design**
- **Orbital maintenance**: Analysis of station-keeping requirements
- **Maneuver planning**: Understanding orbital changes from impulsive maneuvers
- **Mission optimization**: Long-term trajectory optimization

### **Celestial Mechanics Research**
- **Multi-body dynamics**: Analysis of complex gravitational interactions
- **Asteroid dynamics**: Study of asteroid orbital evolution
- **Exoplanet studies**: Analysis of exoplanetary system dynamics

### **Astrodynamics Education**
- **Advanced courses**: Graduate-level orbital mechanics
- **Research applications**: Specialized perturbation analysis
- **Theoretical studies**: Fundamental orbital dynamics research

## Usage Examples

### Basic Orbital Propagation

```julia
using AstroPropagators
using AstroForceModels, AstroCoords
using ComponentArrays
using OrdinaryDiffEq

# Initial conditions (can be any orbit type)
r0 = [6378.137 + 800, 0.0, 0.0]  # Position [km]
v0 = [0.0, 7.450, 0.0]           # Velocity [km/s]
u0_cart = [r0; v0]

# Convert to Milankovich elements
cart_state = Cartesian(u0_cart, μ=3.986004415e5)
mil_state = Milankovich(cart_state)
u0_milan = [mil_state.hx, mil_state.hy, mil_state.hz,
            mil_state.ex, mil_state.ey, mil_state.ez, mil_state.L]

# Parameters
JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
p = ComponentVector(JD=JD, μ=3.986004415e5)  # km³/s²

# Force models
gravity_model = GravityHarmonicsAstroModel(
    gravity_model = gravity_data,
    eop_data = eop_data,
    order = 4, degree = 4
)
models = DynamicsModel(gravity_model)

# Create ODE problem using Milankovich EOM
function dynamics!(du, u, p, t)
    Milankovich_EOM!(du, u, p, t, models)
end

prob = ODEProblem(dynamics!, u0_milan, (0.0, 86400.0), p)
sol = solve(prob, Vern8(), abstol=1e-12, reltol=1e-12)

println("Milankovich propagation completed")
```

## References

[1]: Milanković, M. (1939). "Théorie mathématique des phénomènes thermiques produits par la radiation solaire". Gauthier-Villars, Paris.
[2]: Brouwer, D. (1959). "Solution of the problem of artificial satellite theory without drag". Astronomical Journal, 64, 378-397.
[3]: Kozai, Y. (1959). "The motion of a close earth satellite". Astronomical Journal, 64, 367-377.
[4]: Allan, R. R. (1962). "Third-body perturbations in satellite theory". Royal Aircraft Establishment Technical Report.
[5]:  Cefola, P. J. (1972). "Equinoctial orbital elements - Application to artificial satellite orbits". AIAA/AAS Astrodynamics Conference, Paper 72-937.
[6]: https://www.researchgate.net/publication/263032883_On_the_Milankovitch_orbital_elements_for_perturbed_Keplerian_motion#fullTextFileContent
[7]: https://www.researchgate.net/publication/338021512_Utilizing_Maximum_Entropy_Spectral_Analysis_MESA_to_identify_Milankovitch_cycles_in_Lower_Member_of_Miocene_Zhujiang_Formation_in_north_slope_of_Baiyun_Sag_Pearl_River_Mouth_Basin_South_China_Sea