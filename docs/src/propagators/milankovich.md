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

u0 = Array(Milankovich(Cartesian(u0_cart), μ))
sol = propagate(MilankovichPropagator(), u0, p, models, tspan)
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

u0 = Array(Milankovich(Cartesian(u0_cart), μ))
sol = propagate(MilankovichPropagator(), u0, p, models, tspan)
```

### Using `ODEProblem` Directly

```julia
using OrdinaryDiffEqVerner, SciMLBase

# Using the same force model setup from above...
u0 = Array(Milankovich(Cartesian(u0_cart), μ))

f!(du, u, p, t) = Milankovich_EOM!(du, u, p, t, models)

prob = ODEProblem(f!, u0, tspan, p)
sol = solve(prob, Vern9(); abstol=1e-13, reltol=1e-13)
```

## References

[1]: Milanković, M. (1939). "Théorie mathématique des phénomènes thermiques produits par la radiation solaire". Gauthier-Villars, Paris.
[2]: Brouwer, D. (1959). "Solution of the problem of artificial satellite theory without drag". Astronomical Journal, 64, 378-397.
[3]: Kozai, Y. (1959). "The motion of a close earth satellite". Astronomical Journal, 64, 367-377.
[4]: Allan, R. R. (1962). "Third-body perturbations in satellite theory". Royal Aircraft Establishment Technical Report.
[5]:  Cefola, P. J. (1972). "Equinoctial orbital elements - Application to artificial satellite orbits". AIAA/AAS Astrodynamics Conference, Paper 72-937.
[6]: https://www.researchgate.net/publication/263032883_On_the_Milankovitch_orbital_elements_for_perturbed_Keplerian_motion#fullTextFileContent
[7]: https://www.researchgate.net/publication/338021512_Utilizing_Maximum_Entropy_Spectral_Analysis_MESA_to_identify_Milankovitch_cycles_in_Lower_Member_of_Miocene_Zhujiang_Formation_in_north_slope_of_Baiyun_Sag_Pearl_River_Mouth_Basin_South_China_Sea