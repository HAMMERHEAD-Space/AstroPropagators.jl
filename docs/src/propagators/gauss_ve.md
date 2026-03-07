# Gauss Variational Equations

The Gauss Variational Equations method, also known as Gauss Planetary Equations, is a classical approach to orbital propagation that integrates perturbations directly in Keplerian orbital elements. Developed by Carl Friedrich Gauss, this method is particularly efficient when perturbative forces are small compared to the central gravitational force.

## Physical Description

The Gauss method is based on the principle that in the absence of perturbations, Keplerian orbital elements remain constant. When perturbative forces are present, these elements vary slowly with time according to the Gauss planetary equations:

```
da/dt = (2a²/h) * [e*sin(f)*R + (a*(1-e²)/r)*T]
de/dt = (1/h) * [sin(f)*R + ((cos(f) + cos(E))*T)]
di/dt = (r*cos(θ)/h) * W
dΩ/dt = (r*sin(θ)/(h*sin(i))) * W
dω/dt = (1/(h*e)) * [-cos(f)*R + sin(f)*(1 + r/a)*T] - (r*sin(θ)*cos(i)/(h*sin(i))) * W
df/dt = h/r² + (1/(h*e)) * [cos(f)*R - sin(f)*(1 + r/a)*T]
```

Where:
- `a, e, i, Ω, ω, f` are the six Keplerian elements
- `R, T, W` are radial, tangential, and normal force components
- `h = √(μa(1-e²))` is specific angular momentum
- `r = a(1-e²)/(1+e*cos(f))` is orbital radius
- `θ = ω + f` is argument of latitude

## Coordinate System

The method uses the RTN (Radial-Tangential-Normal) coordinate system:

- **Radial (R)**: Along the position vector from center to spacecraft
- **Tangential (T)**: In the orbital plane, perpendicular to radial, in direction of motion
- **Normal (W)**: Perpendicular to orbital plane, completing right-handed system

## Advantages and Disadvantages

### Advantages
- **Physical Insight**: Direct evolution of meaningful orbital parameters
- **Computational Efficiency**: Fewer force evaluations for small perturbations
- **Element Preservation**: Natural representation of orbit shape and orientation
- **Analytical Understanding**: Clear relationship between forces and orbital changes
- **Perturbation Theory**: Ideal for analytical and semi-analytical studies

### Disadvantages

- **Large Perturbations**: Less accurate when perturbative forces are large
- **Computational Overhead**: Coordinate transformations and trigonometric functions
- **Element Coupling**: Strong coupling between elements can cause numerical issues

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

u0 = Array(Keplerian(Cartesian(u0_cart), μ))
sol = propagate(GaussVEPropagator(), u0, p, models, tspan)
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

u0 = Array(Keplerian(Cartesian(u0_cart), μ))
sol = propagate(GaussVEPropagator(), u0, p, models, tspan)
```

### Using `ODEProblem` Directly

```julia
using OrdinaryDiffEqVerner, SciMLBase

# Using the same force model setup from above...
u0 = Array(Keplerian(Cartesian(u0_cart), μ))

f!(du, u, p, t) = GaussVE_EOM!(du, u, p, t, models)

prob = ODEProblem(f!, u0, tspan, p)
sol = solve(prob, Vern9(); abstol=1e-13, reltol=1e-13)
```

## Applications

The Gauss VE method is particularly suitable for:

### Orbital Mechanics Research
- **Perturbation theory validation**: Comparing numerical with analytical results
- **Resonance studies**: Understanding mean motion and secular resonances
- **Orbital element sensitivity**: Analyzing response to different perturbations

### Mission Analysis
- **Long-term orbit evolution**: Multi-year propagation studies
- **Station-keeping requirements**: Fuel budget estimation for orbit maintenance
- **Orbit lifetime analysis**: Natural decay due to perturbations

### Educational Applications
- **Teaching orbital mechanics**: Direct visualization of element evolution
- **Parameter studies**: Understanding effect of different perturbations
- **Method comparison**: Demonstrating different propagation approaches

## References

[1]: Gauss, C. F. "Theoria motus corporum coelestium in sectionibus conicis solem ambientium." 1809.
[2]: Vallado, David A. "Fundamentals of Astrodynamics and Applications." 4th ed., 2013.
[3]: Battin, Richard H. "An Introduction to the Mathematics and Methods of Astrodynamics." AIAA, 1999.
[4]: Brouwer, D., and G. M. Clemence. "Methods of Celestial Mechanics." Academic Press, 1961.
[5]: Roy, A. E. "Orbital Motion." 4th ed., Institute of Physics Publishing, 2005. 