# Modified Equinoctial Elements

The Modified Equinoctial Elements (MEE) propagator uses the Gauss Variational Equations to integrate perturbations directly in the modified equinoctial orbital element set `[p, f, g, h, k, L]`. This formulation is singularity-free for circular and equatorial orbits, making it the natural choice for low-thrust trajectory optimization and guidance algorithms such as Q-Law.

## Physical Description

The modified equinoctial elements are defined in terms of the classical Keplerian elements:

```
p = a(1 − e²)         semi-latus rectum
f = e cos(ω + Ω)      eccentricity × cos(longitude of periapsis)
g = e sin(ω + Ω)      eccentricity × sin(longitude of periapsis)
h = tan(i/2) cos(Ω)   half-angle inclination × cos(RAAN)
k = tan(i/2) sin(Ω)   half-angle inclination × sin(RAAN)
L = Ω + ω + ν         true longitude
```

The Gauss Variational Equations in MEE under perturbing accelerations `[Fr, Fθ, Fh]` in the RTN frame are:

```
dp/dt = (2p/q) √(p/μ) Fθ
df/dt = √(p/μ) [q sin(L) Fr + ((q+1)cos(L) + f) Fθ − g(h sin(L) − k cos(L)) Fh] / q
dg/dt = √(p/μ) [−q cos(L) Fr + ((q+1)sin(L) + g) Fθ + f(h sin(L) − k cos(L)) Fh] / q
dh/dt = √(p/μ) s² cos(L) Fh / (2q)
dk/dt = √(p/μ) s² sin(L) Fh / (2q)
dL/dt = √(μp) q²/p² + √(p/μ) (h sin(L) − k cos(L)) Fh / q
```

where `q = 1 + f cos(L) + g sin(L)` and `s² = 1 + h² + k²`.

## Advantages and Disadvantages

### Advantages
- **Singularity-free**: No singularities for circular (`e = 0`) or equatorial (`i = 0`) orbits, unlike classical Keplerian elements
- **Natural for low-thrust**: Semi-latus rectum `p` and true longitude `L` vary smoothly under continuous thrust, avoiding the rapid osculations seen in classical elements
- **Compact state**: Only 6 elements (same as classical Keplerian) with no quaternion normalization constraints
- **Efficient GVE**: The `dp/dt` equation depends only on the tangential acceleration, simplifying analytical insight
- **AD-compatible**: All operations are smooth and differentiable

### Disadvantages
- **Retrograde singularity**: `h` and `k` are singular at `i = π` (retrograde orbits); the retrograde factor variant can avoid this at the cost of an `i = 0` singularity
- **Indirect position**: Converting to/from Cartesian requires trigonometric evaluations
- **True longitude**: `L` grows secularly, which may require modular arithmetic for long propagations

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

u0 = Array(ModEq(Cartesian(u0_cart), μ))
sol = propagate(ModifiedEquinoctialPropagator(), u0, p, models, tspan)
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

u0 = Array(ModEq(Cartesian(u0_cart), μ))
sol = propagate(ModifiedEquinoctialPropagator(), u0, p, models, tspan)
```

### Using `ODEProblem` Directly

```julia
using OrdinaryDiffEqVerner, SciMLBase

# Using the same force model setup from above...
u0 = Array(ModEq(Cartesian(u0_cart), μ))

# In-place (what propagate! uses internally)
f!(du, u, p, t) = ModEq_EOM!(du, u, p, t, models)
prob = ODEProblem(f!, u0, tspan, p)
sol = solve(prob, Vern9(); abstol=1e-13, reltol=1e-13)

# Out-of-place (what propagate uses internally)
f(u, p, t) = ModEq_EOM(u, p, t, models)
prob = ODEProblem{false}(f, u0, tspan, p)
sol = solve(prob, Vern9(); abstol=1e-13, reltol=1e-13)
```

## Applications

The Modified Equinoctial Elements propagator is particularly suitable for:

### Low-Thrust Trajectory Optimization
- **Q-Law guidance**: MEE are the native element set for Q-Law Lyapunov feedback control
- **Indirect optimization**: Smooth element variations enable gradient-based optimization
- **Continuous thrust modeling**: No discontinuities from element singularities during spiral maneuvers

### Mission Analysis
- **GEO transfers**: Singularity-free at zero inclination and eccentricity
- **Low-Earth orbit maintenance**: Handles near-circular orbits without numerical issues
- **Constellation design**: Smooth propagation across a range of inclinations and eccentricities

### Cross-Formulation Validation
- **Independent verification**: Compare against Cowell, GaussVE, and other formulations
- **Element-space analysis**: Direct insight into how perturbations affect orbit shape

## References

[1]: Walker, M.J.H., Ireland, B. & Owens, J. "A Set of Modified Equinoctial Orbit Elements." Celestial Mechanics, 36, 409–419, 1985.
[2]: Battin, R.H. "An Introduction to the Mathematics and Methods of Astrodynamics." AIAA Education Series, 1999.
[3]: Broucke, R.A. & Cefola, P.J. "On the Equinoctial Orbit Elements." Celestial Mechanics, 5, 303–310, 1972.
[4]: Varga, G. & Pérez, J.M.S. "Many-Revolution Low-Thrust Orbit Transfer Computation Using Equinoctial Q-Law Including J2 and Eclipse Constraints." ICATT, 2015.
