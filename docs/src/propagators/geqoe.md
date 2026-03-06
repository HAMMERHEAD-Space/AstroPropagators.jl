# Generalized Equinoctial Orbital Elements (GEqOE)

The Generalized Equinoctial Orbital Elements (GEqOE) are a set of six orbital elements that generalize the classical equinoctial elements when some or all of the perturbing forces acting on the propagated body are derived from a disturbing potential. The perturbing potential is embedded directly into the element definitions, producing a non-osculating ellipse that better captures the perturbed dynamics, leading to dramatically improved propagation performance.

## Physical Description

### Element Definitions

The GEqOE state vector consists of 6 elements: `[ν, p₁, p₂, L, q₁, q₂]`

#### **Generalized Mean Motion (ν)**
- Defined as `ν = (1/μ)(-2E)^{3/2}` where `E` is the total energy (including the disturbing potential)
- Reduces to the classical mean motion `n` when perturbations are zero
- Constant along unperturbed Keplerian motion

#### **Generalized Eccentricity Components (p₁, p₂)**
- `p₁ = g sin(Ψ)`: generalized analog of the equinoctial element `h = e sin(ω + Ω)`
- `p₂ = g cos(Ψ)`: generalized analog of the equinoctial element `k = e cos(ω + Ω)`
- Where `g` is the magnitude of the generalized Laplace vector and `Ψ` is the generalized longitude of pericenter
- Non-singular for circular orbits (`g = 0`)

#### **Generalized Mean Longitude (L)**
- `L = M + Ψ` where `M` is the generalized mean anomaly
- Varies linearly with time along Keplerian motion: `dL/dt = ν`

#### **Inclination Components (q₁, q₂)**
- `q₁ = tan(i/2) sin(Ω)`: identical to the classical equinoctial element `p`
- `q₂ = tan(i/2) cos(Ω)`: identical to the classical equinoctial element `q`
- Non-singular for equatorial orbits (`i = 0`)
- Singular for retrograde equatorial orbits (`i = π`)

### Key Properties

- **Non-singular** for circular and equatorial orbits
- **Singular** for retrograde equatorial orbits (`i = π`) and rectilinear motion
- **Defined for negative total energy** (`E < 0`, i.e., bound orbits)
- **Perturbing potential embedded** in element definitions via the generalized angular momentum `c`, generalized semi-major axis `a`, and generalized eccentricity `g`

### Conservation Properties

In the unperturbed two-body problem (with no perturbations):
- **ν** (generalized mean motion) is conserved
- **p₁, p₂** (generalized eccentricity components) are conserved
- **L** (generalized mean longitude) varies linearly: `L̇ = ν`
- **q₁, q₂** (inclination components) are conserved

## Configuration

GEqOE requires a `RegularizedCoordinateConfig` that specifies the perturbing potential `W`:

```julia
config = RegularizedCoordinateConfig(; W=W)
```

For a **Keplerian orbit** (no perturbations), set `W = 0`:
```julia
config = RegularizedCoordinateConfig(; W=0.0)
```

For **perturbed orbits**, compute the initial perturbing potential:
```julia
W = potential(Cartesian(u0), p, 0.0, gravity_model) -
    potential(Cartesian(u0), p, 0.0, KeplerianGravityAstroModel(; μ=μ))
config = RegularizedCoordinateConfig(; W=W)
```

## Usage Examples

### Basic Orbital Propagation

```julia
using AstroPropagators
using AstroForceModels, AstroCoords
using ComponentArrays
using OrdinaryDiffEqAdamsBashforthMoulton

# Initial Cartesian state [km, km/s]
u0_cart = [-1076.225, -6765.896, -332.309, 9.357, -3.312, -1.188]

# Parameters
JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
μ = 398600.4418
p = ComponentVector(; JD=JD, μ=μ)

# Force model
grav_model = KeplerianGravityAstroModel()
model_list = CentralBodyDynamicsModel(grav_model)

# Configuration (W = 0 for Keplerian)
config = RegularizedCoordinateConfig(; W=0.0)

# Convert to GEqOE
u0_geqoe = Array(GEqOE(Cartesian(u0_cart), μ, config))

# Propagate
EOM!(du, u, p, t) = GEqOE_EOM!(du, u, p, t, model_list, config)
prob = ODEProblem(EOM!, u0_geqoe, (0.0, 86400.0), p)
sol = solve(prob, VCABM(); abstol=1e-13, reltol=1e-13)

# Convert back to Cartesian
final_cart = Cartesian(GEqOE(sol.u[end]), μ, config)
```

### High-Fidelity Propagation

```julia
# Gravity model with J2+ harmonics
grav_model = GravityHarmonicsAstroModel(;
    gravity_model=grav_coeffs, eop_data=eop_data, order=36, degree=36
)

# Third-body and surface forces
sun_model = ThirdBodyModel(; body=SunBody(), eop_data=eop_data)
moon_model = ThirdBodyModel(; body=MoonBody(), eop_data=eop_data)
srp_model = SRPAstroModel(; satellite_srp_model=CannonballFixedSRP(0.2), ...)
drag_model = DragAstroModel(; satellite_drag_model=CannonballFixedDrag(0.2), ...)

model_list = CentralBodyDynamicsModel(
    grav_model, (sun_model, moon_model, srp_model, drag_model)
)

# Compute initial perturbing potential
W = potential(Cartesian(u0_cart), p, 0.0, grav_model) -
    potential(Cartesian(u0_cart), p, 0.0, KeplerianGravityAstroModel(; μ=μ))
config = RegularizedCoordinateConfig(; W=W)

# Convert and propagate
u0_geqoe = Array(GEqOE(Cartesian(u0_cart), μ, config))
EOM!(du, u, p, t) = GEqOE_EOM!(du, u, p, t, model_list, config)
prob = ODEProblem(EOM!, u0_geqoe, (0.0, 86400.0), p)
sol = solve(prob, VCABM(); abstol=1e-13, reltol=1e-13)
```

## Performance Characteristics

The GEqOE formulation shows **remarkable improvement** compared to the classical equinoctial elements and Cowell's method, particularly for:

- **Conservative perturbations** (e.g., geopotential harmonics): Several orders of magnitude improvement in position accuracy
- **Mixed perturbations** (geopotential + third-body): Substantial improvement maintained when non-conservative forces are present
- **Long-term propagation**: Slower error growth compared to classical methods

The improvement comes from embedding the disturbing potential directly into the element definitions, reducing the variation rate of the elements under perturbation.

## References

[1] Baù, G., Hernando-Ayuso, J., & Bombardelli, C. (2021). "A generalization of the equinoctial orbital elements." Celestial Mechanics and Dynamical Astronomy, 133(9), 1-32. [arXiv:2105.04424](https://arxiv.org/abs/2105.04424)

[2] Broucke, R. A., & Cefola, P. J. (1972). "On the equinoctial orbit elements." Celestial Mechanics, 5, 303-310.

[3] Walker, M. J. H., Ireland, B., & Owens, J. (1985). "A set of modified equinoctial orbit elements." Celestial Mechanics, 36, 409-419.
