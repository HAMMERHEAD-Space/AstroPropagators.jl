# Unified State Model (USM)

The Unified State Model (USM) is a family of non-singular orbital element sets that combine velocity hodograph components with different attitude parameterizations. Developed by Shepperd in the 1980s, USM provides a unified framework for orbital propagation that avoids the singularities present in classical Keplerian elements.

## Physical Description

USM elements use a hybrid approach combining:
1. **Velocity hodograph components** (C, Rf1, Rf2) that describe the orbital motion
2. **Attitude parameters** that describe the orbital plane orientation

### Element Definitions

#### **USM7 Elements** `[C, Rf1, Rf2, ϵO1, ϵO2, ϵO3, η0]`

**Velocity Hodograph Components:**
- **C**: Velocity hodograph component normal to the radial vector in the orbital plane
  ```julia
  C = √(μ / (a * (1 - e²)))
  ```
- **Rf1**: Velocity hodograph component 90° ahead of eccentricity vector along intermediate rotating frame X-axis
  ```julia
  Rf1 = -e * C * sin(Ω + ω)
  ```
- **Rf2**: Velocity hodograph component 90° ahead of eccentricity vector along intermediate rotating frame Y-axis
  ```julia
  Rf2 = e * C * cos(Ω + ω)
  ```

**Quaternion Components (Euler Parameters):**
- **ϵO1**: First imaginary quaternion component
  ```julia
  ϵO1 = sin(i/2) * cos((Ω - u)/2)
  ```
- **ϵO2**: Second imaginary quaternion component
  ```julia
  ϵO2 = sin(i/2) * sin((Ω - u)/2)
  ```
- **ϵO3**: Third imaginary quaternion component
  ```julia
  ϵO3 = cos(i/2) * sin((Ω + u)/2)
  ```
- **η0**: Real quaternion component
  ```julia
  η0 = cos(i/2) * cos((Ω + u)/2)
  ```

Where `u = ω + f` is the argument of latitude.

#### **USM6 Elements** `[C, Rf1, Rf2, σ1, σ2, σ3]`

USM6 uses the same velocity hodograph components as USM7 but replaces quaternions with Modified Rodrigues Parameters (MRPs):

- **C, Rf1, Rf2**: Same velocity hodograph components as USM7
- **σ1, σ2, σ3**: Modified Rodrigues Parameters for attitude representation

#### **USMEM Elements** `[C, Rf1, Rf2, a1, a2, a3]`

USMEM uses the same velocity hodograph components but with exponential mapping for attitude:

- **C, Rf1, Rf2**: Same velocity hodograph components as USM7
- **a1, a2, a3**: Exponential mapping components for attitude representation

### Mathematical Relationships

#### **Orbital Parameters from USM Elements**

```julia
# From velocity hodograph components
R = √(Rf1² + Rf2²)
e = R / C
a = μ / (2 * C * ve2 - (ve1² + ve2²))

# From quaternion components (USM7)
i = acos(1 - 2 * (ϵO1² + ϵO2²))
Ω = atan((ϵO1 * ϵO3 + ϵO2 * η0) / √((ϵO1² + ϵO2²) * (η0² + ϵO3²)),
         (ϵO1 * η0 - ϵO2 * ϵO3) / √((ϵO1² + ϵO2²) * (η0² + ϵO3²)))
```

#### **Velocity Hodograph Interpretation**

The velocity hodograph components represent the velocity vector in a special coordinate frame:
- **C**: Magnitude of the velocity component perpendicular to the position vector
- **Rf1, Rf2**: Components of the velocity in the orbital plane relative to the eccentricity vector

## Applications

### **General Orbital Propagation**
- **Multi-purpose use**: Suitable for all orbit types and eccentricities
- **Long-term integration**: Stable for extended propagation periods
- **Perturbation analysis**: Well-suited for studying orbital evolution

### **Spacecraft Mission Design**
- **Mission planning**: Robust orbital propagation for mission design
- **Trajectory optimization**: Stable framework for optimization algorithms
- **Attitude coupling**: Natural integration with attitude dynamics

### **Research Applications**
- **Algorithm development**: Testing new propagation techniques
- **Method comparison**: Benchmarking against other element sets
- **Theoretical studies**: Exploring unified orbital representations

## Usage Examples

### USM7 — Keplerian Propagation

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

u0 = Array(USM7(Cartesian(u0_cart), μ))
sol = propagate(USM7Propagator(), u0, p, models, tspan)
```

### USM7 — High-Fidelity Propagation

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

u0 = Array(USM7(Cartesian(u0_cart), μ))
sol = propagate(USM7Propagator(), u0, p, models, tspan)
```

### USM6 and USMEM

USM6 and USMEM share the same velocity-hodograph core but use different attitude
parameterizations. Just swap the coordinate type and propagator:

```julia
# USM6 — Modified Rodrigues Parameters
u0_usm6 = Array(USM6(Cartesian(u0_cart), μ))
sol = propagate(USM6Propagator(), u0_usm6, p, models, tspan)

# USMEM — Exponential Mapping
u0_usmem = Array(USMEM(Cartesian(u0_cart), μ))
sol = propagate(USMEMPropagator(), u0_usmem, p, models, tspan)
```

### Using `ODEProblem` Directly

```julia
using OrdinaryDiffEqVerner, SciMLBase

# Using the same force model setup from above...
u0 = Array(USM7(Cartesian(u0_cart), μ))

f!(du, u, p, t) = USM7_EOM!(du, u, p, t, models)

prob = ODEProblem(f!, u0, tspan, p)
sol = solve(prob, Vern9(); abstol=1e-13, reltol=1e-13)

# Similarly for USM6 and USMEM:
# f!(du, u, p, t) = USM6_EOM!(du, u, p, t, models)
# f!(du, u, p, t) = USMEM_EOM!(du, u, p, t, models)
```

## References

[1]: Shepperd, S. W. (1985). "Universal Keplerian state transition matrix". Celestial Mechanics, 35(2), 129-144.
[2]: Shepperd, S. W. (1978). "Quaternion from rotation matrix". Journal of Guidance and Control, 1(3), 223-224.
[3]: Shuster, M. D. (1993). "A survey of attitude representations". Journal of the Astronautical Sciences, 41(4), 439-517.
[4]: Schaub, H., & Junkins, J. L. (2003). "Analytical mechanics of space systems". AIAA.
[5]: Vallado, D. A. (2013). "Fundamentals of astrodynamics and applications". 4th ed., Microcosm Press.
[6]: Van den Broeck, M. (2017). "An Approach to Generalizing Taylor Series Integration for Low-Thrust Trajectories". Delft University of Technology.
[7]: https://link.springer.com/article/10.1007/s10569-011-9396-5#:~:text=The%20Unified%20State%20Model%20is,a%20set%20of%20seven%20elements.
[8]: https://link.springer.com/article/10.1007/BF01227757