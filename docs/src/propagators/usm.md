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

### USM7 Propagation

```julia
using AstroPropagators
using AstroForceModels, AstroCoords
using ComponentArrays
using OrdinaryDiffEq

# Initial Cartesian state
r0 = [6378.137 + 408, 0.0, 0.0]  # Position [km]
v0 = [0.0, 7.660, 0.0]           # Velocity [km/s]
u0_cart = [r0; v0]

# Convert to USM7
cart_state = Cartesian(u0_cart, μ=3.986004415e5)  # km, km/s, km³/s²
usm_state = USM7(cart_state)
u0_usm7 = [usm_state.C, usm_state.Rf1, usm_state.Rf2,
           usm_state.ϵO1, usm_state.ϵO2, usm_state.ϵO3, usm_state.η0]

# Parameters
JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
p = ComponentVector(JD=JD, μ=3.986004415e5)  # km³/s²

# Force models
gravity_model = KeplerianGravityAstroModel(μ=3.986004415e5)
models = DynamicsModel(gravity_model)

# Create ODE problem using USM7 EOM
function dynamics!(du, u, p, t)
    USM7_EOM!(du, u, p, t, models)
end

prob = ODEProblem(dynamics!, u0_usm7, (0.0, 86400.0), p)
sol = solve(prob, Vern8(), abstol=1e-12, reltol=1e-12)
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