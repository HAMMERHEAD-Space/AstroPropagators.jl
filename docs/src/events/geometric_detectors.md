# Geometric Detectors

Geometric and geographic event detectors trigger on spatial quantities: altitude, geocentric latitude, sub-satellite longitude, and beta angle.

## Altitude Crossing

`altitude_condition` is zero when the spherical altitude crosses a target value. Uses `|r| - (R_body + h_target)` as the g-function.

!!! note
    This is spherical altitude (distance from body center minus body radius), not geodetic altitude. Sufficiently accurate for event detection purposes.

### Terminate at atmospheric re-entry

```julia
cond = altitude_condition(Cartesian, 120.0)
cb = ContinuousCallback(
    cond,
    nothing,                                          # ignore ascent through 120 km
    (integrator) -> terminate!(integrator),            # terminate on descent
)

sol = propagate(CowellPropagator(), u0_cart, p, models, tspan; callback=cb)
```

### Log all passages through ISS altitude

```julia
log = []
cond = altitude_condition(Cartesian, 400.0)
cb = ContinuousCallback(
    cond,
    (integrator) -> push!(log, (t=integrator.t, dir=:ascending)),
    (integrator) -> push!(log, (t=integrator.t, dir=:descending)),
)
```

### Non-Earth body

Override `R_body` for other central bodies:

```julia
# Lunar orbit — detect 100 km altitude crossing (R_Moon ≈ 1737.4 km)
cond = altitude_condition(Cartesian, 100.0; R_body=1737.4)

# Mars orbit (R_Mars ≈ 3396.2 km)
cond = altitude_condition(Cartesian, 200.0; R_body=3396.2)
```

## Geocentric Latitude Crossing

`latitude_condition` is zero when geocentric latitude crosses a target. Uses `asin(r_z / |r|) - lat_target`.

### Detect equator crossings

```julia
cond = latitude_condition(Cartesian, 0.0)
cb = ContinuousCallback(cond, affect!, affect!)
```

### Detect high-latitude passages

```julia
# North and south ±60° latitude
cond_n = latitude_condition(Cartesian, deg2rad(60))
cond_s = latitude_condition(Cartesian, deg2rad(-60))

cb = CallbackSet(
    ContinuousCallback(cond_n, north_affect!, north_affect!),
    ContinuousCallback(cond_s, south_affect!, south_affect!),
)
```

## Sub-Satellite Longitude Crossing

`longitude_condition` is zero when the sub-satellite longitude crosses a target. Uses `sin(lon - lon_target)` to handle the ±π discontinuity. Longitude is computed via approximate GMST rotation from ECI.

!!! note "Requires `integrator.p.JD`"
    The GMST computation needs the epoch Julian Date from `integrator.p.JD`.

### Detect passage over ground target

```julia
# Passage over Washington, DC (77°W)
target_lon = deg2rad(-77.0)
cond = longitude_condition(Cartesian, target_lon)
cb = ContinuousCallback(cond, affect!, affect!)
```

### Detect prime meridian crossing

```julia
cond = longitude_condition(Cartesian, 0.0)
cb = ContinuousCallback(cond, affect!, affect!)
```

### With regularized coordinates

```julia
cond = longitude_condition(EDromo, config, deg2rad(-77.0))
```

## Beta Angle Crossing

`beta_angle_condition` is zero when the beta angle (angle between orbit plane and Sun direction) crosses a target. Uses `asin(ĥ · ŝ) - β_target` where `ĥ` is the orbit normal and `ŝ` is the unit Sun direction.

Beta angle determines eclipse duration and thermal environment:
- `β = 0°` → Sun in orbit plane, maximum eclipse duration
- `β = 90°` → Sun normal to orbit plane, no eclipses possible

!!! note "Requires Sun ephemeris and `integrator.p.JD`"
    Beta angle computation needs a `ThirdBodyModel` for the Sun and the epoch Julian Date.

### Detect zero beta angle

```julia
using AstroForceModels

sun = ThirdBodyModel(; body=SunBody(), eop_data=eop_data)

cond = beta_angle_condition(Cartesian, sun, 0.0)
cb = ContinuousCallback(cond, affect!, affect!)
```

### Terminate when beta drops below critical value

```julia
# ISS eclipse-free boundary ≈ 14.5°
β_crit = deg2rad(14.5)
cond = beta_angle_condition(Cartesian, sun, β_crit)
cb = ContinuousCallback(cond,
    nothing,
    (integrator) -> terminate!(integrator),
)
```

### Dawn-dusk SSO transition detection

```julia
# Dawn-dusk sun-synchronous orbit transition ≈ 75.5°
cond = beta_angle_condition(Cartesian, sun, deg2rad(75.5))
cb = ContinuousCallback(cond, affect!, affect!)
```

## API Reference

```@docs
altitude_condition
latitude_condition
longitude_condition
beta_angle_condition
```
