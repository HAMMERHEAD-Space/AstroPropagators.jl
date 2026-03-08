# Access / Visibility Events

Access and visibility detectors are provided as a **package extension** that loads automatically when `AstroMeasurements.jl` is available. These detectors use `AstroMeasurements` for the underlying geometry computations.

## Setup

```julia
using AstroPropagators
using AstroMeasurements  # triggers the extension
```

Once both packages are loaded, `elevation_condition`, `relative_distance_condition`, and `angular_separation_condition` become available.

## Elevation Above Ground Station

`elevation_condition` is zero when the satellite's elevation above a ground station crosses a minimum threshold.

### Detect visibility windows

```julia
using AstroPropagators, AstroMeasurements, AstroCoords, AstroForceModels
using SatelliteToolboxTransformations, ComponentArrays

eop_data = fetch_iers_eop()

# Define ground station
gs = GroundStation(; lat=deg2rad(38.9), lon=deg2rad(-77.0), alt=0.0)

# Minimum elevation 10° (default)
cond = elevation_condition(Cartesian, gs, eop_data)

passes = []
cb = ContinuousCallback(
    cond,
    (integrator) -> push!(passes, (t=integrator.t, event=:aos)),   # AOS
    (integrator) -> push!(passes, (t=integrator.t, event=:los)),   # LOS
)

sol = propagate(CowellPropagator(), u0_cart, p, models, tspan; callback=cb)
```

### Custom minimum elevation

```julia
# 5° minimum elevation
cond = elevation_condition(Cartesian, gs, eop_data; el_min=deg2rad(5.0))

# 20° minimum elevation (high-quality passes only)
cond = elevation_condition(Cartesian, gs, eop_data; el_min=deg2rad(20.0))
```

### Terminate at AOS

```julia
cond = elevation_condition(Cartesian, gs, eop_data)
cb = ContinuousCallback(cond,
    (integrator) -> terminate!(integrator);  # terminate at AOS
    affect_neg! = nothing,
)
```

### With regularized coordinates

```julia
cond = elevation_condition(EDromo, config, gs, eop_data; el_min=deg2rad(10.0))
```

## Relative Distance

`relative_distance_condition` is zero when the distance from the spacecraft to a target crosses a threshold.

The target position is provided as a function `target_fn(JD) -> [x, y, z]` returning ECI position [km].

### Detect close approach to another object

```julia
# Target position function (e.g., from ephemeris or another propagation)
function target_pos(JD)
    # return ECI position [km] at Julian Date JD
    return [x, y, z]
end

# Detect when distance crosses 100 km
cond = relative_distance_condition(Cartesian, target_pos, 100.0)

approaches = []
cb = ContinuousCallback(
    cond,
    (integrator) -> push!(approaches, (t=integrator.t, dir=:departing)),
    (integrator) -> push!(approaches, (t=integrator.t, dir=:approaching)),
)
```

### Terminate at minimum safe distance

```julia
cond = relative_distance_condition(Cartesian, target_pos, 10.0)  # 10 km
cb = ContinuousCallback(cond,
    nothing,
    (integrator) -> terminate!(integrator),  # stop if getting closer
)
```

## Angular Separation

`angular_separation_condition` is zero when the angular separation between the spacecraft and a beacon, as seen from an observer, crosses a threshold.

Both beacon and observer positions are provided as functions `fn(JD) -> [x, y, z]` returning ECI positions [km].

### Sun exclusion zone

Detect when the spacecraft enters/exits an angular exclusion zone around the Sun as seen from a ground station:

```julia
# Sun position function
sun_pos(JD) = sun_model(JD, Position()) ./ 1E3

# Ground station position function
gs_pos(JD) = get_observer_state(gs, JD, eop_data)[1:3]

# 15° exclusion zone
cond = angular_separation_condition(Cartesian, sun_pos, gs_pos, deg2rad(15.0))

cb = ContinuousCallback(
    cond,
    (integrator) -> push!(log, (t=integrator.t, event=:exit_zone)),
    (integrator) -> push!(log, (t=integrator.t, event=:enter_zone)),
)
```

### With regularized coordinates

```julia
cond = angular_separation_condition(EDromo, config, sun_pos, gs_pos, deg2rad(15.0))
```

## Multi-Station Visibility

Combine multiple ground station elevation conditions to detect visibility from any station:

```julia
gs1 = GroundStation(; lat=deg2rad(38.9), lon=deg2rad(-77.0), alt=0.0)  # DC
gs2 = GroundStation(; lat=deg2rad(51.5), lon=deg2rad(-0.1), alt=0.0)   # London
gs3 = GroundStation(; lat=deg2rad(35.7), lon=deg2rad(139.7), alt=0.0)  # Tokyo

cond1 = elevation_condition(Cartesian, gs1, eop_data)
cond2 = elevation_condition(Cartesian, gs2, eop_data)
cond3 = elevation_condition(Cartesian, gs3, eop_data)

# Detect visibility from any station using or_condition
any_visible = or_condition(or_condition(cond1, cond2), cond3)
cb = ContinuousCallback(any_visible, aos_affect!, los_affect!)
```

## API Reference

```@docs
elevation_condition
relative_distance_condition
angular_separation_condition
```
