# Orbital Detectors

Orbital event detectors trigger on orbital mechanics quantities: apsides, node crossings, and orbital element crossings. All are computed by converting the integrator state to `Keplerian` (or `Cartesian` for geometric conditions) via the state adapter layer.

## Apside Detection

`apside_condition` is zero when `dot(r, v) = 0`, which occurs at both periapsis and apoapsis.

The `ContinuousCallback` semantics distinguish the two:
- **Periapsis** (positive-going crossing): fires `affect!`
- **Apoapsis** (negative-going crossing): fires `affect_neg!`

### Log all apsides

```julia
cond = apside_condition(Cartesian)

peri_log = []
apo_log  = []

cb = ContinuousCallback(
    cond,
    (integrator) -> push!(peri_log, integrator.t),   # periapsis
    (integrator) -> push!(apo_log, integrator.t),     # apoapsis
)

sol = propagate(CowellPropagator(), u0_cart, p, models, tspan; callback=cb)
```

### Terminate at first periapsis

```julia
cond = apside_condition(Cartesian)
cb = ContinuousCallback(cond,
    (integrator) -> terminate!(integrator);
    affect_neg! = nothing,  # ignore apoapsis
)
```

### With EDromo (regularized)

```julia
cond = apside_condition(EDromo, config)
end_cb = build_termination_callback(86400.0, EDromo, config)

cb = CallbackSet(
    ContinuousCallback(cond, peri_affect!, apo_affect!),
    end_cb,
)
```

## Node Crossing Detection

`node_condition` is zero when the z-component of ECI position is zero (equatorial plane crossing).

- **Ascending node** (positive-going): fires `affect!`
- **Descending node** (negative-going): fires `affect_neg!`

### Log all node crossings

```julia
cond = node_condition(Cartesian)

asc_log  = []
desc_log = []

cb = ContinuousCallback(
    cond,
    (integrator) -> push!(asc_log, integrator.t),
    (integrator) -> push!(desc_log, integrator.t),
)
```

### Terminate at ascending node (Keplerian state)

```julia
cond = node_condition(Keplerian)  # works with GaussVE state vectors
cb = ContinuousCallback(cond,
    (integrator) -> terminate!(integrator);
    affect_neg! = nothing,
)
```

## True Anomaly Crossing

`true_anomaly_condition` is zero when true anomaly equals a target value. Uses `sin(f - f_target)` to handle the 2π discontinuity.

### Detect passage through f = π/2

```julia
cond = true_anomaly_condition(Cartesian, π/2)
cb = ContinuousCallback(cond, affect!, affect!)
```

### Trigger burn at specific true anomaly

```julia
trigger = EventTrigger(true_anomaly_condition(Cartesian, π))
effect  = FixedDeltaV([0.1, 0.0, 0.0], RTNFrame())
cb = build_maneuver_callback(trigger, effect, Cartesian)
```

### EDromo (config goes before f_target)

```julia
cond = true_anomaly_condition(EDromo, config, π/2)
```

## Argument of Latitude Crossing

`argument_of_latitude_condition` is zero when `ω + f` equals a target. Useful for near-circular orbits where periapsis is poorly defined.

```julia
# Detect passage through argument of latitude = 90°
cond = argument_of_latitude_condition(Cartesian, π/2)
cb = ContinuousCallback(cond, affect!, affect!)
```

## Mean Anomaly Crossing

`mean_anomaly_condition` is zero when mean anomaly equals a target. Mean anomaly progresses linearly with time for Keplerian motion, making this useful for uniformly-spaced event detection.

```julia
cond = mean_anomaly_condition(Cartesian, 0.0)  # periapsis passage in M
cb = ContinuousCallback(cond, affect!, affect!)
```

## RAAN Crossing

`raan_condition` is zero when the right ascension of the ascending node equals a target. RAAN drifts slowly due to J2, so crossings are infrequent for typical LEO orbits.

```julia
# Detect when RAAN crosses 0°
cond = raan_condition(Cartesian, 0.0)
cb = ContinuousCallback(cond, affect!, affect!)
```

## Combining Multiple Orbital Events

Multiple detectors compose naturally via `CallbackSet`:

```julia
peri_cond = apside_condition(Cartesian)
node_cond = node_condition(Cartesian)

cb = CallbackSet(
    ContinuousCallback(peri_cond, peri_affect!, apo_affect!),
    ContinuousCallback(node_cond, asc_affect!, desc_affect!),
)

sol = propagate(CowellPropagator(), u0_cart, p, models, tspan; callback=cb)
```

## API Reference

```@docs
apside_condition
node_condition
true_anomaly_condition
argument_of_latitude_condition
mean_anomaly_condition
raan_condition
```
