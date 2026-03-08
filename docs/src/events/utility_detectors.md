# Utility Detectors

Utility detectors provide time-based events and tools for composing or modifying other event conditions.

## Date/Time Crossing

`date_condition` is zero when physical time equals a target. For standard coordinate sets, this is simply `t - t_target`. For regularized sets, physical time is extracted from the state vector.

### Terminate at specific time

```julia
cond = date_condition(Cartesian, 43200.0)  # 12 hours
cb = ContinuousCallback(cond, (integrator) -> terminate!(integrator))
```

### With EDromo (regularized)

The condition automatically handles fictitious-to-physical time conversion:

```julia
cond = date_condition(EDromo, config, 43200.0)
cb = ContinuousCallback(cond, (integrator) -> terminate!(integrator))
```

!!! tip
    For simple integration termination, `build_termination_callback` is a shortcut:
    ```julia
    cb = build_termination_callback(43200.0, Cartesian)
    cb = build_termination_callback(43200.0, EDromo, config)
    ```

## Boolean Composition

Combine condition functions using logical AND/OR operations. These use `min`/`max` approximations, which work well for combining conditions with clear sign separation.

### AND — both conditions must be satisfied

`and_condition(g1, g2)` returns `min(g1, g2)`, which is negative only when both conditions are negative.

```julia
# In eclipse AND above 45° latitude
eclipse_cond = eclipse_condition(Cartesian, Conical(), sun_model)
lat_cond     = latitude_condition(Cartesian, deg2rad(45))

combined = and_condition(eclipse_cond, lat_cond)
cb = ContinuousCallback(combined, affect!, affect!)
```

### OR — either condition is sufficient

`or_condition(g1, g2)` returns `max(g1, g2)`, which is positive when either condition is positive.

```julia
# Periapsis OR descending node
apside_cond = apside_condition(Cartesian)
node_cond   = node_condition(Cartesian)

combined = or_condition(apside_cond, node_cond)
cb = ContinuousCallback(combined, affect!, affect!)
```

!!! warning
    `min`/`max` composition can reduce root-finding precision near corners where both conditions are simultaneously near zero. For critical applications, consider using separate callbacks via `CallbackSet` instead.

## Negation

`negate_condition(g)` flips the sign, swapping the positive/negative crossing semantics (`affect!` ↔ `affect_neg!`).

```julia
# Original: affect! fires on ascending node
cond = node_condition(Cartesian)

# Negated: affect! fires on descending node
neg_cond = negate_condition(cond)
```

## Time Shifting

`shift_condition` offsets an event by `Δt` seconds. Positive `Δt` fires *earlier* than the original; negative fires *later*. This is an approximation — it shifts the time argument, not the exact event location.

### Fire 5 minutes before eclipse entry

```julia
eclipse_cond = eclipse_condition(Cartesian, Conical(), sun_model)

# 300 seconds earlier
early_warning = shift_condition(eclipse_cond, 300.0, Cartesian)
cb = ContinuousCallback(early_warning, affect!, affect!)
```

### Delayed event

```julia
# Fire 60 seconds after apoapsis
apside_cond = apside_condition(Cartesian)
delayed = shift_condition(apside_cond, -60.0, Cartesian)
```

!!! note
    For slowly-varying conditions this approximation is very accurate. For fast-varying conditions (e.g., near periapsis on a highly eccentric orbit), the shifted event may not precisely correspond to `Δt` before/after the real event.

## Chaining with Other Detectors

Utility detectors are designed to compose with orbital, eclipse, and geometric detectors:

```julia
# Fire a maneuver 5 minutes before the 3rd periapsis
periapsis_cond = apside_condition(Cartesian)
early_cond = shift_condition(periapsis_cond, 300.0, Cartesian)

trigger = EventTrigger(early_cond, 3)  # 3rd occurrence
effect  = FixedDeltaV([0.1, 0.0, 0.0], RTNFrame())
cb = build_maneuver_callback(trigger, effect, Cartesian)
```

## API Reference

```@docs
date_condition
negate_condition
and_condition
or_condition
shift_condition
```
