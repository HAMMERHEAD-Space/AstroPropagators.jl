# Maneuvers

The maneuver system decouples *when* a maneuver fires (triggers) from *what* happens (effects). This enables flexible composition of burn timing strategies with delta-V computation methods.

## Architecture

```
┌────────────────┐     ┌──────────────────┐
│    Triggers     │     │     Effects       │
│                 │     │                   │
│  TimeTrigger    │     │  FixedDeltaV      │
│  EventTrigger   │     │  ComputedDeltaV   │
└────────┬───────┘     └────────┬──────────┘
         │                      │
         └──────────┬───────────┘
                    │
           build_maneuver_callback
                    │
                    ▼
            ContinuousCallback
```

Alternatively, event conditions can be paired with generic actions (not just maneuvers):

```
┌──────────────────┐     ┌──────────────────┐
│  Any condition    │     │   Event Actions   │
│  (from detectors) │     │                   │
│                   │     │  TerminateAction   │
│                   │     │  LogAction         │
│                   │     │  ContinueAction    │
│                   │     │  ManeuverAction     │
└────────┬──────────┘     └────────┬──────────┘
         │                         │
         └──────────┬──────────────┘
                    │
            build_event_callback
                    │
                    ▼
             ContinuousCallback
```

## Triggers

### TimeTrigger

Fire at a specific physical time:

```julia
trigger = TimeTrigger(43200.0)  # 12 hours into propagation
```

### EventTrigger

Fire when a condition function crosses zero:

```julia
# Fire at first periapsis
trigger = EventTrigger(apside_condition(Cartesian))

# Fire at 3rd periapsis
trigger = EventTrigger(apside_condition(Cartesian), 3)

# Fire at ascending node
trigger = EventTrigger(node_condition(Cartesian))

# Fire at specific true anomaly
trigger = EventTrigger(true_anomaly_condition(Cartesian, π))
```

## Effects

### FixedDeltaV

Apply a predetermined delta-V vector:

```julia
# Inertial frame (default)
effect = FixedDeltaV([0.1, 0.0, 0.0])

# RTN frame (radial-transverse-normal)
effect = FixedDeltaV([0.1, 0.0, 0.0], RTNFrame())

# VNB frame (velocity-normal-binormal)
effect = FixedDeltaV([0.0, 0.0, 0.05], VNBFrame())
```

### ComputedDeltaV

Compute delta-V at event time from the current state. The compute function receives
the Cartesian state, physical time, and integrator parameters:

```julia
# Hohmann-like: compute tangential ΔV based on current state
function hohmann_dv(cart_state, t_phys, p)
    r = sqrt(cart_state[1]^2 + cart_state[2]^2 + cart_state[3]^2)
    v = sqrt(cart_state[4]^2 + cart_state[5]^2 + cart_state[6]^2)
    v_target = sqrt(p.μ * (2/r - 1/a_target))
    return SVector(0.0, v_target - v, 0.0)  # tangential only
end

effect = ComputedDeltaV(hohmann_dv, RTNFrame())
```

```julia
# State-dependent correction
function station_keeping_dv(cart_state, t_phys, p)
    # Compute required correction from current state
    return SVector(Δvx, Δvy, Δvz)
end

effect = ComputedDeltaV(station_keeping_dv, InertialFrame())
```

## Building Callbacks

### Single Maneuver: `build_maneuver_callback`

Combines a trigger and effect into a `ContinuousCallback`:

```julia
# Time-triggered burn in RTN frame
cb = build_maneuver_callback(
    TimeTrigger(43200.0),
    FixedDeltaV([0.1, 0.0, 0.0], RTNFrame()),
    Cartesian,
)

sol = propagate(CowellPropagator(), u0_cart, p, models, tspan; callback=cb)
```

```julia
# Burn at apoapsis (affect_neg! crossing of dot(r,v))
cb = build_maneuver_callback(
    EventTrigger(apside_condition(Cartesian)),
    FixedDeltaV([0.0, 0.5, 0.0], RTNFrame()),
    Cartesian,
)
```

```julia
# Computed ΔV at 3rd periapsis
cb = build_maneuver_callback(
    EventTrigger(apside_condition(Cartesian), 3),
    ComputedDeltaV(my_dv_func, RTNFrame()),
    Cartesian,
)
```

### Multiple Burns: `build_maneuver_schedule`

Combines a vector of `(trigger, effect)` pairs into a `CallbackSet`:

```julia
burns = [
    (TimeTrigger(3600.0),  FixedDeltaV([0.1, 0.0, 0.0], RTNFrame())),
    (TimeTrigger(7200.0),  FixedDeltaV([0.0, 0.05, 0.0], RTNFrame())),
    (TimeTrigger(10800.0), FixedDeltaV([-0.05, 0.0, 0.0], RTNFrame())),
]

cbs = build_maneuver_schedule(burns, Cartesian)
sol = propagate(CowellPropagator(), u0_cart, p, models, tspan; callback=cbs)
```

### Generic Events: `build_event_callback`

Pairs any condition function with an action:

```julia
# Terminate at first eclipse entry
cond = eclipse_condition(Cartesian, Conical(), sun_model)
cb = build_event_callback(cond, TerminateAction(), Cartesian)
```

```julia
# Log all periapsis crossings
log = []
cond = apside_condition(Cartesian)
cb = build_event_callback(cond, LogAction(log), Cartesian)
# log entries: [(t=..., state=Cartesian(...)), ...]
```

```julia
# Apply a maneuver at eclipse exit
cond = eclipse_condition(Cartesian, Conical(), sun_model)
action = ManeuverAction(FixedDeltaV([0.01, 0.0, 0.0], RTNFrame()))
cb = build_event_callback(cond, action, Cartesian)
```

### Termination: `build_termination_callback`

Shortcut for stopping integration at a physical time:

```julia
# Standard coordinates
cb = build_termination_callback(86400.0, Cartesian)

# Regularized coordinates (handles fictitious time automatically)
cb = build_termination_callback(86400.0, EDromo, config)
```

## Regularized Coordinates

All builders have overloads accepting a `RegularizedCoordinateConfig`. The config handles:
- Physical time extraction for `TimeTrigger`
- Coordinate conversion for applying burns
- Time element reset after burn application

```julia
config = RegularizedCoordinateConfig(u0_cart, μ; W=W, t₀=0.0, flag_time=PhysicalTime())

# Time-triggered burn with EDromo
cb = build_maneuver_callback(
    TimeTrigger(43200.0),
    FixedDeltaV([0.1, 0.0, 0.0], RTNFrame()),
    EDromo,
    config,
)

# Event-triggered burn with EDromo
cb = build_maneuver_callback(
    EventTrigger(apside_condition(EDromo, config)),
    FixedDeltaV([0.1, 0.0, 0.0], RTNFrame()),
    EDromo,
    config,
)

# Maneuver schedule with KS
burns = [
    (TimeTrigger(3600.0), FixedDeltaV([0.1, 0.0, 0.0])),
    (TimeTrigger(7200.0), FixedDeltaV([0.0, 0.05, 0.0])),
]
cbs = build_maneuver_schedule(burns, KustaanheimoStiefel, config)
```

## Thrust Frames

Delta-V vectors are specified in a thrust frame, then rotated to inertial for application via `AstroForceModels.transform_to_inertial`:

| Frame | Components | Description |
|---|---|---|
| `InertialFrame()` | `[Δvx, Δvy, Δvz]` | ECI inertial — no rotation |
| `RTNFrame()` | `[ΔvR, ΔvT, ΔvN]` | Radial, Transverse, Normal |
| `VNBFrame()` | `[ΔvV, ΔvN, ΔvB]` | Velocity, Normal, Binormal |

## Common Workflows

### Hohmann Transfer

```julia
# Burn 1: prograde at periapsis → raise apoapsis
cb1 = build_maneuver_callback(
    EventTrigger(apside_condition(Cartesian)),      # periapsis
    FixedDeltaV([0.0, Δv1, 0.0], RTNFrame()),      # prograde
    Cartesian,
)

# Burn 2: prograde at new apoapsis → circularize
cb2 = build_maneuver_callback(
    EventTrigger(apside_condition(Cartesian), 2),   # next apoapsis (2nd apside)
    FixedDeltaV([0.0, Δv2, 0.0], RTNFrame()),      # prograde
    Cartesian,
)

sol = propagate(CowellPropagator(), u0_cart, p, models, tspan;
    callback=CallbackSet(cb1, cb2),
)
```

### Plane Change at Node

```julia
cb = build_maneuver_callback(
    EventTrigger(node_condition(Cartesian)),         # ascending node
    FixedDeltaV([0.0, 0.0, Δv_normal], RTNFrame()), # normal component
    Cartesian,
)
```

### Station-Keeping with Computed Burns

```julia
function sk_correction(cart, t_phys, p)
    kep = Keplerian(Cartesian(cart), p.μ)
    Δa = a_target - kep.a
    Δv = Δa * p.μ / (2 * kep.a^2 * sqrt(p.μ / kep.a))
    return SVector(0.0, Δv, 0.0)
end

# Correct at every periapsis
cb = build_maneuver_callback(
    EventTrigger(apside_condition(Cartesian)),
    ComputedDeltaV(sk_correction, RTNFrame()),
    Cartesian,
)
```

## Low-Level: `impulsive_burn!`

For direct integrator manipulation (e.g., inside custom callbacks), `impulsive_burn!` applies a delta-V without the trigger/effect abstraction:

```julia
using AstroPropagators: impulsive_burn!

cb = ContinuousCallback(
    my_condition,
    (integrator) -> impulsive_burn!(
        integrator, [0.1, 0.0, 0.0];
        coordinate_set=Cartesian, frame=RTNFrame(),
    ),
)
```

## API Reference

```@docs
build_maneuver_callback
build_maneuver_schedule
build_event_callback
build_termination_callback
impulsive_burn!
```
