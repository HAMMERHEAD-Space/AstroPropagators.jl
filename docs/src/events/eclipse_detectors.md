# Eclipse Detectors

Eclipse event detection reuses the shadow model infrastructure from AstroForceModels.jl. This guarantees consistency between the SRP force model and event detection — the same shadow geometry, Sun ephemeris, and body radii are used.

## Shadow Factor

`AstroForceModels.shadow_model` returns a continuous shadow factor:
- `0.0` = full shadow (umbra)
- `1.0` = full sunlight
- Between 0 and 1 = partial shadow (penumbra)

The condition function is `shadow_factor - threshold`, so:
- **Eclipse entry** (sunlight → shadow): negative-going crossing (`affect_neg!`)
- **Eclipse exit** (shadow → sunlight): positive-going crossing (`affect!`)

## Threshold Selection

The `threshold` keyword controls which transition is detected:

| Threshold | Detects |
|---|---|
| `≈ 0.5` (default) | General shadow entry/exit |
| `≈ 0.01` | Umbra entry/exit |
| `≈ 0.99` | Penumbra entry/exit |

## Shadow Model Types

AstroForceModels provides three shadow geometries:

| Type | Description |
|---|---|
| `Cylindrical()` | Simplest model, no penumbra |
| `Conical()` | Physically accurate, partial shadow region |
| `SmoothedConical()` | Smoothed version for better numerical behavior |

## Usage Examples

### Detect all eclipse transitions

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

sun_model = ThirdBodyModel(; body=SunBody(), eop_data=eop_data)

cond = eclipse_condition(Cartesian, Conical(), sun_model)

eclipse_log = []
cb = ContinuousCallback(
    cond,
    (integrator) -> push!(eclipse_log, (t=integrator.t, event=:exit)),
    (integrator) -> push!(eclipse_log, (t=integrator.t, event=:entry)),
)

sol = propagate(CowellPropagator(), u0_cart, p, models, (0.0, 86400.0); callback=cb)
```

### Terminate at umbra entry

```julia
cond = eclipse_condition(Cartesian, Conical(), sun_model; threshold=0.01)
cb = ContinuousCallback(cond,
    nothing,  # ignore exit
    (integrator) -> terminate!(integrator),  # terminate on entry
)
```

### Reuse existing SRP model

If you already have an `SRPAstroModel` configured, the convenience overload ensures the eclipse detector uses the exact same shadow configuration:

```julia
srp = SRPAstroModel(;
    satellite_srp_model=CannonballFixedSRP(0.2),
    sun_data=sun_model, eop_data=eop_data, shadow_model=SmoothedConical(),
)

cond = eclipse_condition(Cartesian, srp; threshold=0.5)
```

### With regularized coordinates (EDromo)

```julia
cond = eclipse_condition(EDromo, config, Conical(), sun_model; threshold=0.5)

end_cb = build_termination_callback(86400.0, EDromo, config)
cb = CallbackSet(
    ContinuousCallback(cond, exit_affect!, entry_affect!),
    end_cb,
)
```

### Eclipse duration computation

Combine entry/exit logging to compute eclipse durations:

```julia
cond = eclipse_condition(Cartesian, Conical(), sun_model)

transitions = []
cb = ContinuousCallback(
    cond,
    (integrator) -> push!(transitions, (t=integrator.t, event=:exit)),
    (integrator) -> push!(transitions, (t=integrator.t, event=:entry)),
)

sol = propagate(CowellPropagator(), u0_cart, p, models, (0.0, 86400.0); callback=cb)

# Compute durations from entry/exit pairs
for i in 1:2:length(transitions)-1
    if transitions[i].event == :entry && transitions[i+1].event == :exit
        duration = transitions[i+1].t - transitions[i].t
        println("Eclipse duration: $(duration / 60.0) minutes")
    end
end
```

## With ODEProblem Directly

Eclipse conditions work identically when bypassing `propagate()`:

```julia
f(u, p, t) = Cowell_EOM(u, p, t, models)
prob = ODEProblem{false}(f, u0_cart, (0.0, 86400.0), p)

cond = eclipse_condition(Cartesian, Conical(), sun_model)
cb = ContinuousCallback(cond, exit_affect!, entry_affect!)

sol = solve(prob, Vern9(); abstol=1e-13, reltol=1e-13, callback=cb)
```

## API Reference

```@docs
eclipse_condition
```
