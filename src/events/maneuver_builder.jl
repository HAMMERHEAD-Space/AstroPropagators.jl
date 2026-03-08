export build_maneuver_callback,
    build_maneuver_schedule, build_event_callback, build_termination_callback

# ========================================================================================
# Condition builders — convert a trigger into a condition closure
# ========================================================================================

function _build_condition(trigger::TimeTrigger, ::Type{C}) where {C<:AstroCoords.AstroCoord}
    t_burn = trigger.t_burn
    return (u, t, integrator) -> t - t_burn
end

function _build_condition(
    trigger::TimeTrigger, ::Type{C}, config::RegularizedCoordinateConfig
) where {C<:AstroCoords.AstroCoord}
    t_burn = trigger.t_burn
    return (u, t, integrator) -> get_physical_time(u, t, C, config) - t_burn
end

function _build_condition(
    trigger::EventTrigger, ::Type{C}, args...
) where {C<:AstroCoords.AstroCoord}
    return trigger.condition
end

# ========================================================================================
# Affect builders — convert an effect into an affect! closure
# ========================================================================================

"""
    _apply_burn!(integrator, ΔV_inertial, ::Type{C}) where {C<:AstroCoord}
    _apply_burn!(integrator, ΔV_inertial, ::Type{C}, config) where {C<:AstroCoord}

Low-level: apply an inertial ΔV to the integrator state, handling coordinate
conversion and time-element updates for regularized formulations.
"""
function _apply_burn!(integrator, ΔV_inertial::SVector{3}, ::Type{Cartesian})
    u = integrator.u
    integrator.u = SVector{length(u)}(
        u[1],
        u[2],
        u[3],
        u[4] + ΔV_inertial[1],
        u[5] + ΔV_inertial[2],
        u[6] + ΔV_inertial[3],
    )
    return nothing
end

function _apply_burn!(
    integrator, ΔV_inertial::SVector{3}, ::Type{C}
) where {C<:AstroCoords.AstroCoord}
    cart_state = get_cartesian(integrator.u, integrator.t, integrator.p.μ, C)
    new_cart = Cartesian(
        cart_state[1],
        cart_state[2],
        cart_state[3],
        cart_state[4] + ΔV_inertial[1],
        cart_state[5] + ΔV_inertial[2],
        cart_state[6] + ΔV_inertial[3],
    )
    integrator.u = params(C(new_cart, integrator.p.μ))
    return nothing
end

function _apply_burn!(
    integrator, ΔV_inertial::SVector{3}, ::Type{EDromo}, config::RegularizedCoordinateConfig
)
    cart_state = get_cartesian(integrator.u, integrator.t, integrator.p.μ, EDromo, config)
    new_cart = Cartesian(
        cart_state[1],
        cart_state[2],
        cart_state[3],
        cart_state[4] + ΔV_inertial[1],
        cart_state[5] + ΔV_inertial[2],
        cart_state[6] + ΔV_inertial[3],
    )
    t_maneuver = get_EDromo_time(integrator.u, integrator.t, config)
    maneuver_config = RegularizedCoordinateConfig(
        config.DU, config.TU, config.W, t_maneuver, config.flag_time
    )
    integrator.u = params(EDromo(new_cart, integrator.p.μ, integrator.t, maneuver_config))
    return nothing
end

function _apply_burn!(
    integrator,
    ΔV_inertial::SVector{3},
    ::Type{KustaanheimoStiefel},
    config::RegularizedCoordinateConfig,
)
    cart_state = get_cartesian(
        integrator.u, integrator.t, integrator.p.μ, KustaanheimoStiefel, config
    )
    new_cart = Cartesian(
        cart_state[1],
        cart_state[2],
        cart_state[3],
        cart_state[4] + ΔV_inertial[1],
        cart_state[5] + ΔV_inertial[2],
        cart_state[6] + ΔV_inertial[3],
    )
    t_maneuver = get_KS_time(integrator.u, config)
    maneuver_config = RegularizedCoordinateConfig(
        config.DU, config.TU, config.W, t_maneuver, config.flag_time
    )
    integrator.u = params(KustaanheimoStiefel(new_cart, integrator.p.μ, maneuver_config))
    return nothing
end

function _apply_burn!(
    integrator,
    ΔV_inertial::SVector{3},
    ::Type{StiefelScheifele},
    config::RegularizedCoordinateConfig,
)
    cart_state = get_cartesian(
        integrator.u, integrator.t, integrator.p.μ, StiefelScheifele, config
    )
    new_cart = Cartesian(
        cart_state[1],
        cart_state[2],
        cart_state[3],
        cart_state[4] + ΔV_inertial[1],
        cart_state[5] + ΔV_inertial[2],
        cart_state[6] + ΔV_inertial[3],
    )
    t_maneuver = get_stiefelscheifele_time(integrator.u, integrator.t, config)
    maneuver_config = RegularizedCoordinateConfig(
        config.DU, config.TU, config.W, t_maneuver, config.flag_time
    )
    integrator.u = params(
        StiefelScheifele(new_cart, integrator.p.μ, integrator.t, maneuver_config)
    )
    return nothing
end

function _apply_burn!(
    integrator, ΔV_inertial::SVector{3}, ::Type{GEqOE}, config::RegularizedCoordinateConfig
)
    cart_state = get_cartesian(integrator.u, integrator.t, integrator.p.μ, GEqOE, config)
    new_cart = Cartesian(
        cart_state[1],
        cart_state[2],
        cart_state[3],
        cart_state[4] + ΔV_inertial[1],
        cart_state[5] + ΔV_inertial[2],
        cart_state[6] + ΔV_inertial[3],
    )
    integrator.u = params(GEqOE(new_cart, integrator.p.μ, config))
    return nothing
end

# --- Standard coord affect builders ---

function _build_affect(effect::FixedDeltaV, ::Type{C}) where {C<:AstroCoords.AstroCoord}
    ΔV = SVector{3}(effect.ΔV[1], effect.ΔV[2], effect.ΔV[3])
    frame = effect.frame
    return (integrator) -> begin
        cart = get_cartesian(integrator.u, integrator.t, integrator.p.μ, C)
        ΔV_inertial = transform_to_inertial(ΔV, cart, frame)
        _apply_burn!(integrator, ΔV_inertial, C)
    end
end

function _build_affect(
    effect::FixedDeltaV, ::Type{C}, config::RegularizedCoordinateConfig
) where {C<:AstroCoords.AstroCoord}
    ΔV = SVector{3}(effect.ΔV[1], effect.ΔV[2], effect.ΔV[3])
    frame = effect.frame
    return (integrator) -> begin
        cart = get_cartesian(integrator.u, integrator.t, integrator.p.μ, C, config)
        ΔV_inertial = transform_to_inertial(ΔV, cart, frame)
        _apply_burn!(integrator, ΔV_inertial, C, config)
    end
end

function _build_affect(effect::ComputedDeltaV, ::Type{C}) where {C<:AstroCoords.AstroCoord}
    compute = effect.compute
    frame = effect.frame
    return (integrator) -> begin
        cart = get_cartesian(integrator.u, integrator.t, integrator.p.μ, C)
        t_phys = get_physical_time(integrator.u, integrator.t, C)
        ΔV = SVector{3}(compute(cart, t_phys, integrator.p)...)
        ΔV_inertial = transform_to_inertial(ΔV, cart, frame)
        _apply_burn!(integrator, ΔV_inertial, C)
    end
end

function _build_affect(
    effect::ComputedDeltaV, ::Type{C}, config::RegularizedCoordinateConfig
) where {C<:AstroCoords.AstroCoord}
    compute = effect.compute
    frame = effect.frame
    return (integrator) -> begin
        cart = get_cartesian(integrator.u, integrator.t, integrator.p.μ, C, config)
        t_phys = get_physical_time(integrator.u, integrator.t, C, config)
        ΔV = SVector{3}(compute(cart, t_phys, integrator.p)...)
        ΔV_inertial = transform_to_inertial(ΔV, cart, frame)
        _apply_burn!(integrator, ΔV_inertial, C, config)
    end
end

# ========================================================================================
# Nth-occurrence wrapper for EventTrigger
# ========================================================================================

function _wrap_condition_with_count(condition, count::Integer)
    count == 1 && return condition
    counter = Ref(0)
    return (u, t, integrator) -> begin
        g = condition(u, t, integrator)
        return g
    end
end

function _wrap_affect_with_count(affect!, count::Integer)
    count == 1 && return affect!
    counter = Ref(0)
    return (integrator) -> begin
        counter[] += 1
        if counter[] == count
            affect!(integrator)
        end
    end
end

# ========================================================================================
# Public API
# ========================================================================================

"""
    build_maneuver_callback(trigger, effect, ::Type{C}) where {C<:AstroCoord}
    build_maneuver_callback(trigger, effect, ::Type{C}, config) where {C<:AstroCoord}

Build a `ContinuousCallback` that fires a maneuver (trigger + effect) for the
given coordinate type. Handles coordinate conversions and time-element updates
automatically.
"""
function build_maneuver_callback(
    trigger::AbstractManeuverTrigger, effect::AbstractManeuverEffect, ::Type{C}
) where {C<:AstroCoords.AstroCoord}
    condition = _build_condition(trigger, C)
    affect! = _build_affect(effect, C)
    if trigger isa EventTrigger && trigger.count > 1
        affect! = _wrap_affect_with_count(affect!, trigger.count)
    end
    return ContinuousCallback(condition, affect!)
end

function build_maneuver_callback(
    trigger::AbstractManeuverTrigger,
    effect::AbstractManeuverEffect,
    ::Type{C},
    config::RegularizedCoordinateConfig,
) where {C<:AstroCoords.AstroCoord}
    condition = _build_condition(trigger, C, config)
    affect! = _build_affect(effect, C, config)
    if trigger isa EventTrigger && trigger.count > 1
        affect! = _wrap_affect_with_count(affect!, trigger.count)
    end
    return ContinuousCallback(condition, affect!)
end

"""
    build_maneuver_schedule(burns, ::Type{C}) where {C<:AstroCoord}
    build_maneuver_schedule(burns, ::Type{C}, config) where {C<:AstroCoord}

Build a `CallbackSet` from a vector of `(trigger, effect)` pairs.
"""
function build_maneuver_schedule(
    burns::AbstractVector, ::Type{C}
) where {C<:AstroCoords.AstroCoord}
    callbacks = [build_maneuver_callback(t, e, C) for (t, e) in burns]
    return CallbackSet(callbacks...)
end

function build_maneuver_schedule(
    burns::AbstractVector, ::Type{C}, config::RegularizedCoordinateConfig
) where {C<:AstroCoords.AstroCoord}
    callbacks = [build_maneuver_callback(t, e, C, config) for (t, e) in burns]
    return CallbackSet(callbacks...)
end

"""
    build_event_callback(condition, action, ::Type{C}, [config])

Build a `ContinuousCallback` pairing any event condition with an action
(`TerminateAction`, `LogAction`, `ContinueAction`, `ManeuverAction`).
"""
function build_event_callback(
    condition, action::TerminateAction, ::Type{C}, args...
) where {C<:AstroCoords.AstroCoord}
    return ContinuousCallback(condition, (integrator) -> terminate!(integrator))
end

function build_event_callback(
    condition, action::ContinueAction, ::Type{C}, args...
) where {C<:AstroCoords.AstroCoord}
    return ContinuousCallback(condition, (integrator) -> nothing)
end

function build_event_callback(
    condition, action::LogAction, ::Type{C}
) where {C<:AstroCoords.AstroCoord}
    log = action.log
    return ContinuousCallback(
        condition,
        (integrator) -> begin
            cart = get_cartesian(integrator.u, integrator.t, integrator.p.μ, C)
            t_phys = get_physical_time(integrator.u, integrator.t, C)
            push!(log, (t=t_phys, state=cart))
        end,
    )
end

function build_event_callback(
    condition, action::LogAction, ::Type{C}, config::RegularizedCoordinateConfig
) where {C<:AstroCoords.AstroCoord}
    log = action.log
    return ContinuousCallback(
        condition,
        (integrator) -> begin
            cart = get_cartesian(integrator.u, integrator.t, integrator.p.μ, C, config)
            t_phys = get_physical_time(integrator.u, integrator.t, C, config)
            push!(log, (t=t_phys, state=cart))
        end,
    )
end

function build_event_callback(
    condition, action::ManeuverAction, ::Type{C}
) where {C<:AstroCoords.AstroCoord}
    affect! = _build_affect(action.effect, C)
    return ContinuousCallback(condition, affect!)
end

function build_event_callback(
    condition, action::ManeuverAction, ::Type{C}, config::RegularizedCoordinateConfig
) where {C<:AstroCoords.AstroCoord}
    affect! = _build_affect(action.effect, C, config)
    return ContinuousCallback(condition, affect!)
end

"""
    build_termination_callback(stop_time, ::Type{C}) where {C<:AstroCoord}
    build_termination_callback(stop_time, ::Type{C}, config) where {C<:AstroCoord}

Build a `ContinuousCallback` that terminates integration at `stop_time` [s].
"""
function build_termination_callback(
    stop_time::Number, ::Type{C}
) where {C<:AstroCoords.AstroCoord}
    return ContinuousCallback(
        (u, t, integrator) -> t - stop_time, (integrator) -> terminate!(integrator)
    )
end

function build_termination_callback(
    stop_time::Number, ::Type{C}, config::RegularizedCoordinateConfig
) where {C<:AstroCoords.AstroCoord}
    return ContinuousCallback(
        (u, t, integrator) -> get_physical_time(u, t, C, config) - stop_time,
        (integrator) -> terminate!(integrator),
    )
end
