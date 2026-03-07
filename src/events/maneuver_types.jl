export AbstractManeuverTrigger,
    TimeTrigger,
    EventTrigger,
    AbstractManeuverEffect,
    FixedDeltaV,
    ComputedDeltaV,
    AbstractEventAction,
    TerminateAction,
    LogAction,
    ContinueAction,
    ManeuverAction

# ========================================================================================
# Maneuver Triggers — define WHEN a maneuver fires
# ========================================================================================

abstract type AbstractManeuverTrigger end

"""
    TimeTrigger(t_burn::Number)

Trigger a maneuver at physical time `t_burn` [s].
"""
struct TimeTrigger{T<:Number} <: AbstractManeuverTrigger
    t_burn::T
end

"""
    EventTrigger(condition; count=1)

Trigger a maneuver when `condition` crosses zero. `condition` must have signature
`(u, t, integrator) -> Float64`. Optional `count` fires on the Nth occurrence.
"""
struct EventTrigger{F,CT<:Integer} <: AbstractManeuverTrigger
    condition::F
    count::CT
end
EventTrigger(condition) = EventTrigger(condition, 1)

# ========================================================================================
# Maneuver Effects — define WHAT happens when a maneuver fires
# ========================================================================================

abstract type AbstractManeuverEffect end

"""
    FixedDeltaV(ΔV::AbstractVector; frame=InertialFrame())

Apply a fixed impulsive ΔV [km/s] in the specified thrust frame.
"""
struct FixedDeltaV{VT<:AbstractVector,FT<:AbstractThrustFrame} <: AbstractManeuverEffect
    ΔV::VT
    frame::FT
end
FixedDeltaV(ΔV::AbstractVector) = FixedDeltaV(ΔV, InertialFrame())

"""
    ComputedDeltaV(compute; frame=InertialFrame())

Apply an impulsive ΔV computed from current state at event time.
`compute(cart_state, t_physical, p) -> SVector{3}` returns ΔV [km/s]
in the specified thrust frame.
"""
struct ComputedDeltaV{F,FT<:AbstractThrustFrame} <: AbstractManeuverEffect
    compute::F
    frame::FT
end
ComputedDeltaV(compute) = ComputedDeltaV(compute, InertialFrame())

# ========================================================================================
# Event Actions — define how a generic event (non-maneuver) is handled
# ========================================================================================

abstract type AbstractEventAction end

"""
    TerminateAction()

Terminate integration when the event fires.
"""
struct TerminateAction <: AbstractEventAction end

"""
    LogAction(log::Vector)

Record event time and Cartesian state into `log`, then continue integration.
"""
struct LogAction{LT<:AbstractVector} <: AbstractEventAction
    log::LT
end

"""
    ContinueAction()

Continue integration without side effects. Default/placeholder action.
"""
struct ContinueAction <: AbstractEventAction end

"""
    ManeuverAction(effect::AbstractManeuverEffect)

Apply a maneuver effect when the event fires.
"""
struct ManeuverAction{ET<:AbstractManeuverEffect} <: AbstractEventAction
    effect::ET
end
