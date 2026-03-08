export date_condition, negate_condition, and_condition, or_condition, shift_condition

# ========================================================================================
# Utility and composite event condition functions
# ========================================================================================

"""
    date_condition(::Type{C}, t_target) where {C<:AstroCoord}
    date_condition(::Type{C}, config, t_target) where {C<:AstroCoord}

Condition function zero when physical time equals `t_target` [s].
For regularized sets, extracts physical time from the state vector.
"""
function date_condition(::Type{C}, t_target::Number) where {C<:AstroCoords.AstroCoord}
    return (u, t, integrator) -> t - t_target
end

function date_condition(
    ::Type{C}, config::RegularizedCoordinateConfig, t_target::Number
) where {C<:AstroCoords.AstroCoord}
    return (u, t, integrator) -> get_physical_time(u, t, C, config) - t_target
end

"""
    negate_condition(condition)

Flip the sign of a condition function, swapping positive/negative crossing semantics.
"""
function negate_condition(condition)
    return (u, t, integrator) -> -condition(u, t, integrator)
end

"""
    and_condition(g1, g2)

Logical AND of two condition functions via `min(g1, g2)`. Negative only when
both are negative.
"""
function and_condition(g1, g2)
    return (u, t, integrator) -> min(g1(u, t, integrator), g2(u, t, integrator))
end

"""
    or_condition(g1, g2)

Logical OR of two condition functions via `max(g1, g2)`. Positive when either
is positive.
"""
function or_condition(g1, g2)
    return (u, t, integrator) -> max(g1(u, t, integrator), g2(u, t, integrator))
end

"""
    shift_condition(condition, Δt, ::Type{C}) where {C<:AstroCoord}
    shift_condition(condition, Δt, ::Type{C}, config) where {C<:AstroCoord}

Shift an event `Δt` seconds earlier (positive) or later (negative) by offsetting
the time argument. Approximate — shifts the time, not the exact event location.
"""
function shift_condition(condition, Δt::Number, ::Type{C}) where {C<:AstroCoords.AstroCoord}
    return (u, t, integrator) -> condition(u, t + Δt, integrator)
end

function shift_condition(
    condition, Δt::Number, ::Type{C}, config::RegularizedCoordinateConfig
) where {C<:AstroCoords.AstroCoord}
    return (u, t, integrator) -> condition(u, t + Δt, integrator)
end
