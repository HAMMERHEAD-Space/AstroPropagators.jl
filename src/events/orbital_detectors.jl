export apside_condition,
    node_condition,
    true_anomaly_condition,
    argument_of_latitude_condition,
    mean_anomaly_condition,
    raan_condition

# ========================================================================================
# Orbital event condition functions
#
# Each returns a closure `(u, t, integrator) -> Float64` suitable for use as the
# `condition` argument of a DifferentialEquations.jl `ContinuousCallback`.
# Zero-crossings of the returned value trigger the event.
#
# First argument is always the AstroCoord type that the integrator state vector
# represents. For regularized sets, pass RegularizedCoordinateConfig as the
# second argument (before any detector-specific parameters).
# ========================================================================================

"""
    apside_condition(::Type{C}) where {C<:AstroCoord}
    apside_condition(::Type{C}, config::RegularizedCoordinateConfig) where {C<:AstroCoord}

Condition function zero at apsides (periapsis and apoapsis). The g-function is
`dot(r, v)`. Use `affect!` / `affect_neg!` to distinguish periapsis from apoapsis.
"""
function apside_condition(::Type{C}) where {C<:AstroCoords.AstroCoord}
    return (u, t, integrator) -> begin
        cart = get_cartesian(u, t, integrator.p.μ, C)
        r = SVector{3}(cart[1], cart[2], cart[3])
        v = SVector{3}(cart[4], cart[5], cart[6])
        return dot(r, v)
    end
end

function apside_condition(
    ::Type{C}, config::RegularizedCoordinateConfig
) where {C<:AstroCoords.AstroCoord}
    return (u, t, integrator) -> begin
        cart = get_cartesian(u, t, integrator.p.μ, C, config)
        r = SVector{3}(cart[1], cart[2], cart[3])
        v = SVector{3}(cart[4], cart[5], cart[6])
        return dot(r, v)
    end
end

"""
    node_condition(::Type{C}) where {C<:AstroCoord}
    node_condition(::Type{C}, config::RegularizedCoordinateConfig) where {C<:AstroCoord}

Condition function zero at ascending/descending node crossings. The g-function is
the z-component of ECI position. Ascending node → `affect!`, descending → `affect_neg!`.
"""
function node_condition(::Type{C}) where {C<:AstroCoords.AstroCoord}
    return (u, t, integrator) -> begin
        cart = get_cartesian(u, t, integrator.p.μ, C)
        return cart[3]
    end
end

function node_condition(
    ::Type{C}, config::RegularizedCoordinateConfig
) where {C<:AstroCoords.AstroCoord}
    return (u, t, integrator) -> begin
        cart = get_cartesian(u, t, integrator.p.μ, C, config)
        return cart[3]
    end
end

"""
    true_anomaly_condition(::Type{C}, f_target::Number) where {C<:AstroCoord}
    true_anomaly_condition(::Type{C}, config, f_target) where {C<:AstroCoord}

Condition function zero when true anomaly equals `f_target` [rad].
Uses `sin(f - f_target)` to handle the 2π discontinuity.
"""
function true_anomaly_condition(
    ::Type{C}, f_target::Number
) where {C<:AstroCoords.AstroCoord}
    return (u, t, integrator) -> begin
        kep = get_keplerian(u, t, integrator.p.μ, C)
        return sin(kep.f - f_target)
    end
end

function true_anomaly_condition(
    ::Type{C}, config::RegularizedCoordinateConfig, f_target::Number
) where {C<:AstroCoords.AstroCoord}
    return (u, t, integrator) -> begin
        kep = get_keplerian(u, t, integrator.p.μ, C, config)
        return sin(kep.f - f_target)
    end
end

"""
    argument_of_latitude_condition(::Type{C}, u_target::Number) where {C<:AstroCoord}
    argument_of_latitude_condition(::Type{C}, config, u_target) where {C<:AstroCoord}

Condition function zero when argument of latitude `ω + f` equals `u_target` [rad].
Uses `sin(ω + f - u_target)` to handle the 2π discontinuity.
"""
function argument_of_latitude_condition(
    ::Type{C}, u_target::Number
) where {C<:AstroCoords.AstroCoord}
    return (u, t, integrator) -> begin
        kep = get_keplerian(u, t, integrator.p.μ, C)
        return sin(kep.ω + kep.f - u_target)
    end
end

function argument_of_latitude_condition(
    ::Type{C}, config::RegularizedCoordinateConfig, u_target::Number
) where {C<:AstroCoords.AstroCoord}
    return (u, t, integrator) -> begin
        kep = get_keplerian(u, t, integrator.p.μ, C, config)
        return sin(kep.ω + kep.f - u_target)
    end
end

"""
    mean_anomaly_condition(::Type{C}, M_target::Number) where {C<:AstroCoord}
    mean_anomaly_condition(::Type{C}, config, M_target) where {C<:AstroCoord}

Condition function zero when mean anomaly equals `M_target` [rad].
Uses `sin(M - M_target)` to handle the 2π discontinuity.
"""
function mean_anomaly_condition(
    ::Type{C}, M_target::Number
) where {C<:AstroCoords.AstroCoord}
    return (u, t, integrator) -> begin
        kep = get_keplerian(u, t, integrator.p.μ, C)
        return sin(kep.M - M_target)
    end
end

function mean_anomaly_condition(
    ::Type{C}, config::RegularizedCoordinateConfig, M_target::Number
) where {C<:AstroCoords.AstroCoord}
    return (u, t, integrator) -> begin
        kep = get_keplerian(u, t, integrator.p.μ, C, config)
        return sin(kep.M - M_target)
    end
end

"""
    raan_condition(::Type{C}, Ω_target::Number) where {C<:AstroCoord}
    raan_condition(::Type{C}, config, Ω_target) where {C<:AstroCoord}

Condition function zero when RAAN equals `Ω_target` [rad].
Uses `sin(Ω - Ω_target)` to handle the 2π discontinuity.
"""
function raan_condition(::Type{C}, Ω_target::Number) where {C<:AstroCoords.AstroCoord}
    return (u, t, integrator) -> begin
        kep = get_keplerian(u, t, integrator.p.μ, C)
        return sin(kep.Ω - Ω_target)
    end
end

function raan_condition(
    ::Type{C}, config::RegularizedCoordinateConfig, Ω_target::Number
) where {C<:AstroCoords.AstroCoord}
    return (u, t, integrator) -> begin
        kep = get_keplerian(u, t, integrator.p.μ, C, config)
        return sin(kep.Ω - Ω_target)
    end
end
