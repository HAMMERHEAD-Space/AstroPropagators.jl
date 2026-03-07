module AstroPropagatorsMeasurementsExt

using AstroPropagators
using AstroMeasurements
using AstroCoords
using AstroForceModels: Position
using StaticArraysCore
using LinearAlgebra

import AstroPropagators: get_cartesian, get_physical_time

"""
    elevation_condition(
        ::Type{C},
        ground_station::AbstractGroundStation,
        eop_data;
        el_min=deg2rad(10.0),
    ) where {C<:AstroCoord}
    elevation_condition(::Type{C}, config, ground_station, eop_data; el_min=...)

Return a condition function that is zero when the satellite elevation above
a ground station crosses `el_min` [rad].

Uses `AstroMeasurements.compute_elevation` for the geometry.
"""
function AstroPropagators.elevation_condition(
    ::Type{C},
    ground_station::AstroMeasurements.AbstractGroundStation,
    eop_data;
    el_min::Number=deg2rad(10.0),
) where {C<:AstroCoords.AstroCoord}
    return (u, t, integrator) -> begin
        cart = get_cartesian(u, t, integrator.p.μ, C)
        JD = integrator.p.JD + t / 86400.0
        el = AstroMeasurements.compute_elevation(
            SVector{6}(cart[1], cart[2], cart[3], cart[4], cart[5], cart[6]),
            ground_station,
            JD,
            eop_data,
        )
        return el - el_min
    end
end

function AstroPropagators.elevation_condition(
    ::Type{C},
    config::AstroCoords.RegularizedCoordinateConfig,
    ground_station::AstroMeasurements.AbstractGroundStation,
    eop_data;
    el_min::Number=deg2rad(10.0),
) where {C<:AstroCoords.AstroCoord}
    return (u, t, integrator) -> begin
        cart = get_cartesian(u, t, integrator.p.μ, C, config)
        t_phys = get_physical_time(u, t, C, config)
        JD = integrator.p.JD + t_phys / 86400.0
        el = AstroMeasurements.compute_elevation(
            SVector{6}(cart[1], cart[2], cart[3], cart[4], cart[5], cart[6]),
            ground_station,
            JD,
            eop_data,
        )
        return el - el_min
    end
end

"""
    relative_distance_condition(
        ::Type{C},
        target_state_fn,
        d_target::Number,
    ) where {C<:AstroCoord}
    relative_distance_condition(::Type{C}, config, target_state_fn, d_target)

Return a condition function that is zero when the distance from the spacecraft
to a target crosses `d_target` [km].

`target_state_fn(JD)` must return the target's ECI position as a 3-element vector [km].
"""
function AstroPropagators.relative_distance_condition(
    ::Type{C}, target_state_fn, d_target::Number
) where {C<:AstroCoords.AstroCoord}
    return (u, t, integrator) -> begin
        cart = get_cartesian(u, t, integrator.p.μ, C)
        r_sc = SVector{3}(cart[1], cart[2], cart[3])
        JD = integrator.p.JD + t / 86400.0
        r_target = SVector{3}(target_state_fn(JD)...)
        return norm(r_sc - r_target) - d_target
    end
end

function AstroPropagators.relative_distance_condition(
    ::Type{C},
    config::AstroCoords.RegularizedCoordinateConfig,
    target_state_fn,
    d_target::Number,
) where {C<:AstroCoords.AstroCoord}
    return (u, t, integrator) -> begin
        cart = get_cartesian(u, t, integrator.p.μ, C, config)
        t_phys = get_physical_time(u, t, C, config)
        JD = integrator.p.JD + t_phys / 86400.0
        r_sc = SVector{3}(cart[1], cart[2], cart[3])
        r_target = SVector{3}(target_state_fn(JD)...)
        return norm(r_sc - r_target) - d_target
    end
end

"""
    angular_separation_condition(
        ::Type{C},
        beacon_fn,
        observer_fn,
        threshold::Number,
    ) where {C<:AstroCoord}
    angular_separation_condition(::Type{C}, config, beacon_fn, observer_fn, threshold)

Return a condition function that is zero when the angular separation between the
spacecraft and a beacon, as seen from an observer, crosses `threshold` [rad].

- `beacon_fn(JD)` returns the beacon's ECI position [km] (e.g., Sun).
- `observer_fn(JD)` returns the observer's ECI position [km] (e.g., ground station).
"""
function AstroPropagators.angular_separation_condition(
    ::Type{C}, beacon_fn, observer_fn, threshold::Number
) where {C<:AstroCoords.AstroCoord}
    return (u, t, integrator) -> begin
        cart = get_cartesian(u, t, integrator.p.μ, C)
        r_sc = SVector{3}(cart[1], cart[2], cart[3])
        JD = integrator.p.JD + t / 86400.0
        r_beacon = SVector{3}(beacon_fn(JD)...)
        r_observer = SVector{3}(observer_fn(JD)...)
        d_sc = r_sc - r_observer
        d_beacon = r_beacon - r_observer
        cosθ = dot(d_sc, d_beacon) / (norm(d_sc) * norm(d_beacon))
        cosθ = clamp(cosθ, -one(cosθ), one(cosθ))
        return acos(cosθ) - threshold
    end
end

function AstroPropagators.angular_separation_condition(
    ::Type{C},
    config::AstroCoords.RegularizedCoordinateConfig,
    beacon_fn,
    observer_fn,
    threshold::Number,
) where {C<:AstroCoords.AstroCoord}
    return (u, t, integrator) -> begin
        cart = get_cartesian(u, t, integrator.p.μ, C, config)
        t_phys = get_physical_time(u, t, C, config)
        JD = integrator.p.JD + t_phys / 86400.0
        r_sc = SVector{3}(cart[1], cart[2], cart[3])
        r_beacon = SVector{3}(beacon_fn(JD)...)
        r_observer = SVector{3}(observer_fn(JD)...)
        d_sc = r_sc - r_observer
        d_beacon = r_beacon - r_observer
        cosθ = dot(d_sc, d_beacon) / (norm(d_sc) * norm(d_beacon))
        cosθ = clamp(cosθ, -one(cosθ), one(cosθ))
        return acos(cosθ) - threshold
    end
end

end # module
