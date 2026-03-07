export altitude_condition, latitude_condition, longitude_condition, beta_angle_condition

# ========================================================================================
# Geometric / geographic event condition functions
#
# Detectors for geometric quantities: altitude, latitude, longitude, beta angle.
# All return closures `(u, t, integrator) -> Float64` for ContinuousCallback.
# ========================================================================================

"""
    altitude_condition(::Type{C}, h_target; R_body=6378.137) where {C<:AstroCoord}
    altitude_condition(::Type{C}, config, h_target; R_body=6378.137) where {C<:AstroCoord}

Condition function zero when spherical altitude crosses `h_target` [km].
g-function: `|r| - (R_body + h_target)`.
"""
function altitude_condition(
    ::Type{C}, h_target::Number; R_body::Number=6378.137
) where {C<:AstroCoords.AstroCoord}
    r_target = R_body + h_target
    return (u, t, integrator) -> begin
        cart = get_cartesian(u, t, integrator.p.μ, C)
        r = sqrt(cart[1]^2 + cart[2]^2 + cart[3]^2)
        return r - r_target
    end
end

function altitude_condition(
    ::Type{C},
    config::RegularizedCoordinateConfig,
    h_target::Number;
    R_body::Number=6378.137,
) where {C<:AstroCoords.AstroCoord}
    r_target = R_body + h_target
    return (u, t, integrator) -> begin
        cart = get_cartesian(u, t, integrator.p.μ, C, config)
        r = sqrt(cart[1]^2 + cart[2]^2 + cart[3]^2)
        return r - r_target
    end
end

"""
    latitude_condition(::Type{C}, lat_target) where {C<:AstroCoord}
    latitude_condition(::Type{C}, config, lat_target) where {C<:AstroCoord}

Condition function zero when geocentric latitude crosses `lat_target` [rad].
g-function: `asin(r_z / |r|) - lat_target`.
"""
function latitude_condition(::Type{C}, lat_target::Number) where {C<:AstroCoords.AstroCoord}
    return (u, t, integrator) -> begin
        cart = get_cartesian(u, t, integrator.p.μ, C)
        r = sqrt(cart[1]^2 + cart[2]^2 + cart[3]^2)
        return asin(cart[3] / r) - lat_target
    end
end

function latitude_condition(
    ::Type{C}, config::RegularizedCoordinateConfig, lat_target::Number
) where {C<:AstroCoords.AstroCoord}
    return (u, t, integrator) -> begin
        cart = get_cartesian(u, t, integrator.p.μ, C, config)
        r = sqrt(cart[1]^2 + cart[2]^2 + cart[3]^2)
        return asin(cart[3] / r) - lat_target
    end
end

@inline function _gmst(JD::Number)
    T = (JD - 2451545.0) / 36525.0
    θ =
        280.46061837 +
        360.98564736629 * (JD - 2451545.0) +
        T^2 * (0.000387933 - T / 38710000.0)
    return deg2rad(mod(θ, 360.0))
end

"""
    longitude_condition(::Type{C}, lon_target) where {C<:AstroCoord}
    longitude_condition(::Type{C}, config, lon_target) where {C<:AstroCoord}

Condition function zero when sub-satellite longitude crosses `lon_target` [rad].
Uses `sin(lon - lon_target)` to handle the ±π discontinuity. Longitude is
computed via approximate GMST rotation from ECI. Requires `integrator.p.JD`.
"""
function longitude_condition(
    ::Type{C}, lon_target::Number
) where {C<:AstroCoords.AstroCoord}
    return (u, t, integrator) -> begin
        cart = get_cartesian(u, t, integrator.p.μ, C)
        JD = integrator.p.JD + t / 86400.0
        gmst = _gmst(JD)
        lon = atan(cart[2], cart[1]) - gmst
        return sin(lon - lon_target)
    end
end

function longitude_condition(
    ::Type{C}, config::RegularizedCoordinateConfig, lon_target::Number
) where {C<:AstroCoords.AstroCoord}
    return (u, t, integrator) -> begin
        cart = get_cartesian(u, t, integrator.p.μ, C, config)
        t_phys = get_physical_time(u, t, C, config)
        JD = integrator.p.JD + t_phys / 86400.0
        gmst = _gmst(JD)
        lon = atan(cart[2], cart[1]) - gmst
        return sin(lon - lon_target)
    end
end

"""
    beta_angle_condition(::Type{C}, sun_data, β_target) where {C<:AstroCoord}
    beta_angle_condition(::Type{C}, config, sun_data, β_target) where {C<:AstroCoord}

Condition function zero when the beta angle (angle between orbit plane and Sun
direction) crosses `β_target` [rad]. g-function: `asin(ĥ · ŝ) - β_target`.
Requires `integrator.p.JD`.
"""
function beta_angle_condition(
    ::Type{C}, sun_data::AstroForceModels.ThirdBodyModel, β_target::Number
) where {C<:AstroCoords.AstroCoord}
    return (u, t, integrator) -> begin
        cart = get_cartesian(u, t, integrator.p.μ, C)
        r = SVector{3}(cart[1], cart[2], cart[3])
        v = SVector{3}(cart[4], cart[5], cart[6])
        h = cross(r, v)
        h_hat = h / norm(h)
        JD = integrator.p.JD + t / 86400.0
        sun_pos = sun_data(JD, Position()) ./ 1E3
        s_hat = SVector{3}(sun_pos[1], sun_pos[2], sun_pos[3]) / norm(sun_pos)
        β = asin(dot(h_hat, s_hat))
        return β - β_target
    end
end

function beta_angle_condition(
    ::Type{C},
    config::RegularizedCoordinateConfig,
    sun_data::AstroForceModels.ThirdBodyModel,
    β_target::Number,
) where {C<:AstroCoords.AstroCoord}
    return (u, t, integrator) -> begin
        cart = get_cartesian(u, t, integrator.p.μ, C, config)
        r = SVector{3}(cart[1], cart[2], cart[3])
        v = SVector{3}(cart[4], cart[5], cart[6])
        h = cross(r, v)
        h_hat = h / norm(h)
        t_phys = get_physical_time(u, t, C, config)
        JD = integrator.p.JD + t_phys / 86400.0
        sun_pos = sun_data(JD, Position()) ./ 1E3
        s_hat = SVector{3}(sun_pos[1], sun_pos[2], sun_pos[3]) / norm(sun_pos)
        β = asin(dot(h_hat, s_hat))
        return β - β_target
    end
end
