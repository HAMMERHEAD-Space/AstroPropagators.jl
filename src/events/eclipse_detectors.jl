export eclipse_condition

# ========================================================================================
# Eclipse event condition functions
#
# Detects eclipse transitions using AstroForceModels.shadow_model. The shadow
# factor is continuous: 0.0 = full shadow, 1.0 = full sunlight.
# ========================================================================================

"""
    eclipse_condition(::Type{C}, shadow_type, sun_data; threshold=0.5, R_Sun=..., R_Occulting=...)
    eclipse_condition(::Type{C}, config, shadow_type, sun_data; kwargs...)

Condition function zero at eclipse transitions. Returns `shadow_factor - threshold`.
Eclipse entry is a negative-going crossing; exit is positive-going.

Threshold controls what transition is detected: `≈0.5` for general shadow
entry/exit, `≈0.01` for umbra, `≈0.99` for penumbra.
"""
function eclipse_condition(
    ::Type{C},
    shadow_type::AstroForceModels.ShadowModelType,
    sun_data::AstroForceModels.ThirdBodyModel;
    threshold::Number=0.5,
    R_Sun::Number=AstroForceModels.R_SUN,
    R_Occulting::Number=AstroForceModels.R_EARTH,
) where {C<:AstroCoords.AstroCoord}
    return (u, t, integrator) -> begin
        cart = get_cartesian(u, t, integrator.p.μ, C)
        sat_pos = SVector{3}(cart[1], cart[2], cart[3])
        JD = integrator.p.JD + t / 86400.0
        sun_pos = sun_data(JD, Position()) ./ 1E3
        sf = AstroForceModels.shadow_model(
            sat_pos, sun_pos, shadow_type; R_Sun=R_Sun, R_Occulting=R_Occulting
        )
        return sf - threshold
    end
end

function eclipse_condition(
    ::Type{C},
    config::RegularizedCoordinateConfig,
    shadow_type::AstroForceModels.ShadowModelType,
    sun_data::AstroForceModels.ThirdBodyModel;
    threshold::Number=0.5,
    R_Sun::Number=AstroForceModels.R_SUN,
    R_Occulting::Number=AstroForceModels.R_EARTH,
) where {C<:AstroCoords.AstroCoord}
    return (u, t, integrator) -> begin
        cart = get_cartesian(u, t, integrator.p.μ, C, config)
        sat_pos = SVector{3}(cart[1], cart[2], cart[3])
        t_phys = get_physical_time(u, t, C, config)
        JD = integrator.p.JD + t_phys / 86400.0
        sun_pos = sun_data(JD, Position()) ./ 1E3
        sf = AstroForceModels.shadow_model(
            sat_pos, sun_pos, shadow_type; R_Sun=R_Sun, R_Occulting=R_Occulting
        )
        return sf - threshold
    end
end

"""
    eclipse_condition(::Type{C}, srp_model::SRPAstroModel; threshold=0.5)
    eclipse_condition(::Type{C}, config, srp_model; threshold=0.5)

Convenience overload that extracts shadow model parameters from an existing
`SRPAstroModel`, ensuring consistency with the force model configuration.
"""
function eclipse_condition(
    ::Type{C}, srp_model::AstroForceModels.SRPAstroModel; threshold::Number=0.5
) where {C<:AstroCoords.AstroCoord}
    return eclipse_condition(
        C,
        srp_model.shadow_model,
        srp_model.sun_data;
        threshold=threshold,
        R_Sun=srp_model.R_Sun,
        R_Occulting=srp_model.R_Occulting,
    )
end

function eclipse_condition(
    ::Type{C},
    config::RegularizedCoordinateConfig,
    srp_model::AstroForceModels.SRPAstroModel;
    threshold::Number=0.5,
) where {C<:AstroCoords.AstroCoord}
    return eclipse_condition(
        C,
        config,
        srp_model.shadow_model,
        srp_model.sun_data;
        threshold=threshold,
        R_Sun=srp_model.R_Sun,
        R_Occulting=srp_model.R_Occulting,
    )
end
