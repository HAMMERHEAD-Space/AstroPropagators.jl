export KS_EOM, KS_EOM!
export KS_time_condition, end_KS_integration, impulsive_burn_ks!, KS_burn

function KS_EOM(
    u::AbstractArray,
    p::ComponentVector,
    ϕ::Number,
    models::AbstractDynamicsModel,
    config::RegularizedCoordinateConfig,
)
    KSp1, KSp2, KSp3, KSp4, KSv1, KSv2, KSv3, KSv4, h, τ = u
    μ::Number = p.μ

    # Extract parameters from config
    DU, TU = config.DU, config.TU

    ##################################################
    #* 1. Auxiliary Quantities
    ##################################################
    KSp = SVector{4}(KSp1, KSp2, KSp3, KSp4)
    KSv = SVector{4}(KSv1, KSv2, KSv3, KSv4)
    r_mag = KSp1^2 + KSp2^2 + KSp3^2 + KSp4^2
    L = KS_matrix(KSp)

    ##################################################
    #* 2. Cartesian State and Time
    ##################################################
    u_cart = Cartesian(KustaanheimoStiefel(u), μ, config)
    tt = get_KS_time(u, config)

    ##################################################
    #* 3. Potential Based Perturbations
    ##################################################
    kep_model = KeplerianGravityAstroModel(; μ=μ)
    U =
        (
            potential(u_cart, p, tt, models.gravity_model) -
            potential(u_cart, p, tt, kep_model)
        ) * (TU^2 / DU^2)
    ∇Uᵣ_inertial =
        (
            acceleration(u_cart, p, tt, models.gravity_model) -
            acceleration(u_cart, p, tt, kep_model)
        ) * (TU^2 / DU)
    ∇Uₜ = potential_time_derivative(u_cart, p, tt, models.gravity_model) * (TU^3 / DU^2)

    ##################################################
    #* 4. Non-Potential Based Perturbations
    ##################################################
    p_inertial =
        AstroForceModels.sum_accelerations(u_cart, p, tt, models.perturbing_models) *
        (TU^2 / DU)

    ##################################################
    #* 5. Equations of Motion
    ##################################################
    # Total perturbation in R^4
    F_4d = SVector{4}(
        ∇Uᵣ_inertial[1] + p_inertial[1],
        ∇Uᵣ_inertial[2] + p_inertial[2],
        ∇Uᵣ_inertial[3] + p_inertial[3],
        0.0,
    )
    p_inertial_4d = SVector{4}(p_inertial[1], p_inertial[2], p_inertial[3], 0.0)

    # Derivatives
    dKSp = KSv
    dKSv = -0.5 * (KSp * (h + U) - r_mag * (L' * F_4d))
    dh = -r_mag * ∇Uₜ - 2.0 * dot(KSv, L' * p_inertial_4d)

    if config.flag_time isa PhysicalTime
        dτ = r_mag
    elseif config.flag_time isa LinearTime
        μ_non_dim = μ * TU^2 / DU^3
        lte1 = (μ_non_dim - 2.0 * r_mag * U) / (2.0 * h)
        lte20 = 2.0 * (L' * -F_4d)
        lte2 = (r_mag / (4.0 * h)) * dot(KSp, lte20)
        lte3 = dh / (h^2) * dot(KSp, KSv)
        dτ = lte1 - lte2 - lte3
    else
        error(
            "Time flag in RegularizedCoordinateConfig not supported by Kustaanheimo-Stiefel formulation.",
        )
    end

    return SVector{10}(
        dKSp[1], dKSp[2], dKSp[3], dKSp[4], dKSv[1], dKSv[2], dKSv[3], dKSv[4], dh, dτ
    )
end

function KS_EOM!(
    du::AbstractArray,
    u::AbstractArray,
    p::ComponentVector,
    ϕ::Number,
    models::AbstractDynamicsModel,
    config::RegularizedCoordinateConfig,
)
    du .= KS_EOM(u, p, ϕ, models, config)

    return nothing
end

function KS_matrix(u::AbstractVector)
    return SMatrix{4,4}(
        u[1],
        u[2],
        u[3],
        u[4],
        -u[2],
        u[1],
        u[4],
        -u[3],
        -u[3],
        -u[4],
        u[1],
        u[2],
        u[4],
        -u[3],
        u[2],
        -u[1],
    )
end

"""
    KS_time_condition(u, ϕ, integrator, config; event_time=0.0)

Event function for `DifferentialEquations.jl` callbacks which triggers
when the current integration time equals a specified `event_time`.

"""
function KS_time_condition(
    u::AbstractVector,
    ϕ::Number,
    integrator::T,
    config::RegularizedCoordinateConfig;
    event_time::Number=0.0,
) where {T<:SciMLBase.DEIntegrator}
    t = get_KS_time(u, config)
    return t - event_time
end

"""
    end_KS_integration(stop_time, config)

Returns a `ContinuousCallback` which terminates a `DifferentialEquations.jl`
integration when the current time equals `event_time`.

"""
function end_KS_integration(stop_time::Number, config::RegularizedCoordinateConfig)
    ContinuousCallback(
        (u, ϕ, integrator) ->
            KS_time_condition(u, ϕ, integrator, config; event_time=stop_time),
        (integrator) -> terminate!(integrator);
    )
end

"""
    impulsive_burn_ks!(integrator, ΔV, config)

Applies an impulsive maneuver `ΔV` to a KS state within a
`DifferentialEquations.jl` integrator.

"""
function impulsive_burn_ks!(
    integrator::T, ΔV::AbstractVector, config::RegularizedCoordinateConfig
) where {T<:SciMLBase.DEIntegrator}
    # Pass phi separately to coordinate transformations
    cart_state = Cartesian(KustaanheimoStiefel(integrator.u), integrator.p.μ, config)

    new_state = cart_state + SVector{6}(0, 0, 0, ΔV[1], ΔV[2], ΔV[3])
    new_cart_state = Cartesian(new_state...)

    t_maneuver = get_KS_time(integrator.u, config)

    # Create new config for the maneuver with updated t₀
    maneuver_config = RegularizedCoordinateConfig(
        config.DU, config.TU, config.W, t_maneuver, config.flag_time
    )

    # Convert back to KS using the phi from the maneuver time
    integrator.u = params(
        KustaanheimoStiefel(new_cart_state, integrator.p.μ, maneuver_config)
    )

    return nothing
end

"""
    KS_burn(burn_time, ΔV, config)

Returns a `ContinuousCallback` which triggers an `impulsive_burn_ks!`
maneuver at a specified `burn_time`.

"""
function KS_burn(burn_time::Number, ΔV::AbstractVector, config::RegularizedCoordinateConfig)
    ContinuousCallback(
        (u, ϕ, integrator) ->
            KS_time_condition(u, ϕ, integrator, config; event_time=burn_time),
        (integrator) -> impulsive_burn_ks!(integrator, ΔV, config),
    )
end
