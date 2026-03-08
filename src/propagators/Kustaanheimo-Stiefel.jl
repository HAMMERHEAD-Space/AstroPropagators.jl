export KS_EOM, KS_EOM!

"""
    KS_EOM(u, p, ϕ, models, config)

Equations of motion for the Kustaanheimo-Stiefel (KS) regularized formulation.

# Arguments
- `u::AbstractArray`: The KS state vector `[KSp₁..₄, KSv₁..₄, h, τ]`.
- `p::ComponentVector`: Parameter vector containing `μ` and `JD`.
- `ϕ::Number`: The independent variable (fictitious time).
- `models::AbstractDynamicsModel`: Force model composition.
- `config::RegularizedCoordinateConfig`: Regularization configuration (DU, TU, time element type).

# Returns
- `SVector{10}`: Instantaneous rate of change of the KS state.
"""
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

"""
    KS_EOM!(du, u, p, ϕ, models, config)

In-place version of [`KS_EOM`](@ref).
"""
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
