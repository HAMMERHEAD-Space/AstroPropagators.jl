export StiSche_EOM, StiSche_EOM!

"""
    StiSche_EOM(u, p, ϕ, models, config)

Equations of motion for the Stiefel-Scheifele regularized formulation.

# Arguments
- `u::AbstractArray`: The Stiefel-Scheifele state vector `[α₁..₄, β₁..₄, ω, τ]`.
- `p::ComponentVector`: Parameter vector containing `μ` and `JD`.
- `ϕ::Number`: The independent variable (fictitious time).
- `models::AbstractDynamicsModel`: Force model composition.
- `config::RegularizedCoordinateConfig`: Regularization configuration (DU, TU, time element type).

# Returns
- `SVector{10}`: Instantaneous rate of change of the Stiefel-Scheifele state.
"""
function StiSche_EOM(
    u::AbstractArray,
    p::ComponentVector,
    ϕ::Number,
    models::AbstractDynamicsModel,
    config::RegularizedCoordinateConfig,
)
    α1, α2, α3, α4, β1, β2, β3, β4, ω, _ = u
    μ::Number = p.μ

    # Extract parameters from config
    DU, TU = config.DU, config.TU

    ##################################################
    #* 1. Auxiliary Quantities (1)
    ##################################################
    sϕ2, cϕ2 = sincos(0.5 * ϕ)
    α = SVector{4}(α1, α2, α3, α4)
    β = SVector{4}(β1, β2, β3, β4)

    ##################################################
    #* 2. Position and Time in Inertial Frame
    ##################################################
    KSp = α * cϕ2 + β * sϕ2
    KSv = 0.5 * (-α * sϕ2 + β * cϕ2)
    r_mag = dot(KSp, KSp)

    u_cart = Cartesian(StiefelScheifele(u), μ, ϕ, config)
    tt = get_stiefelscheifele_time(u, ϕ, config)

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
    #* 5. Auxiliary Quantities (2)
    ##################################################
    p_inertial_4d = SVector{4}(p_inertial[1], p_inertial[2], p_inertial[3], 0.0)
    Lp = KS_matrix(KSp)' * p_inertial_4d

    ∇Uᵣ_inertial_4d = SVector{4}(∇Uᵣ_inertial[1], ∇Uᵣ_inertial[2], ∇Uᵣ_inertial[3], 0.0)
    ∇U_u = -2.0 * KS_matrix(KSp)' * ∇Uᵣ_inertial_4d

    ##################################################
    #* 6. Equations of Motion
    ##################################################
    dω = -r_mag / (8.0 * ω^2) * ∇Uₜ - (0.5 / ω) * dot(KSv, Lp)

    aux =
        (0.5 / ω^2) * (0.5 * U * KSp + r_mag / 4.0 * (∇U_u - 2.0 * Lp)) + 2.0 / ω * dω * KSv

    dα = aux * sϕ2
    dβ = -aux * cϕ2

    if config.flag_time isa PhysicalTime
        dτ = 0.5 * r_mag / ω
    elseif config.flag_time isa LinearTime
        μ_non_dim = μ * TU^2 / DU^3
        lte1 = (μ_non_dim - 2.0 * r_mag * U) / (8.0 * ω^3)
        lte2 = (r_mag / (16.0 * ω^3)) * dot(KSp, ∇U_u - 2.0 * Lp)
        lte3 = (2.0 / ω^2) * dω * dot(KSp, KSv)
        dτ = lte1 - lte2 - lte3
    else
        error(
            "Time flag in RegularizedCoordinateConfig not supported by Stiefel-Scheifele formulation.",
        )
    end

    return SVector{10}(dα[1], dα[2], dα[3], dα[4], dβ[1], dβ[2], dβ[3], dβ[4], dω, dτ)
end

"""
    StiSche_EOM!(du, u, p, ϕ, models, config)

In-place version of [`StiSche_EOM`](@ref).
"""
function StiSche_EOM!(
    du::AbstractArray,
    u::AbstractArray,
    p::ComponentVector,
    ϕ::Number,
    models::AbstractDynamicsModel,
    config::RegularizedCoordinateConfig,
)
    du .= StiSche_EOM(u, p, ϕ, models, config)
    return nothing
end
