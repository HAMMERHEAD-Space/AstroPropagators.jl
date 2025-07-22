export StiSche_EOM, StiSche_EOM!
export StiSche_time_condition,
    end_StiSche_integration, impulsive_burn_stische!, StiSche_burn

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
        μ_nd = μ * TU^2 / DU^3
        lte1 = (μ_nd - 2.0 * r_mag * U) / (8.0 * ω^3)
        lte2 = (r_mag / (16.0 * ω^3)) * dot(KSp, ∇U_u - 2.0 * Lp)
        lte3 = (2.0 / ω^2) * dω * dot(KSp, KSv)
        dτ = lte1 - lte2 - lte3
    end

    return SVector{10}(dα[1], dα[2], dα[3], dα[4], dβ[1], dβ[2], dβ[3], dβ[4], dω, dτ)
end

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

"""
    StiSche_time_condition(u, ϕ, integrator, config; event_time=0.0)

Event function for `DifferentialEquations.jl` callbacks which triggers
when the current integration time equals a specified `event_time`.

"""
function StiSche_time_condition(
    u::AbstractVector,
    ϕ::Number,
    integrator::T,
    config::RegularizedCoordinateConfig;
    event_time::Number=0.0,
) where {T<:SciMLBase.DEIntegrator}
    t = get_stiefelscheifele_time(u, ϕ, config)
    return t - event_time
end

"""
    end_StiSche_integration(stop_time, config)

Returns a `ContinuousCallback` which terminates a `DifferentialEquations.jl`
integration when the current time equals `event_time`.

"""
function end_StiSche_integration(stop_time::Number, config::RegularizedCoordinateConfig)
    ContinuousCallback(
        (u, ϕ, integrator) ->
            StiSche_time_condition(u, ϕ, integrator, config; event_time=stop_time),
        (integrator) -> terminate!(integrator);
    )
end

"""
    impulsive_burn_stische!(integrator, ΔV, config)

Applies an impulsive maneuver `ΔV` to a Stiefel-Scheifele state within a
`DifferentialEquations.jl` integrator.

"""
function impulsive_burn_stische!(
    integrator::T, ΔV::AbstractVector, config::RegularizedCoordinateConfig
) where {T<:SciMLBase.DEIntegrator}
    # Pass phi separately to coordinate transformations
    cart_state = Cartesian(
        StiefelScheifele(integrator.u), integrator.p.μ, integrator.t, config
    )

    new_state = cart_state + SVector{6}(0, 0, 0, ΔV[1], ΔV[2], ΔV[3])
    new_cart_state = Cartesian(new_state...)

    t_maneuver = get_stiefelscheifele_time(integrator.u, integrator.t, config)

    # Create new config for the maneuver with updated t₀
    maneuver_config = RegularizedCoordinateConfig(
        config.DU, config.TU, config.W, t_maneuver, config.flag_time
    )

    # Convert back to Stiefel-Scheifele using the phi from the maneuver time
    integrator.u = params(
        StiefelScheifele(new_cart_state, integrator.p.μ, integrator.t, maneuver_config)
    )

    return nothing
end

"""
    StiSche_burn(burn_time, ΔV, config)

Returns a `ContinuousCallback` which triggers an `impulsive_burn_stische!`
maneuver at a specified `burn_time`.

"""
function StiSche_burn(
    burn_time::Number, ΔV::AbstractVector, config::RegularizedCoordinateConfig
)
    ContinuousCallback(
        (u, ϕ, integrator) ->
            StiSche_time_condition(u, ϕ, integrator, config; event_time=burn_time),
        (integrator) -> impulsive_burn_stische!(integrator, ΔV, config),
    )
end
