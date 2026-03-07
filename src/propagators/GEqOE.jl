export GEqOE_EOM, GEqOE_EOM!
export impulsive_burn_geqoe!, GEqOE_burn

"""
    GEqOE_EOM(
        u::AbstractVector,
        p::ComponentVector,
        t::Number,
        models::AbstractDynamicsModel,
        config::RegularizedCoordinateConfig,
    )

Computes the equations of motion for the Generalized Equinoctial Orbital Elements (GEqOE)
formulation.

The GEqOE state vector is `[ν, p₁, p₂, L, q₁, q₂]` where ν is the generalized mean motion,
(p₁, p₂) are generalized eccentricity components, L is the generalized mean longitude,
and (q₁, q₂) are inclination components.

# Arguments
- `u::AbstractVector`: The GEqOE state vector `[ν, p₁, p₂, L, q₁, q₂]`.
- `p::ComponentVector`: Parameter vector containing gravitational parameter `μ` and Julian Date `JD`.
- `t::Number`: The current physical time.
- `models::AbstractDynamicsModel`: Force model composition.
- `config::RegularizedCoordinateConfig`: Configuration containing the perturbing potential `W`.

# Returns
- `SVector{6}`: Time derivatives `[dν, dp₁, dp₂, dL, dq₁, dq₂]`.
"""
function GEqOE_EOM(
    u::AbstractVector,
    ps::ComponentVector,
    t::Number,
    models::AbstractDynamicsModel,
    config::RegularizedCoordinateConfig,
)
    ν, p₁, p₂, L_elem, q₁, q₂ = u
    μ::Number = ps.μ

    ##################################################
    # 1. GEqOE → Cartesian state
    ##################################################
    u_cart = Cartesian(GEqOE(u), μ, config)

    r_vec = SVector{3}(u_cart[1], u_cart[2], u_cart[3])
    v_vec = SVector{3}(u_cart[4], u_cart[5], u_cart[6])

    r_mag = √(sum(abs2.(r_vec)))
    ṙ_dot = dot(r_vec, v_vec) / r_mag
    h_vec = cross(r_vec, v_vec)
    h = √(sum(abs2.(h_vec)))

    eᵣ = r_vec / r_mag
    eₕ = h_vec / h
    eₐ = cross(eₕ, eᵣ)

    ##################################################
    # 2. Derived quantities from GEqOE elements
    ##################################################
    a = (μ / ν^2)^(1.0 / 3.0)
    g² = p₁^2 + p₂^2
    ρ = a * (1.0 - g²)
    c = √(μ * ρ)
    β = √(1.0 - g²)
    α = 1.0 / (1.0 + β)

    ς = r_mag / ρ
    ς̃ = 1.0 + ς

    γ = 1.0 + q₁^2 + q₂^2
    eₓ = SVector{3}((1.0 - q₁^2 + q₂^2) / γ, (2.0 * q₁ * q₂) / γ, (-2.0 * q₁) / γ)
    eᵧ = SVector{3}((2.0 * q₁ * q₂) / γ, (1.0 + q₁^2 - q₂^2) / γ, (2.0 * q₂) / γ)

    sL = dot(eᵣ, eᵧ)
    cL = dot(eᵣ, eₓ)

    ŵh = q₁ * cL - q₂ * sL

    ##################################################
    # 3. Force evaluation
    ##################################################
    kep_model = KeplerianGravityAstroModel(; μ=μ)

    U = potential(u_cart, ps, t, models.gravity_model) - potential(u_cart, ps, t, kep_model)

    U_t = potential_time_derivative(u_cart, ps, t, models.gravity_model)

    a_grav_pert =
        acceleration(u_cart, ps, t, models.gravity_model) -
        acceleration(u_cart, ps, t, kep_model)

    P_vec = AstroForceModels.sum_accelerations(u_cart, ps, t, models.perturbing_models)

    F_vec = a_grav_pert + P_vec

    Fr = dot(F_vec, eᵣ)
    Fh = dot(F_vec, eₕ)
    Pr = dot(P_vec, eᵣ)
    Pf = dot(P_vec, eₐ)

    ##################################################
    # 4. Energy rate (Eq. 43)
    ##################################################
    Ė = U_t + ṙ_dot * Pr + (h / r_mag) * Pf

    ##################################################
    # 5. Equations of motion (Eqs. 42, 45, 46, 48, 49, 50)
    ##################################################
    hc_r2 = (h - c) / r_mag^2
    rh_whFh = (r_mag / h) * ŵh * Fh
    twoU_rFr = 2.0 * U - r_mag * Fr

    # ν̇ (Eq. 42)
    dν = -3.0 * (ν / μ^2)^(1.0 / 3.0) * Ė

    # ṗ₁ (Eq. 45)
    dp₁ =
        p₂ * (hc_r2 - rh_whFh) +
        (r_mag * ṙ_dot / c^2) * (p₁ + ς̃ * p₂ + ς * cL) * twoU_rFr +
        (r_mag / μ) * (ς * p₁ + ς̃ * sL) * Ė

    # ṗ₂ (Eq. 46)
    dp₂ =
        p₁ * (rh_whFh - hc_r2) +
        (r_mag * ṙ_dot / c^2) * (p₂ - ς̃ * p₁ - ς * sL) * twoU_rFr +
        (r_mag / μ) * (ς * p₂ + ς̃ * cL) * Ė

    # L̇ (Eq. 48)
    dL =
        ν + hc_r2 - rh_whFh +
        (r_mag * ṙ_dot * c * ς̃ * α / μ^2) * Ė +
        (1.0 / c) * (1.0 / α + α * (1.0 - r_mag / a)) * twoU_rFr

    # q̇₁ (Eq. 49)
    dq₁ = (r_mag * γ * Fh * sL) / (2.0 * h)

    # q̇₂ (Eq. 50)
    dq₂ = (r_mag * γ * Fh * cL) / (2.0 * h)

    return SVector{6}(dν, dp₁, dp₂, dL, dq₁, dq₂)
end

"""
    GEqOE_EOM!(
        du::AbstractVector,
        u::AbstractVector,
        p::ComponentVector,
        t::Number,
        models::AbstractDynamicsModel,
        config::RegularizedCoordinateConfig,
    )

In-place version of [`GEqOE_EOM`](@ref).

# Arguments
- `du::AbstractVector`: In-place vector to store the time derivatives.
- `u::AbstractVector`: The GEqOE state vector `[ν, p₁, p₂, L, q₁, q₂]`.
- `p::ComponentVector`: Parameter vector containing gravitational parameter `μ` and Julian Date `JD`.
- `t::Number`: The current physical time.
- `models::AbstractDynamicsModel`: Force model composition.
- `config::RegularizedCoordinateConfig`: Configuration containing the perturbing potential `W`.

# Returns
- `nothing`
"""
function GEqOE_EOM!(
    du::AbstractVector,
    u::AbstractVector,
    p::ComponentVector,
    t::Number,
    models::AbstractDynamicsModel,
    config::RegularizedCoordinateConfig,
)
    du .= GEqOE_EOM(u, p, t, models, config)

    return nothing
end

"""
    impulsive_burn_geqoe!(integrator, ΔV, config; frame=InertialFrame())

Apply an instantaneous velocity change to a GEqOE state within a
`DifferentialEquations.jl` integrator.

The `ΔV` vector is interpreted in the reference frame specified by `frame`
(see `InertialFrame`, `RTNFrame`, `VNBFrame`).
"""
function impulsive_burn_geqoe!(
    integrator::T,
    ΔV::AbstractVector,
    config::RegularizedCoordinateConfig;
    frame::AbstractThrustFrame=InertialFrame(),
) where {T<:SciMLBase.DEIntegrator}
    cart_state = Cartesian(GEqOE(integrator.u), integrator.p.μ, config)
    ΔV_inertial = transform_to_inertial(SVector{3}(ΔV[1], ΔV[2], ΔV[3]), cart_state, frame)

    new_state =
        cart_state + SVector{6}(0, 0, 0, ΔV_inertial[1], ΔV_inertial[2], ΔV_inertial[3])
    new_cart_state = Cartesian(new_state...)

    integrator.u = params(GEqOE(new_cart_state, integrator.p.μ, config))

    return nothing
end

"""
    GEqOE_burn(burn_time, ΔV, config; frame=InertialFrame())

Returns a `ContinuousCallback` which triggers an [`impulsive_burn_geqoe!`](@ref)
maneuver at a specified `burn_time`.

The `ΔV` is interpreted in the given `frame` (default: `InertialFrame`).
"""
function GEqOE_burn(
    burn_time::Number,
    ΔV::AbstractVector,
    config::RegularizedCoordinateConfig;
    frame::AbstractThrustFrame=InertialFrame(),
)
    ContinuousCallback(
        (u, t, integrator) -> t - burn_time,
        (integrator) -> impulsive_burn_geqoe!(integrator, ΔV, config; frame),
    )
end
