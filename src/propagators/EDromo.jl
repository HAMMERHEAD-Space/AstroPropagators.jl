export EDromo_EOM, EDromo_EOM!

"""
    EDromo_EOM(u, p, ϕ, models, config)

Equations of motion for the EDromo formulation.

# Arguments
- `u::AbstractArray`: The EDromo state vector `[ζ₁, ζ₂, ζ₃, ζ₄, ζ₅, ζ₆, ζ₇, ζ₈]`.
- `p::ComponentVector`: Parameter vector containing `μ` and `JD`.
- `ϕ::Number`: The independent variable (fictitious time).
- `models::AbstractDynamicsModel`: Force model composition.
- `config::RegularizedCoordinateConfig`: Regularization configuration (DU, TU, time element type).

# Returns
- `SVector{8}`: Instantaneous rate of change of the EDromo state.
"""
function EDromo_EOM(
    u::AbstractArray,
    p::ComponentVector,
    ϕ::Number,
    models::AbstractDynamicsModel,
    config::RegularizedCoordinateConfig,
)
    ζ1, ζ2, ζ3, ζ4, ζ5, ζ6, ζ7, _ = u

    # Extract parameters from config
    DU, TU, flag_time = config.DU, config.TU, config.flag_time

    ##################################################
    #* 1. Auxiliary Quantities (1)
    ##################################################
    sϕ, cϕ = sincos(ϕ)

    ρ = 1.0 - ζ1 * cϕ - ζ2 * sϕ
    r_mag = ζ3 * ρ
    χ = ζ1 * sϕ - ζ2 * cϕ
    ϵ = √(1.0 - ζ1^2 - ζ2^2)

    cν = (cϕ - ζ1 + (χ * ζ2) / (ϵ + 1.0)) / ρ
    sν = (sϕ - ζ2 - (χ * ζ1) / (ϵ + 1.0)) / ρ

    μ::Number = p.μ

    ##################################################
    #* 2. Inertial State and Time
    ##################################################
    # Pass phi separately to coordinate transformations
    u_cart = Cartesian(EDromo(u), μ, ϕ, config)
    tt = get_EDromo_time(u, ϕ, config)

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

    # Transform to orbital frame
    ∇Uᵣ_orbital = inert2orb_edromo(∇Uᵣ_inertial, u, cν, sν)

    ##################################################
    #* 4. Non-Potential Based Perturbations
    ##################################################
    p_inertial =
        AstroForceModels.sum_accelerations(u_cart, p, tt, models.perturbing_models) *
        (TU^2 / DU)

    # Transform to orbital frame
    p_orbital = inert2orb_edromo(p_inertial, u, cν, sν)

    ##################################################
    #* 5. Auxiliary Quantities (2)
    ##################################################
    F = p_orbital + ∇Uᵣ_orbital

    η = √(ϵ^2 - 2.0 * ζ3 * ρ^2 * U)
    aux0 = p_orbital[1] * χ + p_orbital[2] * η + ∇Uₜ * √(ζ3) * ρ

    dζ3 = 2.0 * ζ3^3 * aux0
    L3 = dζ3 / (2.0 * ζ3)

    aux1 = ((2.0 * U - F[1] * r_mag) * (2.0 - ρ + ϵ) * r_mag) / (ϵ * (ϵ + 1.0))
    aux2 = (L3 * χ * (ρ - ϵ)) / (ϵ * (ϵ + 1.0))
    wz = (η - ϵ) / ρ + aux1 + aux2

    ##################################################
    #* 6. Equations of Motion
    ##################################################
    aux3 = (F[1] * r_mag - 2.0 * U) * r_mag
    dζ1 = aux3 * sϕ + L3 * ((1.0 + ρ) * cϕ - ζ1)
    dζ2 = -aux3 * cϕ + L3 * ((1.0 + ρ) * sϕ - ζ2)

    aux4 = (F[3] * r_mag^2) / (2.0 * η)
    dζ4 = aux4 * (ζ7 * cν - ζ6 * sν) + 0.5 * wz * ζ5
    dζ5 = aux4 * (ζ6 * cν + ζ7 * sν) - 0.5 * wz * ζ4
    dζ6 = -aux4 * (ζ5 * cν - ζ4 * sν) + 0.5 * wz * ζ7
    dζ7 = -aux4 * (ζ4 * cν + ζ5 * sν) - 0.5 * wz * ζ6

    if flag_time isa PhysicalTime
        dζ8 = √(ζ3) * r_mag
    elseif flag_time isa ConstantTime
        dζ8 = ζ3^1.5 * (aux3 + (χ - 1.5 * ϕ) * dζ3 / ζ3)
    elseif flag_time isa LinearTime
        dζ8 = ζ3^1.5 * (1.0 + aux3 + 2.0 * L3 * χ)
    end

    return SVector{8}(dζ1, dζ2, dζ3, dζ4, dζ5, dζ6, dζ7, dζ8)
end

"""
    EDromo_EOM!(du, u, p, ϕ, models, config)

In-place version of [`EDromo_EOM`](@ref).
"""
function EDromo_EOM!(
    du::AbstractVector,
    u::AbstractVector,
    p::ComponentVector,
    ϕ::Number,
    models::AbstractDynamicsModel,
    config::RegularizedCoordinateConfig,
)
    du .= EDromo_EOM(u, p, ϕ, models, config)

    return nothing
end

"""
    inert2orb_edromo(v_inertial, u, cν, sν) -> SVector{3}

Transforms a vector from the inertial frame to the EDROMO orbital frame.

The orbital frame is defined by the radial, tangential, and normal directions.

# Args

- `v_inertial::AbstractVector`: The vector in the inertial reference frame.
- `u::AbstractArray`: The EDromo state vector `[ζ₁, ζ₂, ζ₃, ζ₄, ζ₅, ζ₆, ζ₇, ζ₈]`.
- `cν::Number`: The cosine of the auxiliary angle `ν`.
- `sν::Number`: The sine of the auxiliary angle `ν`.

# Returns

- `SVector{3}`: The transformed vector in the EDROMO orbital frame.
"""
@inline function inert2orb_edromo(
    v_inertial::AbstractVector, u::AbstractArray, cν::Number, sν::Number
)
    _, _, _, ζ4, ζ5, ζ6, ζ7, _ = u

    # Inertial -> Intermediate rotation matrix (transposed)
    R1_T =
        2.0 * SMatrix{3,3}(
            0.5-ζ5^2-ζ6^2,
            ζ4*ζ5-ζ6*ζ7,
            ζ4*ζ6+ζ5*ζ7,
            ζ4*ζ5+ζ6*ζ7,
            0.5-ζ4^2-ζ6^2,
            ζ5*ζ6-ζ4*ζ7,
            ζ4*ζ6-ζ5*ζ7,
            ζ5*ζ6+ζ4*ζ7,
            0.5-ζ4^2-ζ5^2,
        )

    # Intermediate -> Orbital rotation matrix
    R2 = SMatrix{3,3}(cν, -sν, 0.0, sν, cν, 0.0, 0.0, 0.0, 1.0)

    return R2 * R1_T * v_inertial
end

