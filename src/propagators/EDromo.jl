export EDromo_EOM, EDromo_EOM!
export EDromo_time_condition, end_EDromo_integration
export EDromo_burn, impulsive_burn_edromo!

"""
function EDromo_EOM(
    u::AbstractArray,
    p::ComponentVector,
    ϕ::Number,
    models::AbstractDynamicsModel,
    config::RegularizedCoordinateConfig,
)

Computes the equations of motion for the EDromo formulation. The formulation is based on [1] and [2].

# Arguments
-`u::AbstractArray`: The EDromo state vector `[ζ₁, ζ₂, ζ₃, ζ₄, ζ₅, ζ₆, ζ₇, ζ₈]`.
-`p::ComponentVector`: A `ComponentVector` containing parameters for the propagation,
  including the gravitational parameter `μ` and the initial Julian Date `JD`.
-`ϕ::Number`: The independent variable, which is the fictitious time.
-`models::AbstractDynamicsModel`: Tuple of the acceleration models.
-`config::RegularizedCoordinateConfig`: Configuration struct containing EDromo formulation parameters.

# References

[1] Baù, G., Bombardelli, C., Peláez, J., and Lorenzini, E., "Nonsingular
    orbital elements for special perturbations in the two-body problem".
    MNRAS 454(3), pp. 2890-2908. 2015.
[2] Amato, Davide. "THALASSA: Orbit propagator for near-Earth and cislunar space." Astrophysics Source Code Library (2019): ascl-1905.
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
function EDromo_EOM!(
    du::AbstractVector,
    u::AbstractArray,
    p::ComponentVector,
    ϕ::Number,
    models::AbstractDynamicsModel,
    config::RegularizedCoordinateConfig,
)

Computes the equations of motion for the EDromo formulation. The formulation is based on [1] and [2].

# Arguments
-`du::AbstractVector`: In-place vector to store the instantenous rate of change of the current state with respect to time.
-`u::AbstractArray`: The EDromo state vector `[ζ₁, ζ₂, ζ₃, ζ₄, ζ₅, ζ₆, ζ₇, ζ₈]`.
-`p::ComponentVector`: A `ComponentVector` containing parameters for the propagation,
  including the gravitational parameter `μ` and the initial Julian Date `JD`.
-`ϕ::Number`: The independent variable, which is the fictitious time.
-`models::AbstractDynamicsModel`: Tuple of the acceleration models.
-`config::RegularizedCoordinateConfig`: Configuration struct containing EDromo formulation parameters.

# References

[1] Baù, G., Bombardelli, C., Peláez, J., and Lorenzini, E., "Nonsingular
    orbital elements for special perturbations in the two-body problem".
    MNRAS 454(3), pp. 2890-2908. 2015.
[2] Amato, Davide. "THALASSA: Orbit propagator for near-Earth and cislunar space." Astrophysics Source Code Library (2019): ascl-1905.
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

"""
    EDromo_time_condition(u, ϕ, integrator, config; event_time=0.0)

Event function for `DifferentialEquations.jl` callbacks which triggers
when the current integration time equals a specified `event_time`.

"""
function EDromo_time_condition(
    u::AbstractVector,
    ϕ::Number,
    integrator::T,
    config::RegularizedCoordinateConfig;
    event_time::Number=0.0,
) where {T<:SciMLBase.DEIntegrator}
    # Pass phi separately
    t = get_EDromo_time(u, ϕ, config)
    return t - event_time
end

"""
    end_EDromo_integration(stop_time, config)

Returns a `ContinuousCallback` which terminates a `DifferentialEquations.jl`
integration when the current time equals `event_time`.

"""
function end_EDromo_integration(stop_time::Number, config::RegularizedCoordinateConfig)
    ContinuousCallback(
        (u, ϕ, integrator) ->
            EDromo_time_condition(u, ϕ, integrator, config; event_time=stop_time),
        (integrator) -> terminate!(integrator);
    )
end

"""
    impulsive_burn_edromo!(integrator, ΔV, config; frame=InertialFrame())

Apply an instantaneous velocity change to an EDromo state within a
`DifferentialEquations.jl` integrator.

The `ΔV` vector is interpreted in the reference frame specified by `frame`
(see [`InertialFrame`](@ref), [`RTNFrame`](@ref), [`VNBFrame`](@ref)).
"""
function impulsive_burn_edromo!(
    integrator::T,
    ΔV::AbstractVector,
    config::RegularizedCoordinateConfig;
    frame::AbstractThrustFrame=InertialFrame(),
) where {T<:SciMLBase.DEIntegrator}
    cart_state = Cartesian(EDromo(integrator.u), integrator.p.μ, integrator.t, config)
    ΔV_inertial = transform_to_inertial(SVector{3}(ΔV[1], ΔV[2], ΔV[3]), cart_state, frame)

    new_state = cart_state + SVector{6}(0, 0, 0, ΔV_inertial[1], ΔV_inertial[2], ΔV_inertial[3])
    new_cart_state = Cartesian(new_state...)

    t_maneuver = get_EDromo_time(integrator.u, integrator.t, config)

    maneuver_config = RegularizedCoordinateConfig(
        config.DU, config.TU, config.W, t_maneuver, config.flag_time
    )

    integrator.u = params(
        EDromo(new_cart_state, integrator.p.μ, integrator.t, maneuver_config)
    )

    return nothing
end

"""
    EDromo_burn(burn_time, ΔV, config; frame=InertialFrame())

Returns a `ContinuousCallback` which triggers an `impulsive_burn_edromo!`
maneuver at a specified `burn_time`.

The `ΔV` is interpreted in the given `frame` (default: [`InertialFrame`](@ref)).
"""
function EDromo_burn(
    burn_time::Number,
    ΔV::AbstractVector,
    config::RegularizedCoordinateConfig;
    frame::AbstractThrustFrame=InertialFrame(),
)
    ContinuousCallback(
        (u, ϕ, integrator) ->
            EDromo_time_condition(u, ϕ, integrator, config; event_time=burn_time),
        (integrator) -> impulsive_burn_edromo!(integrator, ΔV, config; frame),
    )
end
