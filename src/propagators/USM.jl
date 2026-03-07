export USM7_EOM, USM7_EOM!
"""
    USM7_EOM(u, p, t, models)

Unified State Model (quaternions) propagation schema for orbital trajectories.

# Arguments
- `u::AbstractVector`: The current USM7 state `[C, Rf1, Rf2, ϵO1, ϵO2, ϵO3, η0]`.
- `p::ComponentVector`: The parameter vector containing `μ` and `JD`.
- `t::Number`: The current time.
- `models::AbstractDynamicsModel`: Force model composition.

# Returns
- `SVector{7}`: Instantaneous rate of change of the USM7 state.
"""
function USM7_EOM(
    u::AbstractVector, p::ComponentVector, t::Number, models::AbstractDynamicsModel
)
    C, Rf1, Rf2, ϵO1, ϵO2, ϵO3, η0 = u

    μ::Number = p.μ

    sinλ = (2 * ϵO3 * η0) / (ϵO3^2 + η0^2)
    cosλ = (η0^2 - ϵO3^2) / (ϵO3^2 + η0^2)
    l = (ϵO1 * ϵO3 - ϵO2 * η0) / (ϵO3^2 + η0^2)
    ve2 = C - Rf1 * sinλ + Rf2 * cosλ
    ω3 = (C * ve2^2) / μ
    ρ = C / ve2

    u_cart = Cartesian(USM7(u), p.μ)

    fe = inertial_to_RTN(
        build_dynamics_model(u_cart, p, t, models) -
        acceleration(u_cart, p, t, KeplerianGravityAstroModel(; μ=μ)),
        u_cart,
    )

    ω1 = fe[3] / ve2

    return SVector{7}(
        -ρ * fe[2],
        fe[1] * cosλ - fe[2] * (1.0 + ρ) * sinλ - fe[3] * l * (Rf2 / ve2),
        fe[1] * sinλ + fe[2] * (1.0 + ρ) * cosλ + fe[3] * l * (Rf1 / ve2),
        0.5 * (ω3 * ϵO2 + ω1 * η0),
        0.5 * (-ω3 * ϵO1 + ω1 * ϵO3),
        0.5 * (-ω1 * ϵO2 + ω3 * η0),
        0.5 * (-ω1 * ϵO1 - ω3 * ϵO3),
    )
end

"""
    USM7_EOM!(du, u, p, t, models)

In-place version of [`USM7_EOM`](@ref).
"""
function USM7_EOM!(
    du::AbstractVector,
    u::AbstractVector,
    p::ComponentVector,
    t::Number,
    models::AbstractDynamicsModel,
)
    du .= USM7_EOM(u, p, t, models)

    return nothing
end

export USM6_EOM, USM6_EOM!
"""
    USM6_EOM(u, p, t, models)

Unified State Model (MRPs) propagation schema for orbital trajectories.

# Arguments
- `u::AbstractVector`: The current USM6 state `[C, Rf1, Rf2, σ1, σ2, σ3]`.
- `p::ComponentVector`: The parameter vector containing `μ` and `JD`.
- `t::Number`: The current time.
- `models::AbstractDynamicsModel`: Force model composition.

# Returns
- `SVector{6}`: Instantaneous rate of change of the USM6 state.
"""
function USM6_EOM(
    u::AbstractVector, p::ComponentVector, t::Number, models::AbstractDynamicsModel
)
    C, Rf1, Rf2, σ1, σ2, σ3 = u
    μ::Number = p.μ

    σ² = σ1^2 + σ2^2 + σ3^2
    D = 4.0 * σ3^2 + (1.0 - σ²)^2

    u_cart = Cartesian(USM6(u), μ)

    sinλ = (4.0 * σ3 * (1.0 - σ²)) / D
    cosλ = ((1.0 - σ²)^2 - 4.0 * σ3^2) / D
    l = (4.0 * σ1 * σ3 - 2.0 * σ2 * (1.0 - σ²)) / D

    ve2 = C - Rf1 * sinλ + Rf2 * cosλ

    ω3 = (C * ve2^2) / μ

    ρ = C / ve2

    fe = inertial_to_RTN(
        build_dynamics_model(u_cart, p, t, models) -
        acceleration(u_cart, p, t, KeplerianGravityAstroModel(; μ=μ)),
        u_cart,
    )

    ω1 = fe[3] / ve2

    return SVector{6}(
        -ρ * fe[2],
        fe[1] * cosλ - fe[2] * (1.0 + ρ) * sinλ - fe[3] * l * (Rf2 / ve2),
        fe[1] * sinλ + fe[2] * (1.0 + ρ) * cosλ + fe[3] * l * (Rf1 / ve2),
        0.25 * ((1.0 - σ² + 2.0 * σ1^2) * ω1 + 2.0 * (σ1 * σ3 + σ2) * ω3),
        0.25 * (2.0 * (σ2 * σ1 + σ3) * ω1 + 2.0 * (σ2 * σ3 - σ1) * ω3),
        0.25 * (2.0 * (σ3 * σ1 - σ2) * ω1 + (1.0 - σ² + 2.0 * σ3^2) * ω3),
    )
end

"""
    USM6_EOM!(du, u, p, t, models)

In-place version of [`USM6_EOM`](@ref).
"""
function USM6_EOM!(
    du::AbstractVector,
    u::AbstractVector,
    p::ComponentVector,
    t::Number,
    models::AbstractDynamicsModel,
)
    du .= USM6_EOM(u, p, t, models)

    return nothing
end

export USMEM_EOM, USMEM_EOM!
"""
    USMEM_EOM(u, p, t, models; Φ_tol=1E-8)

Unified State Model (exponential mapping) propagation schema for orbital trajectories.

# Arguments
- `u::AbstractVector`: The current USMEM state `[C, Rf1, Rf2, a1, a2, a3]`.
- `p::ComponentVector`: The parameter vector containing `μ` and `JD`.
- `t::Number`: The current time.
- `models::AbstractDynamicsModel`: Force model composition.

# Keyword Arguments
- `Φ_tol::Float64`: Threshold to switch to Taylor series expansion to avoid singularity.

# Returns
- `SVector{6}`: Instantaneous rate of change of the USMEM state.
"""
function USMEM_EOM(
    u::AbstractVector,
    p::ComponentVector,
    t::Number,
    models::AbstractDynamicsModel;
    Φ_tol::Float64=1E-8,
)
    C, Rf1, Rf2, a1, a2, a3 = u

    a = SVector{3}(a1, a2, a3)
    Φ = √(sum(abs2.(a)))
    a_cross = skew_sym(a)

    μ::Number = p.μ

    usm7 = USM7(USMEM(u), μ)
    _, _, _, ϵO1, ϵO2, ϵO3, η0 = usm7
    u_cart = Cartesian(usm7, μ)

    sinλ = (2 * ϵO3 * η0) / (ϵO3^2 + η0^2)
    cosλ = (η0^2 - ϵO3^2) / (ϵO3^2 + η0^2)
    l = (ϵO1 * ϵO3 - ϵO2 * η0) / (ϵO3^2 + η0^2)
    ve2 = C - Rf1 * sinλ + Rf2 * cosλ
    ω3 = (C * ve2^2) / μ
    ρ = C / ve2

    fe = inertial_to_RTN(
        build_dynamics_model(u_cart, p, t, models) -
        acceleration(u_cart, p, t, KeplerianGravityAstroModel(; μ=μ)),
        u_cart,
    )

    ω1 = fe[3] / ve2

    ω = SVector{3}(ω1, 0.0, ω3)

    dC = -ρ * fe[2]
    dRf1 = fe[1] * cosλ - fe[2] * (1.0 + ρ) * sinλ - fe[3] * l * (Rf2 / ve2)
    dRf2 = fe[1] * sinλ + fe[2] * (1.0 + ρ) * cosλ + fe[3] * l * (Rf1 / ve2)

    if Φ > Φ_tol
        da =
            (
                SMatrix{3,3}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0) +
                a_cross / 2.0 +
                1.0 / (Φ^2.0) * (1.0 - Φ / 2.0 * cot(Φ / 2.0)) * a_cross * a_cross
            ) * ω
    else
        da =
            0.5 .* (
                (12.0 - Φ^2.0) / 6.0 * ω - cross(ω, a) -
                dot(ω, a) * ((60.0 + Φ^2.0) / 360.0) * a
            )
    end

    return SVector{6}(dC, dRf1, dRf2, da[1], da[2], da[3])
end

"""
    USMEM_EOM!(du, u, p, t, models; Φ_tol=1E-8)

In-place version of [`USMEM_EOM`](@ref).
"""
function USMEM_EOM!(
    du::AbstractVector,
    u::AbstractVector,
    p::ComponentVector,
    t::Number,
    models::AbstractDynamicsModel;
    Φ_tol::Float64=1E-8,
)
    du .= USMEM_EOM(u, p, t, models; Φ_tol=Φ_tol)

    return nothing
end
