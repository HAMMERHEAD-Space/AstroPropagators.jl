export USM7_EOM, USM7_EOM!
"""
    function USM7_EOM(
        u::AbstractVector,
        p::ComponentVector,
        t::Number,
        models::NTuple{N,AstroForceModels.AbstractAstroForceModel},
    ) where {N}

Unified State Model (quaternions) propagation schema for orbital trajectories

Arguments:
-`u::AbstractVector`: The current USM7 state.
-`p::ComponentVector`: The parameter vector, the simulation start date JD and the central body gravitational parameter.
-`t::Number`: The current time.
-`models::NTuple{N,AstroForceModels.AbstractAstroForceModel}`: Tuple of the acceleration models.

Returns:
-`du::AbstractVector`: Instantenous rate of change of the current state with respect to time.
"""
function USM7_EOM(
    u::AbstractVector,
    p::ComponentVector,
    t::Number,
    models::NTuple{N,AstroForceModels.AbstractAstroForceModel},
) where {N}
    C, Rf1, Rf2, ϵO1, ϵO2, ϵO3, η0 = u

    μ::Number = p.μ

    sinλ = (2 * ϵO3 * η0) / (ϵO3^2 + η0^2)
    cosλ = (η0^2 - ϵO3^2) / (ϵO3^2 + η0^2)
    l = (ϵO1 * ϵO3 - ϵO2 * η0) / (ϵO3^2 + η0^2)
    ve2 = C - Rf1 * sinλ + Rf2 * cosλ
    ω3 = (C * ve2^2) / μ
    ρ = C / ve2

    u_cart = Cartesian(USM7(u), p.μ)

    fe =
        RTN_frame(u_cart) * (
            build_dynamics_model(u_cart, p, t, models) -
            acceleration(u_cart, p, t, KeplerianGravityAstroModel(; μ=μ))
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
    function USM7_EOM!(
        du::AbstractVector,
        u::AbstractVector,
        p::ComponentVector,
        t::Number,
        models::NTuple{N,AstroForceModels.AbstractAstroForceModel},
    ) where {N}

Unified State Model (quaternions) propagation schema for orbital trajectories

Arguments:
-`du::AbstractVector`: In-place vector to store the instantenous rate of change of the current state with respect to time.
-`u::AbstractVector`: The current USM7 state.
-`p::ComponentVector`: The parameter vector, the simulation start date JD and the central body gravitational parameter.
-`t::Number`: The current time.
-`models::NTuple{N,AstroForceModels.AbstractAstroForceModel}`: Tuple of the acceleration models.

Returns:
- `nothing`
"""
function USM7_EOM!(
    du::AbstractVector,
    u::AbstractVector,
    p::ComponentVector,
    t::Number,
    models::NTuple{N,AstroForceModels.AbstractAstroForceModel},
) where {N}
    du .= USM7_EOM(u, p, t, models)

    return nothing
end

export USM6_EOM, USM6_EOM!
"""
    function USM6_EOM(
        u::AbstractVector,
        p::ComponentVector,
        t::Number,
        models::NTuple{N,AstroForceModels.AbstractAstroForceModel},
    ) where {N}

Unified State Model (MRP's) propagation schema for orbital trajectories

Arguments:
-`u::AbstractVector`: The current USM6 state.
-`p::ComponentVector`: The parameter vector, the simulation start date JD and the central body gravitational parameter.
-`t::Number`: The current time.
-`models::NTuple{N,AstroForceModels.AbstractAstroForceModel}`: Tuple of the acceleration models.

Returns:
-`du::AbstractVector`: Instantenous rate of change of the current state with respect to time.
"""
function USM6_EOM(
    u::AbstractVector,
    p::ComponentVector,
    t::Number,
    models::NTuple{N,AstroForceModels.AbstractAstroForceModel},
) where {N}
    σ = SVector{3}(u[4], u[5], u[6])
    σ_norm = √(sum(abs2.(σ)))

    #TODO: SHADOW SET SHOULD PROBABLY BE EVENT
    if σ_norm > 1.0
        u[4:6] .= -σ ./ σ_norm
    end

    C, Rf1, Rf2, σ1, σ2, σ3 = u
    μ::Number = p.μ

    σ = SVector{3}(u[4], u[5], u[6])
    σ_norm = √(sum(abs2.(σ)))

    _, _, _, ϵO1, ϵO2, ϵO3, η0 = USM7(USM6(u), μ)
    u_cart = Cartesian(USM6(u), μ)

    sinλ = (4.0 * σ3 * (1.0 - σ_norm^2)) / (4.0 * σ3^2 + (1 - σ_norm^2)^2)
    cosλ = ((1.0 - σ_norm^2)^2 - 4.0 * σ3^2) / (4.0 * σ3^2 + (1 - σ_norm^2)^2)

    l = (ϵO1 * ϵO3 - ϵO2 * η0) / (ϵO3^2 + η0^2)

    ve2 = C - Rf1 * sinλ + Rf2 * cosλ

    ω3 = (C * ve2^2) / μ

    ρ = C / ve2

    fe =
        RTN_frame(u_cart) * (
            build_dynamics_model(u_cart, p, t, models) -
            acceleration(u_cart, p, t, KeplerianGravityAstroModel(; μ=μ))
        )

    ω1 = fe[3] / ve2

    return SVector{6}(
        -ρ * fe[2],
        fe[1] * cosλ - fe[2] * (1.0 + ρ) * sinλ - fe[3] * l * (Rf2 / ve2),
        fe[1] * sinλ + fe[2] * (1.0 + ρ) * cosλ + fe[3] * l * (Rf1 / ve2),
        0.25 * ((1.0 - σ_norm^2 + 2.0 * σ1^2) * ω1 + 2.0 * (σ1 * σ3 + σ2) * ω3),
        0.25 * (2.0 * (σ2 * σ1 + σ3) * ω1 + 2.0 * (σ2 * σ3 - σ1) * ω3),
        0.25 * (2.0 * (σ3 * σ1 - σ2) * ω1 + (1.0 - σ_norm^2 + 2.0 * σ3^2) * ω3),
    )
end

"""
    function USM6_EOM!(
        du::AbstractVector,
        u::AbstractVector,
        p::ComponentVector,
        t::Number,
        models::NTuple{N,AstroForceModels.AbstractAstroForceModel},
    ) where {N}

Unified State Model (MRP's) propagation schema for orbital trajectories

Arguments:
-`du::AbstractVector`: In-place vector to store the instantenous rate of change of the current state with respect to time.
-`u::AbstractVector`: The current USM6 state.
-`p::ComponentVector`: The parameter vector, the simulation start date JD and the central body gravitational parameter.
-`t::Number`: The current time.
-`models::NTuple{N,AstroForceModels.AbstractAstroForceModel}`: Tuple of the acceleration models.

Returns:
- `nothing`
"""
function USM6_EOM!(
    du::AbstractVector,
    u::AbstractVector,
    p::ComponentVector,
    t::Number,
    models::NTuple{N,AstroForceModels.AbstractAstroForceModel},
) where {N}
    du .= USM6_EOM(u, p, t, models)

    return nothing
end

export USMEM_EOM, USMEM_EOM!
"""
    function USMEM_EOM(
        u::AbstractVector,
        p::ComponentVector,
        t::Number,
        models::NTuple{N,AstroForceModels.AbstractAstroForceModel},
    ) where {N}

Unified State Model (exponential mapping) propagation schema for orbital trajectories

Arguments:
-`u::AbstractVector`: The current USMEM state.
-`p::ComponentVector`: The parameter vector, the simulation start date JD and the central body gravitational parameter.
-`t::Number`: The current time.
-`models::NTuple{N,AstroForceModels.AbstractAstroForceModel}`: Tuple of the acceleration models.

# Keyword Arguments"
-`Φ_tol::Float64`: The value to switch to the Taylor series expansion to avoid singularity.

Returns:
-`du::AbstractVector`: Instantenous rate of change of the current state with respect to time.
"""
function USMEM_EOM(
    u::AbstractVector,
    p::ComponentVector,
    t::Number,
    models::NTuple{N,AstroForceModels.AbstractAstroForceModel};
    Φ_tol::Float64=1E-8,
) where {N}
    C, Rf1, Rf2, a1, a2, a3 = u

    a = SVector{3}(a1, a2, a3)
    Φ = √(sum(abs2.(a)))
    a_cross = skew_sym(a)

    μ::Number = p.μ

    (_, _, _, ϵO1, ϵO2, ϵO3, η0) = USM7(USMEM(u), μ)
    u_cart = Cartesian(USMEM(u), μ)

    sinλ = (2 * ϵO3 * η0) / (ϵO3^2 + η0^2)
    cosλ = (η0^2 - ϵO3^2) / (ϵO3^2 + η0^2)
    l = (ϵO1 * ϵO3 - ϵO2 * η0) / (ϵO3^2 + η0^2)
    ve2 = C - Rf1 * sinλ + Rf2 * cosλ
    ω3 = (C * ve2^2) / μ
    ρ = C / ve2

    fe =
        RTN_frame(u_cart) * (
            build_dynamics_model(u_cart, p, t, models) -
            acceleration(u_cart, p, t, KeplerianGravityAstroModel(; μ=μ))
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
    function USMEM_EOM!(
        du::AbstractVector,
        u::AbstractVector,
        p::ComponentVector,
        t::Number,
        models::NTuple{N,AstroForceModels.AbstractAstroForceModel},
    ) where {N}

Unified State Model (exponential mapping) propagation schema for orbital trajectories

Arguments:
-`du::AbstractVector`: In-place vector to store the instantenous rate of change of the current state with respect to time.
-`u::AbstractVector`: The current USM7 state.
-`p::ComponentVector`: The parameter vector, the simulation start date JD and the central body gravitational parameter.
-`t::Number`: The current time.
-`models::NTuple{N,AstroForceModels.AbstractAstroForceModel}`: Tuple of the acceleration models.

# Keyword Arguments"
-`Φ_tol::Float64`: The value to switch to the Taylor series expansion to avoid singularity.

Returns:
- `nothing`
"""
function USMEM_EOM!(
    du::AbstractVector,
    u::AbstractVector,
    p::ComponentVector,
    t::Number,
    models::NTuple{N,AstroForceModels.AbstractAstroForceModel};
    Φ_tol::Float64=1E-8,
) where {N}
    du .= USMEM_EOM(u, p, t, models; Φ_tol=Φ_tol)

    return nothing
end
