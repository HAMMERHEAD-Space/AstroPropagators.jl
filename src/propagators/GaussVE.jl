export GaussVE_EOM, GaussVE_EOM!
"""
    function GaussVE_EOM(
        u::AbstractVector,
        p::ComponentVector,
        t::Number,
        models::NTuple{N,AstroForceModels.AbstractAstroForceModel},
    ) where {N}

Gauss variational propagation schema for orbital trajectories

Arguments:
-`u::AbstractVector`: The current Keplerian state.
-`p::ComponentVector`: The parameter vector, the simulation start date JD and the central body gravitational parameter.
-`t::Number`: The current time.
-`models::NTuple{N,AstroForceModels.AbstractAstroForceModel}`: Tuple of the acceleration models.

Returns:
-`du::AbstractVector`: Instantenous rate of change of the current state with respect to time.
"""
function GaussVE_EOM(
    u::AbstractVector,
    ps::ComponentVector,
    t::Number,
    models::NTuple{N,AstroForceModels.AbstractAstroForceModel},
) where {N}
    a, e, i, _, ω, f = u
    μ = ps.μ

    u_cart = Cartesian(Keplerian(u), μ)
    acc =
        RTN_frame(u_cart) * (
            build_dynamics_model(u_cart, ps, t, models) -
            acceleration(u_cart, ps, t, KeplerianGravityAstroModel(; μ=μ))
        )

    sf, cf = sincos(f)

    p = a * (1.0 - e^2)
    r = p / (1.0 + e * cf)
    θ = ω + f
    h = √(μ * p)

    sθ, cθ = sincos(θ)
    si, ci = sincos(i)

    ar = acc[1]
    aθ = acc[2] * sf
    ah = acc[3]

    da = (2.0 * a^2) / h * (e * sf * ar + p / r * aθ)
    de = 1.0 / h * (p * sf * ar + ((p + r) * cf + r * e) * aθ)
    di = (r * cθ / h) * ah
    dΩ = (r * sθ) / (h * si) * ah
    dω = 1.0 / (h * e) * (-p * cf * ar + (p + r) * sf * aθ) - (r * sθ * ci) / (h * si) * ah
    df = h / (r^2) + 1.0 / (h * e) * (p * cf * ar - (p + r) * sf * aθ)

    return SVector{6}(da, de, di, dΩ, dω, df)
end

"""
    function GaussVE_EOM!(
        du::AbstractVector,
        u::AbstractVector,
        p::ComponentVector,
        t::Number,
        models::NTuple{N,AstroForceModels.AbstractAstroForceModel},
    ) where {N}

Gauss variational propagation schema for orbital trajectories

Arguments:
-`du::AbstractVector`: In-place vector to store the instantenous rate of change of the current state with respect to time.
-`u::AbstractVector`: The current Keplerian state.
-`p::ComponentVector`: The parameter vector, the simulation start date JD and the central body gravitational parameter.
-`t::Number`: The current time.
-`models::NTuple{N,AstroForceModels.AbstractAstroForceModel}`: Tuple of the acceleration models.

Returns:
- `nothing`
"""
function GaussVE_EOM!(
    du::AbstractVector,
    u::AbstractVector,
    p::ComponentVector,
    t::Number,
    models::NTuple{N,AstroForceModels.AbstractAstroForceModel},
) where {N}
    du .= GaussVE_EOM(u, p, t, models)

    return nothing
end
