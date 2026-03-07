export GaussVE_EOM, GaussVE_EOM!
"""
    GaussVE_EOM(
        u::AbstractVector,
        p::ComponentVector,
        t::Number,
        models::AbstractDynamicsModel,
    )

Gauss variational equations propagation schema for orbital trajectories.

# Arguments
- `u::AbstractVector`: The current Keplerian state `[a, e, i, Ω, ω, f]`.
- `p::ComponentVector`: The parameter vector containing `μ` and `JD`.
- `t::Number`: The current time.
- `models::AbstractDynamicsModel`: Force model composition.

# Returns
- `SVector{6}`: Instantaneous rate of change of the Keplerian state.
"""
function GaussVE_EOM(
    u::AbstractVector, ps::ComponentVector, t::Number, models::AbstractDynamicsModel
)
    a, e, i, _, ω, f = u
    μ::Number = ps.μ

    u_cart = Cartesian(Keplerian(u), μ)
    acc = inertial_to_RTN(
        build_dynamics_model(u_cart, ps, t, models) -
        acceleration(u_cart, ps, t, KeplerianGravityAstroModel(; μ=μ)),
        u_cart,
    )

    sf, cf = sincos(f)

    p = a * (1.0 - e^2)
    r = p / (1.0 + e * cf)
    θ = ω + f
    h = √(μ * p)

    sθ, cθ = sincos(θ)
    si, ci = sincos(i)

    ar = acc[1]
    aθ = acc[2]
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
    GaussVE_EOM!(du, u, p, t, models)

In-place version of [`GaussVE_EOM`](@ref).
"""
function GaussVE_EOM!(
    du::AbstractVector,
    u::AbstractVector,
    p::ComponentVector,
    t::Number,
    models::AbstractDynamicsModel,
)
    du .= GaussVE_EOM(u, p, t, models)

    return nothing
end
