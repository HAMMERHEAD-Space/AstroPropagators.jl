export Milankovich_EOM, Milankovich_EOM!
"""
    Milankovich_EOM(u, p, t, models)

Milankovich propagation schema for orbital trajectories.

# Arguments
- `u::AbstractVector`: The current Milankovich state `[Hx, Hy, Hz, ex, ey, ez, L]`.
- `p::ComponentVector`: The parameter vector containing `μ` and `JD`.
- `t::Number`: The current time.
- `models::AbstractDynamicsModel`: Force model composition.

# Returns
- `SVector{7}`: Instantaneous rate of change of the Milankovich state.
"""
function Milankovich_EOM(
    u::AbstractVector, p::ComponentVector, t::Number, models::AbstractDynamicsModel
)
    Hx, Hy, Hz, _, _, _, _ = u
    μ::Number = p.μ

    u_cart = Cartesian(Milankovich(u), μ)

    r = SVector{3}(u_cart[1], u_cart[2], u_cart[3])
    v = SVector{3}(u_cart[4], u_cart[5], u_cart[6])

    ad =
        build_dynamics_model(u_cart, p, t, models) -
        acceleration(u_cart, p, t, KeplerianGravityAstroModel(; μ=μ))

    H = SVector{3}(Hx, Hy, Hz)
    ẑ = SVector{3}(0.0, 0.0, 1.0)

    dH = skew_sym(r) * ad
    de = 1 / μ * (skew_sym(v) * skew_sym(r) - skew_sym(H)) * ad
    dL =
        (dot(ẑ, r) / (norm(H) * (norm(H) + dot(ẑ, H)))) * dot(H, ad) +
        norm(H) / (norm(r)^2)

    return SVector{7}(dH[1], dH[2], dH[3], de[1], de[2], de[3], dL)
end

"""
    Milankovich_EOM!(du, u, p, t, models)

In-place version of [`Milankovich_EOM`](@ref).
"""
function Milankovich_EOM!(
    du::AbstractVector,
    u::AbstractVector,
    p::ComponentVector,
    t::Number,
    models::AbstractDynamicsModel,
)
    du .= Milankovich_EOM(u, p, t, models)

    return nothing
end
