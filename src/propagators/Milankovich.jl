export Milankovich_EOM, Milankovich_EOM!
"""
    function Milankovich_EOM(
        u::AbstractVector,
        p::ComponentVector,
        t::Number,
        models::NTuple{N,AstroForceModels.AbstractAstroForceModel},
    ) where {N}

Milankovich propagation schema for orbital trajectories

Arguments:
-`u::AbstractVector`: The current Milankovich state.
-`p::ComponentVector`: The parameter vector, the simulation start date JD and the central body gravitational parameter.
-`t::Number`: The current time.
-`models::NTuple{N,AstroForceModels.AbstractAstroForceModel}`: Tuple of the acceleration models.

Returns:
-`du::AbstractVector`: Instantenous rate of change of the current state with respect to time.
"""
function Milankovich_EOM(
    u::AbstractVector,
    p::ComponentVector,
    t::Number,
    models::NTuple{N,AstroForceModels.AbstractAstroForceModel},
) where {N}
    Hx, Hy, Hz, _, _, _, _ = u
    μ = p.μ

    u_cart = AstroCoords.Mil2cart(u, μ)

    r = SVector{3}(u_cart[1], u_cart[2], u_cart[3])
    v = SVector{3}(u_cart[4], u_cart[5], u_cart[6])

    ad =
        build_dynamics_model(u_cart, p, t, models) -
        acceleration(u_cart, p, t, KeplerianGravityAstroModel(μ))

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
    function Milankovich_EOM!(
        du::AbstractVector,
        u::AbstractVector,
        p::ComponentVector,
        t::Number,
        models::NTuple{N,AstroForceModels.AbstractAstroForceModel},
    ) where {N}

Milankovich propagation schema for orbital trajectories

Arguments:
-`du::AbstractVector`: In-place vector to store the instantenous rate of change of the current state with respect to time.
-`u::AbstractVector`: The current Milankovich state.
-`p::ComponentVector`: The parameter vector, the simulation start date JD and the central body gravitational parameter.
-`t::Number`: The current time.
-`models::NTuple{N,AstroForceModels.AbstractAstroForceModel}`: Tuple of the acceleration models.

Returns:
- `nothing`
"""
function Milankovich_EOM!(
    du::AbstractVector,
    u::AbstractVector,
    p::ComponentVector,
    t::Number,
    models::NTuple{N,AstroForceModels.AbstractAstroForceModel},
) where {N}
    du .= Milankovich_EOM(u, p, t, models)

    return nothing
end
