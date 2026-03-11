export modified_equinoctial_gve, ModEq_EOM, ModEq_EOM!

"""
    modified_equinoctial_gve(p, f, g, h, k, L, μ) -> SMatrix{6,3}

Gauss Variational Equations matrix for Modified Equinoctial Elements.

Returns the 6×3 matrix `A` such that `d[p,f,g,h,k,L]/dt = A * [Fr,Fθ,Fh]`, where
`[Fr,Fθ,Fh]` is the perturbing acceleration in the RTN frame.

The first row gives dp/dF (semi-latus rectum). Rows 2–6 give
df/dF, dg/dF, dh/dF, dk/dF, dL/dF respectively. The Keplerian
contribution to dL/dt (= √(μp) q²/p²) is **not** included; callers
must add it separately.

# References
- Walker, M.J.H., Ireland, B. & Owens, J. (1985). "A Set of Modified Equinoctial Orbit
  Elements." Celestial Mechanics, 36, 409–419.
- Battin, R.H. (1999). An Introduction to the Mathematics and Methods of Astrodynamics.
  AIAA Education Series, Eq. 9.92.
"""
@inline function modified_equinoctial_gve(
    p::Number, f::Number, g::Number, h::Number, k::Number, L::Number, μ::Number
)
    T = promote_type(
        typeof(p), typeof(f), typeof(g), typeof(h), typeof(k), typeof(L), typeof(μ)
    )

    sL, cL = sincos(L)
    q = one(T) + f * cL + g * sL
    sqrt_pmu = sqrt(T(p) / T(μ))
    c = sqrt_pmu / q   # common factor  √(p/μ) / q

    # Row 1: dp/dF  —  dp/dt = 2p/q · √(p/μ) · Fθ
    A11 = zero(T)
    A12 = T(2) * T(p) * c
    A13 = zero(T)

    # Row 2: df/dF
    A21 = c * q * sL
    A22 = c * ((q + one(T)) * cL + T(f))
    A23 = -c * T(g) * (T(h) * sL - T(k) * cL)

    # Row 3: dg/dF
    A31 = -c * q * cL
    A32 = c * ((q + one(T)) * sL + T(g))
    A33 = c * T(f) * (T(h) * sL - T(k) * cL)

    # Row 4: dh/dF
    s_sq = one(T) + T(h)^2 + T(k)^2
    A41 = zero(T)
    A42 = zero(T)
    A43 = c * s_sq * cL / T(2)

    # Row 5: dk/dF
    A51 = zero(T)
    A52 = zero(T)
    A53 = c * s_sq * sL / T(2)

    # Row 6: dL/dF (perturbation only; Keplerian part added by caller)
    A61 = zero(T)
    A62 = zero(T)
    A63 = c * (T(h) * sL - T(k) * cL)

    return SMatrix{6,3,T}(
        A11,
        A21,
        A31,
        A41,
        A51,
        A61,
        A12,
        A22,
        A32,
        A42,
        A52,
        A62,
        A13,
        A23,
        A33,
        A43,
        A53,
        A63,
    )
end

"""
    ModEq_EOM(
        u::AbstractVector,
        p::ComponentVector,
        t::Number,
        models::AbstractDynamicsModel,
    )

Modified Equinoctial Elements propagation using Gauss Variational Equations.

# Arguments
- `u::AbstractVector`: The current MEE state `[p, f, g, h, k, L]`.
- `p::ComponentVector`: The parameter vector containing `μ` and `JD`.
- `t::Number`: The current time.
- `models::AbstractDynamicsModel`: Force model composition.

# Returns
- `SVector{6}`: Instantaneous rate of change of the MEE state.
"""
function ModEq_EOM(
    u::AbstractVector, ps::ComponentVector, t::Number, models::AbstractDynamicsModel
)
    p_el, f, g, h, k, L = u
    μ::Number = ps.μ

    u_cart = Cartesian(ModEq(u), μ)
    acc = inertial_to_RTN(
        build_dynamics_model(u_cart, ps, t, models) -
        acceleration(u_cart, ps, t, KeplerianGravityAstroModel(; μ=μ)),
        u_cart,
    )

    A = modified_equinoctial_gve(p_el, f, g, h, k, L, μ)
    oe_dot = A * acc

    # Keplerian true-longitude rate:  dL/dt = q² √(μp) / p²
    q = one(eltype(u)) + f * cos(L) + g * sin(L)
    dL_kep = q^2 * sqrt(μ * p_el) / p_el^2

    return SVector{6}(
        oe_dot[1], oe_dot[2], oe_dot[3], oe_dot[4], oe_dot[5], oe_dot[6] + dL_kep
    )
end

"""
    ModEq_EOM!(du, u, p, t, models)

In-place version of [`ModEq_EOM`](@ref).
"""
function ModEq_EOM!(
    du::AbstractVector,
    u::AbstractVector,
    p::ComponentVector,
    t::Number,
    models::AbstractDynamicsModel,
)
    du .= ModEq_EOM(u, p, t, models)

    return nothing
end
