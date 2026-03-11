export inertial_to_RTN

"""
    skew_sym(x::AbstractVector) -> SMatrix{3,3}

Construct the 3×3 skew-symmetric (cross-product) matrix of a 3-vector `x`,
such that `skew_sym(x) * y == cross(x, y)` for any 3-vector `y`.
"""
function skew_sym(x::AbstractVector{<:Number})
    z = zero(eltype(x))
    return SMatrix{3,3}(z, x[3], -x[2], -x[3], z, x[1], x[2], -x[1], z)
end

"""
    inertial_to_RTN(a_inertial::AbstractVector, u::AbstractVector) -> SVector{3}

Project an inertial-frame vector into the RTN (Radial–Transverse–Normal) frame.

The RTN basis is constructed from a Cartesian state `u = [r; v]`:
- **R̂**: `r / |r|`
- **N̂**: `(r × v) / |r × v|`
- **T̂**: `N̂ × R̂`

Returns `SVector{3}(a⋅R̂, a⋅T̂, a⋅N̂)`.
"""
@inline function inertial_to_RTN(
    a_inertial::AbstractVector{AT}, u::AbstractVector{UT}
) where {AT,UT}
    r = SVector{3}(u[1], u[2], u[3])
    v = SVector{3}(u[4], u[5], u[6])

    r_norm = norm(r)
    R̂ = SVector{3}(r[1] / r_norm, r[2] / r_norm, r[3] / r_norm)

    h = cross(r, v)
    h_norm = norm(h)
    N̂ = SVector{3}(h[1] / h_norm, h[2] / h_norm, h[3] / h_norm)

    T̂ = cross(N̂, R̂)

    a = SVector{3,AT}(a_inertial[1], a_inertial[2], a_inertial[3])
    return SVector{3}(dot(a, R̂), dot(a, T̂), dot(a, N̂))
end
