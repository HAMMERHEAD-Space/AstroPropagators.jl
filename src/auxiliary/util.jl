function skew_sym(x::AbstractVector{<:Number})
    return SMatrix{3,3}(0.0, x[3], -x[2], -x[3], 0.0, x[1], x[2], -x[1], 0.0)
end

function quaternions2DCM(u::AbstractVector{<:Number})
    q0, q1, q2, q3 = u

    return [
        q0^2 + q1^2 - q2^2-q3^2 2*(q1 * q2 + q0 * q3) 2*(q1 * q3 - q0 * q2)
        2*(q1 * q2 - q0 * q3) q0^2 - q1^2 + q2^2-q3^2 2*(q2 * q3 + q0 * q1)
        2*(q1 * q3 + q0 * q2) 2*(q2 * q3 - q0 * q1) q0^2 - q1^2 - q2^2+q3^2
    ]
end

function RTN_frame(u::AbstractVector{<:Number})
    R = SVector{3}(u[1], u[2], u[3])
    V = SVector{3}(u[4], u[5], u[6])

    R̂ = R ./ √(sum(abs2.(R)))
    V̂ = V ./ √(sum(abs2.(V)))

    N = cross(R̂, V̂)
    N̂ = N ./ √(sum(abs2.(N)))
    T = cross(N, R)
    T̂ = T ./ √(sum(abs2.(T)))

    return [R̂'; T̂'; N̂']
end
