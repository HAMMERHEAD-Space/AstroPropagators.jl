export Cowell_EOM, Cowell_EOM!
"""
    Cowell_EOM(
        u::AbstractVector,
        p::ComponentVector,
        t::Number,
        models::AbstractDynamicsModel,
    )

Cowell propagation schema for orbital trajectories.

# Arguments
- `u::AbstractVector`: The current Cartesian state.
- `p::ComponentVector`: The parameter vector containing `μ` and `JD`.
- `t::Number`: The current time.
- `models::AbstractDynamicsModel`: Force model composition.

# Returns
- `SVector{6}`: Instantaneous rate of change of the Cartesian state.
"""
function Cowell_EOM(
    u::AbstractVector, p::ComponentVector, t::Number, models::AbstractDynamicsModel
)
    accel = build_dynamics_model(u, p, t, models)
    return SVector{6}(u[4], u[5], u[6], accel[1], accel[2], accel[3])
end

"""
    Cowell_EOM!(du, u, p, t, models)

In-place version of [`Cowell_EOM`](@ref).
"""
function Cowell_EOM!(
    du::AbstractVector,
    u::AbstractVector,
    p::ComponentVector,
    t::Number,
    models::AbstractDynamicsModel,
)
    du .= Cowell_EOM(u, p, t, models)

    return nothing
end
