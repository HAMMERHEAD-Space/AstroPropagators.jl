export Cowell_EOM, Cowell_EOM!
"""
    function Cowell_EOM(
        u::AbstractVector,
        p::ComponentVector,
        t::Number,
        models::NTuple{N,AstroForceModels.AbstractAstroForceModel},
    ) where {N}

Cowell propagation schema for orbital trajectories

Arguments:
-`u::AbstractVector`: The current Cartesian state.
-`p::ComponentVector`: The parameter vector, only the simulation start date JD is provided.
-`t::Number`: The current time.
-`models::NTuple{N,AstroForceModels.AbstractAstroForceModel}`: Tuple of the acceleration models.

Returns:
-`du::AbstractVector`: Instantenous rate of change of the current state with respect to time.
"""
function Cowell_EOM(
    u::AbstractVector,
    p::ComponentVector,
    t::Number,
    models::NTuple{N,AstroForceModels.AbstractAstroForceModel},
) where {N}
    accel = build_dynamics_model(u, p, t, models)
    return SVector{6}(u[4], u[5], u[6], accel[1], accel[2], accel[3])
end

"""
    function Cowell_EOM!(
        du::AbstractVector,
        u::AbstractVector,
        p::ComponentVector,
        t::Number,
        models::NTuple{N,AstroForceModels.AbstractAstroForceModel},
    ) where {N}

Cowell propagation schema for orbital trajectories

Arguments:
-`du::AbstractVector`: In-place vector to store the instantenous rate of change of the current state with respect to time.
-`u::AbstractVector`: The current Cartesian state.
-`p::ComponentVector`: The parameter vector, only the simulation start date JD is provided.
-`t::Number`: The current time.
-`models::NTuple{N,AstroForceModels.AbstractAstroForceModel}`: Tuple of the acceleration models.

Returns:
- `nothing`
"""
function Cowell_EOM!(
    du::AbstractVector,
    u::AbstractVector,
    p::ComponentVector,
    t::Number,
    models::NTuple{N,AstroForceModels.AbstractAstroForceModel},
) where {N}
    du .= Cowell_EOM(u, p, t, models)

    return nothing
end
