export propagate, propagate!, eom, eom!

export AbstractPropagator,
    AbstractStandardPropagator,
    AbstractRegularizedPropagator,
    CowellPropagator,
    EDromoPropagator,
    GaussVEPropagator,
    GEqOEPropagator,
    KSPropagator,
    MilankovichPropagator,
    StiSchePropagator,
    USM7Propagator,
    USM6Propagator,
    USMEMPropagator

# ========================================================================================
# Type Hierarchy
# ========================================================================================

abstract type AbstractPropagator end
abstract type AbstractStandardPropagator <: AbstractPropagator end
abstract type AbstractRegularizedPropagator <: AbstractPropagator end

struct CowellPropagator <: AbstractStandardPropagator end
struct GaussVEPropagator <: AbstractStandardPropagator end
struct MilankovichPropagator <: AbstractStandardPropagator end
struct USM7Propagator <: AbstractStandardPropagator end
struct USM6Propagator <: AbstractStandardPropagator end
struct USMEMPropagator <: AbstractStandardPropagator end

struct EDromoPropagator <: AbstractRegularizedPropagator end
struct KSPropagator <: AbstractRegularizedPropagator end
struct StiSchePropagator <: AbstractRegularizedPropagator end
struct GEqOEPropagator <: AbstractRegularizedPropagator end

# ========================================================================================
# EOM Dispatch — In-place (eom!)
# ========================================================================================

"""
    eom!(prop::AbstractStandardPropagator, du, u, p, t, models)
    eom!(prop::AbstractRegularizedPropagator, du, u, p, t, models, config)

In-place equations of motion. Computes derivatives and writes them into `du`.
Dispatches to the underlying EOM function for the given propagator type
(e.g., `eom!(CowellPropagator(), ...)` calls `Cowell_EOM!(...)`).
"""
@inline eom!(::CowellPropagator, du, u, p, t, models) = Cowell_EOM!(du, u, p, t, models)
@inline eom!(::GaussVEPropagator, du, u, p, t, models) = GaussVE_EOM!(du, u, p, t, models)
@inline eom!(::MilankovichPropagator, du, u, p, t, models) = Milankovich_EOM!(
    du, u, p, t, models
)
@inline eom!(::USM7Propagator, du, u, p, t, models) = USM7_EOM!(du, u, p, t, models)
@inline eom!(::USM6Propagator, du, u, p, t, models) = USM6_EOM!(du, u, p, t, models)
@inline eom!(::USMEMPropagator, du, u, p, t, models) = USMEM_EOM!(du, u, p, t, models)

@inline eom!(::EDromoPropagator, du, u, p, t, models, config) = EDromo_EOM!(
    du, u, p, t, models, config
)
@inline eom!(::KSPropagator, du, u, p, t, models, config) = KS_EOM!(
    du, u, p, t, models, config
)
@inline eom!(::StiSchePropagator, du, u, p, t, models, config) = StiSche_EOM!(
    du, u, p, t, models, config
)
@inline eom!(::GEqOEPropagator, du, u, p, t, models, config) = GEqOE_EOM!(
    du, u, p, t, models, config
)

# ========================================================================================
# EOM Dispatch — Out-of-place (eom)
# ========================================================================================

"""
    eom(prop::AbstractStandardPropagator, u, p, t, models)
    eom(prop::AbstractRegularizedPropagator, u, p, t, models, config)

Out-of-place equations of motion. Returns the derivative vector as an `SVector`.
Dispatches to the underlying EOM function for the given propagator type
(e.g., `eom(CowellPropagator(), ...)` calls `Cowell_EOM(...)`).
"""
@inline eom(::CowellPropagator, u, p, t, models) = Cowell_EOM(u, p, t, models)
@inline eom(::GaussVEPropagator, u, p, t, models) = GaussVE_EOM(u, p, t, models)
@inline eom(::MilankovichPropagator, u, p, t, models) = Milankovich_EOM(u, p, t, models)
@inline eom(::USM7Propagator, u, p, t, models) = USM7_EOM(u, p, t, models)
@inline eom(::USM6Propagator, u, p, t, models) = USM6_EOM(u, p, t, models)
@inline eom(::USMEMPropagator, u, p, t, models) = USMEM_EOM(u, p, t, models)

@inline eom(::EDromoPropagator, u, p, t, models, config) = EDromo_EOM(
    u, p, t, models, config
)
@inline eom(::KSPropagator, u, p, t, models, config) = KS_EOM(u, p, t, models, config)
@inline eom(::StiSchePropagator, u, p, t, models, config) = StiSche_EOM(
    u, p, t, models, config
)
@inline eom(::GEqOEPropagator, u, p, t, models, config) = GEqOE_EOM(u, p, t, models, config)

# ========================================================================================
# propagate — Out-of-place (ODEProblem{false})
# ========================================================================================

"""
    propagate(
        prop::AbstractStandardPropagator,
        u0::AbstractArray,
        p::ComponentArray,
        models::AbstractDynamicsModel,
        tspan::Tuple;
        solver=Vern9(),
        abstol=1e-13,
        reltol=1e-13,
        kwargs...,
    )

Propagate an orbit using a standard (non-regularized) formulation with an out-of-place
ODE problem. Suitable for use with `StaticArrays` or small state vectors.

See also [`propagate!`](@ref) for the in-place variant.

# Arguments
- `prop`: Propagator type (e.g., `CowellPropagator()`, `USM7Propagator()`).
- `u0`: Initial state vector in the propagator's coordinate system.
- `p::ComponentArray`: Parameter vector (must contain `μ` and `JD`).
- `models::AbstractDynamicsModel`: Force model composition.
- `tspan::Tuple`: Integration time span `(t0, tf)` in seconds.

# Keyword Arguments
- `solver`: ODE solver algorithm (default: `Vern9()`).
- `abstol`: Absolute tolerance (default: `1e-13`).
- `reltol`: Relative tolerance (default: `1e-13`).
- `kwargs...`: Additional keyword arguments forwarded to `OrdinaryDiffEq.solve`.

# Returns
- `ODESolution` from DifferentialEquations.jl.
"""
function propagate(
    prop::AbstractStandardPropagator,
    u0::AbstractArray,
    p::ComponentArray,
    models::AstroForceModels.AbstractDynamicsModel,
    tspan::Tuple;
    solver::OrdinaryDiffEqCore.OrdinaryDiffEqAlgorithm=Vern9(),
    abstol::Real=1e-13,
    reltol::Real=1e-13,
    kwargs...,
)
    u0_s = SVector{length(u0)}(u0)
    f(u, p, t) = eom(prop, u, p, t, models)
    prob = ODEProblem{false}(f, u0_s, tspan, p)
    return solve(prob, solver; reltol, abstol, kwargs...)
end

"""
    propagate(
        prop::AbstractRegularizedPropagator,
        u0::AbstractArray,
        p::ComponentArray,
        models::AbstractDynamicsModel,
        tspan::Tuple,
        config::RegularizedCoordinateConfig;
        solver=Vern9(),
        abstol=1e-13,
        reltol=1e-13,
        kwargs...,
    )

Out-of-place propagation for regularized formulations.

See also [`propagate!`](@ref) for the in-place variant.
"""
function propagate(
    prop::AbstractRegularizedPropagator,
    u0::AbstractArray,
    p::ComponentArray,
    models::AstroForceModels.AbstractDynamicsModel,
    tspan::Tuple,
    config::RegularizedCoordinateConfig;
    solver::OrdinaryDiffEqCore.OrdinaryDiffEqAlgorithm=Vern9(),
    abstol::Real=1e-13,
    reltol::Real=1e-13,
    kwargs...,
)
    u0_s = SVector{length(u0)}(u0)
    f(u, p, t) = eom(prop, u, p, t, models, config)
    prob = ODEProblem{false}(f, u0_s, tspan, p)
    return solve(prob, solver; reltol, abstol, kwargs...)
end

# ========================================================================================
# propagate! — In-place (ODEProblem{true})
# ========================================================================================

"""
    propagate!(
        prop::AbstractStandardPropagator,
        u0::AbstractArray,
        p::ComponentArray,
        models::AbstractDynamicsModel,
        tspan::Tuple;
        solver=Vern9(),
        abstol=1e-13,
        reltol=1e-13,
        kwargs...,
    )

Propagate an orbit using a standard (non-regularized) formulation with an in-place
ODE problem. Preferred for mutable state vectors and larger systems.

See also [`propagate`](@ref) for the out-of-place variant.

# Arguments
- `prop`: Propagator type (e.g., `CowellPropagator()`, `USM7Propagator()`).
- `u0`: Initial state vector in the propagator's coordinate system.
- `p::ComponentArray`: Parameter vector (must contain `μ` and `JD`).
- `models::AbstractDynamicsModel`: Force model composition.
- `tspan::Tuple`: Integration time span `(t0, tf)` in seconds.

# Keyword Arguments
- `solver`: ODE solver algorithm (default: `Vern9()`).
- `abstol`: Absolute tolerance (default: `1e-13`).
- `reltol`: Relative tolerance (default: `1e-13`).
- `kwargs...`: Additional keyword arguments forwarded to `OrdinaryDiffEq.solve`.

# Returns
- `ODESolution` from DifferentialEquations.jl.
"""
function propagate!(
    prop::AbstractStandardPropagator,
    u0::AbstractArray,
    p::ComponentArray,
    models::AstroForceModels.AbstractDynamicsModel,
    tspan::Tuple;
    solver::OrdinaryDiffEqCore.OrdinaryDiffEqAlgorithm=Vern9(),
    abstol::Real=1e-13,
    reltol::Real=1e-13,
    kwargs...,
)
    f!(du, u, p, t) = eom!(prop, du, u, p, t, models)
    prob = ODEProblem{true}(f!, u0, tspan, p)
    return solve(prob, solver; reltol, abstol, kwargs...)
end

"""
    propagate!(
        prop::AbstractRegularizedPropagator,
        u0::AbstractArray,
        p::ComponentArray,
        models::AbstractDynamicsModel,
        tspan::Tuple,
        config::RegularizedCoordinateConfig;
        solver=Vern9(),
        abstol=1e-13,
        reltol=1e-13,
        kwargs...,
    )

In-place propagation for regularized formulations.

See also [`propagate`](@ref) for the out-of-place variant.
"""
function propagate!(
    prop::AbstractRegularizedPropagator,
    u0::AbstractArray,
    p::ComponentArray,
    models::AstroForceModels.AbstractDynamicsModel,
    tspan::Tuple,
    config::RegularizedCoordinateConfig;
    solver::OrdinaryDiffEqCore.OrdinaryDiffEqAlgorithm=Vern9(),
    abstol::Real=1e-13,
    reltol::Real=1e-13,
    kwargs...,
)
    f!(du, u, p, t) = eom!(prop, du, u, p, t, models, config)
    prob = ODEProblem{true}(f!, u0, tspan, p)
    return solve(prob, solver; reltol, abstol, kwargs...)
end
