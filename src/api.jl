
# ========================================================================================
# Type Hierarchy
# ========================================================================================

abstract type AbstractPropagator end
abstract type AbstractStandardPropagator <: AbstractPropagator end
abstract type AbstractRegularizedPropagator <: AbstractPropagator end

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
# EOM Dispatch
#
# New propagators extend the API by defining:
#   eom!(::MyPropagator, du, u, p, t, models)           for standard
#   eom!(::MyPropagator, du, u, p, t, models, config)   for regularized
# ========================================================================================

export eom!

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
# propagate — Standard Propagators
# ========================================================================================

"""
    propagate(
        prop::AbstractStandardPropagator,
        u0::AbstractArray,
        p::ComponentArray,
        models::AbstractDynamicsModel,
        tspan::Tuple;
        solver=VCABM(),
        abstol=1e-13,
        reltol=1e-13,
        kwargs...,
    )

Propagate an orbit using a standard (non-regularized) formulation.

# Arguments
- `prop`: Propagator type (e.g., `CowellPropagator()`, `USM7Propagator()`).
- `u0`: Initial state vector in the propagator's coordinate system.
- `p::ComponentArray`: Parameter vector (must contain `μ` and `JD`).
- `models::AbstractDynamicsModel`: Force model composition.
- `tspan::Tuple`: Integration time span `(t0, tf)` in seconds.

# Keyword Arguments
- `solver`: ODE solver algorithm (default: `VCABM()`).
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
    solver::OrdinaryDiffEqCore.OrdinaryDiffEqAlgorithm=VCABM(),
    abstol::Real=1e-13,
    reltol::Real=1e-13,
    kwargs...,
)
    f!(du, u, p, t) = eom!(prop, du, u, p, t, models)
    prob = ODEProblem{true}(f!, u0, tspan, p)
    return solve(prob, solver; reltol, abstol, kwargs...)
end

# ========================================================================================
# propagate — Regularized Propagators
# ========================================================================================

"""
    propagate(
        prop::AbstractRegularizedPropagator,
        u0::AbstractArray,
        p::ComponentArray,
        models::AbstractDynamicsModel,
        tspan::Tuple,
        config::RegularizedCoordinateConfig;
        solver=VCABM(),
        abstol=1e-13,
        reltol=1e-13,
        kwargs...,
    )

Propagate an orbit using a regularized formulation that requires a
[`RegularizedCoordinateConfig`](@ref).

# Arguments
- `prop`: Propagator type (e.g., `EDromoPropagator()`, `GEqOEPropagator()`).
- `u0`: Initial state vector in the propagator's coordinate system.
- `p::ComponentArray`: Parameter vector (must contain `μ` and `JD`).
- `models::AbstractDynamicsModel`: Force model composition.
- `tspan::Tuple`: Integration time span `(t0, tf)` in seconds.
- `config::RegularizedCoordinateConfig`: Regularization configuration (perturbing
  potential, characteristic scales, time element type).

# Keyword Arguments
- `solver`: ODE solver algorithm (default: `VCABM()`).
- `abstol`: Absolute tolerance (default: `1e-13`).
- `reltol`: Relative tolerance (default: `1e-13`).
- `kwargs...`: Additional keyword arguments forwarded to `OrdinaryDiffEq.solve`.

# Returns
- `ODESolution` from DifferentialEquations.jl.
"""
function propagate(
    prop::AbstractRegularizedPropagator,
    u0::AbstractArray,
    p::ComponentArray,
    models::AstroForceModels.AbstractDynamicsModel,
    tspan::Tuple,
    config::RegularizedCoordinateConfig;
    solver::OrdinaryDiffEqCore.OrdinaryDiffEqAlgorithm=VCABM(),
    abstol::Real=1e-13,
    reltol::Real=1e-13,
    kwargs...,
)
    f!(du, u, p, t) = eom!(prop, du, u, p, t, models, config)
    prob = ODEProblem{true}(f!, u0, tspan, p)
    return solve(prob, solver; reltol, abstol, kwargs...)
end
