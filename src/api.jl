
abstract type AbstractPropType end

export CowellPropagator,
    GaussVEPropagator,
    MilankovichPropagator,
    USM7Propagator,
    USM6Propagator,
    USMEMPropagator

struct CowellPropagator <: AbstractPropType end
struct GaussVEPropagator <: AbstractPropType end
struct MilankovichPropagator <: AbstractPropType end
struct USM7Propagator <: AbstractPropType end
struct USM6Propagator <: AbstractPropType end
struct USMEMPropagator <: AbstractPropType end

function propagate(
    u0::AbstractArray{UT},
    p::ComponentArray,
    models::NTuple{N,AstroForceModels.AbstractAstroForceModel},
    tspan::Tuple{TT,TT};
    prop_type::AbstractPropType=CowellPropagator(),
    tsteps::Union{Vector{<:AbstractFloat},Nothing}=nothing,
    ODE_solver::OrdinaryDiffEqCore.OrdinaryDiffEqAlgorithm=VCABM(),
    abstol::Float64=1E-13,
    reltol::Float64=1E-13,
    output_file::Union{String,Nothing}=nothing,
) where {N,UT<:Number,TT<:Number}

    #TODO: DO THIS MORE INTELLIGENTLY
    EOM!(du, u, p, t) =
        if prop_type == CowellPropagator()
            Cowell_EOM!(du, u, p, t, models)
        elseif prop_type == GaussVEPropagator()
            GaussVE_EOM!(du, u, p, t, models)
        elseif prop_type == MilankovichPropagator()
            Milankovich_EOM!(du, u, p, t, models)
        elseif prop_type == USM7Propagator()
            USM7_EOM!(du, u, p, t, models)
        elseif prop_type == USM6Propagator()
            USM6_EOM!(du, u, p, t, models)
        elseif prop_type == USMEMPropagator()
            USMEM_EOM!(du, u, p, t, models)
        end

    ODE_prob::ODEProblem{UT,TT,true} = ODEProblem{true}(EOM!, u0, tspan, p)

    sol = solve(ODE_prob, ODE_solver; reltol=reltol, abstol=abstol)

    if output_file !== nothing
        times = (tsteps === nothing) ? sol.t : tsteps
        states = Array(sol(times))

        #TODO: print to file
        #TODO: Support Ephemeris Writing
    end

    return sol
end
