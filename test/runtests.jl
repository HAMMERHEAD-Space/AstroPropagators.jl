using Test

using AstroCoords
using AstroForceModels
using AstroPropagators
using ComponentArrays
using LinearAlgebra
using OrdinaryDiffEqCore
using OrdinaryDiffEqVerner
using SatelliteToolboxGravityModels
using SatelliteToolboxTransformations
using SciMLBase
using SpaceIndices

using Aqua
using JET
using AllocCheck

include("test_helpers.jl")

@testset "AstroPropagators.jl" begin
    include("propagators/test_cowell.jl")
    include("propagators/test_EDromo.jl")
    include("propagators/test_gaussVE.jl")
    include("propagators/test_KS.jl")
    include("propagators/test_milankovich.jl")
    include("propagators/test_StiSche.jl")
    include("propagators/test_USM.jl")
    include("propagators/test_GEqOE.jl")
    include("events/test_impulsive_maneuvers.jl")
    include("test_cross_formulation.jl")
    include("test_roundtrip.jl")
    include("test_api.jl")
end

# Differentiability tests are gated behind an environment variable to keep the default
# test suite fast. Set ASTROPROPAGATORS_TEST_DIFF to run them:
#   "true"  / "all"       → all 5 backends (ForwardDiff, Enzyme, Mooncake, PolyesterForwardDiff, Zygote)
#   "ForwardDiff"         → ForwardDiff only (fast smoke-test for AD compatibility)
#   unset   / "false"     → skip differentiability tests entirely
const _DIFF_ENV = get(ENV, "ASTROPROPAGATORS_TEST_DIFF", "false")

if _DIFF_ENV ∉ ("false", "")
    using DifferentiationInterface
    using FiniteDiff

    _run_all = _DIFF_ENV ∈ ("true", "all")
    _requested = _run_all ? Set{String}() : Set(strip.(split(_DIFF_ENV, ",")))
    _need(name) = _run_all || name ∈ _requested

    _backend_list = Tuple{String,Any}[]

    if _need("ForwardDiff")
        using ForwardDiff
        push!(_backend_list, ("ForwardDiff", AutoForwardDiff()))
    end
    if _need("Enzyme")
        using Enzyme
        push!(
            _backend_list,
            ("Enzyme", AutoEnzyme(; mode=Enzyme.set_runtime_activity(Enzyme.Forward))),
        )
    end
    if _need("Mooncake")
        using Mooncake
        push!(_backend_list, ("Mooncake", AutoMooncake(; config=nothing)))
    end
    if _need("PolyesterForwardDiff")
        using PolyesterForwardDiff
        push!(_backend_list, ("PolyesterForwardDiff", AutoPolyesterForwardDiff()))
    end
    if _need("Zygote")
        using Zygote
        push!(_backend_list, ("Zygote", AutoZygote()))
    end

    if isempty(_backend_list)
        error(
            "ASTROPROPAGATORS_TEST_DIFF=\"$_DIFF_ENV\" did not match any backend. " *
            "Valid names: ForwardDiff, Enzyme, Mooncake, PolyesterForwardDiff, Zygote",
        )
    end

    const _BACKENDS = Tuple(_backend_list)

    @info "Running differentiability tests with backends: $(join([b[1] for b in _BACKENDS], ", "))"

    @testset "Differentiability" begin
        include("differentiability/test_model_parameters.jl")
        include("differentiability/test_cowell.jl")
        include("differentiability/test_edromo.jl")
        include("differentiability/test_gaussVE.jl")
        include("differentiability/test_KS.jl")
        include("differentiability/test_milankovich.jl")
        include("differentiability/test_StiSche.jl")
        include("differentiability/test_USM.jl")
        include("differentiability/test_GEqOE.jl")
    end
else
    @info "Skipping differentiability tests (set ASTROPROPAGATORS_TEST_DIFF to enable)"
end

@testset "Code Performance" begin
    include("test_performance.jl")
end
