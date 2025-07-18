using Test

using AstroCoords
using AstroForceModels
using AstroPropagators
using ComponentArrays
using LinearAlgebra
using OrdinaryDiffEqAdamsBashforthMoulton
using OrdinaryDiffEqCore
using SatelliteToolboxGravityModels
using SatelliteToolboxTransformations
using SciMLBase
using SpaceIndices

using Aqua
using JET
using AllocCheck

using DifferentiationInterface
using FiniteDiff, ForwardDiff, Enzyme, Mooncake, PolyesterForwardDiff, Zygote

@testset "AstroPropagators.jl" begin
    include("propagators/test_cowell.jl")
    include("propagators/test_gaussVE.jl")
    include("propagators/test_milankovich.jl")
    include("propagators/test_USM.jl")
    include("events/test_impulsive_maneuvers.jl")
end

const _BACKENDS = (
    ("ForwardDiff", AutoForwardDiff()),
    ("Enzyme", AutoEnzyme()),
    ("Mooncake", AutoMooncake(; config=nothing)),
    ("PolyesterForwardDiff", AutoPolyesterForwardDiff()),
    ("Zygote", AutoZygote()),
)

#@testset "Differentiability" begin
#    include("differentiability/test_model_parameters.jl")
#    include("differentiability/test_cowell.jl")
#    include("differentiability/test_gaussVE.jl")
#    include("differentiability/test_milankovich.jl")
#    include("differentiability/test_USM.jl")
#end

@testset "Code Performance" begin
    include("test_performance.jl")
end
