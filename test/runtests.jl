using AstroCoords
using AstroForceModels
using AstroPropagators
using Aqua
using ComponentArrays
using LinearAlgebra
using OrdinaryDiffEqAdamsBashforthMoulton
using OrdinaryDiffEqCore
using SatelliteToolboxGravityModels
using SatelliteToolboxTransformations
using SciMLBase
using SpaceIndices
using Test

#using JET
#using AllocCheck

@testset "AstroPropagators.jl" begin
    include("propagators/test_cowell.jl")
    include("propagators/test_gaussVE.jl")
    include("propagators/test_milankovich.jl")
    include("propagators/test_USM.jl")
    include("events/test_impulsive_maneuvers.jl")
end

#TODO: NEED TO FIX IN AstroForceModels & SatelliteToolbox
#@testset "Code Performance" begin
#    include("test_JET.jl")
#    include("test_allocs.jl")
#end

@testset "Aqua.jl" begin
    Aqua.test_all(AstroPropagators; ambiguities=(recursive = false))
end
