module AstroPropagators

using AstroCoords
using AstroForceModels
using ComponentArrays
using LinearAlgebra
using OrdinaryDiffEqAdamsBashforthMoulton
using OrdinaryDiffEqCore
using SciMLBase
using StaticArraysCore

############################################################################################
#                                         Includes
############################################################################################
# Utilities
include("auxiliary/util.jl")

# Propagators
include("propagators/Cowell.jl")
include("propagators/EDromo.jl")
include("propagators/GaussVE.jl")
#include("propagators/Kustaanheimo-Stiefel.jl")
include("propagators/Milankovich.jl")
#include("propagators/Stiefel-Scheifel.jl")
include("propagators/USM.jl")

# Events
include("events/impulsive_maneuvers.jl")

include("api.jl")

end
