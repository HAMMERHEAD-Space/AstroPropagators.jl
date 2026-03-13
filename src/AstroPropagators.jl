module AstroPropagators

using AstroCoords
using AstroForceModels
using ComponentArrays
using LinearAlgebra
using OrdinaryDiffEqCore
using OrdinaryDiffEqVerner
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
include("propagators/Kustaanheimo-Stiefel.jl")
include("propagators/Milankovich.jl")
include("propagators/Stiefel-Scheifele.jl")
include("propagators/USM.jl")
include("propagators/GEqOE.jl")
include("propagators/ModifiedEquinoctial.jl")

# API (type hierarchy, EOM dispatch, propagate)
include("api.jl")

# Events — state adapter (must come after api.jl for propagator types)
include("events/state_adapter.jl")

# Events — detectors
include("events/orbital_detectors.jl")
include("events/eclipse_detectors.jl")
include("events/geometric_detectors.jl")
include("events/utility_detectors.jl")

# Events — maneuver system
include("events/maneuver_types.jl")
include("events/maneuver_builder.jl")

# Events — legacy impulsive burn (kept for low-level use)
include("events/impulsive_maneuvers.jl")

end
