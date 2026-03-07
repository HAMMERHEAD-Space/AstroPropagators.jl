export coord_type, get_cartesian, get_keplerian, get_physical_time

# ========================================================================================
# State Adapter Layer
#
# Conversion between raw integrator state vectors and AstroCoords.jl types.
# All event detectors and maneuver builders use these functions internally.
# Dispatches on Type{<:AstroCoord}, not propagator types.
# ========================================================================================

"""
    coord_type(::AbstractPropagator) -> Type{<:AstroCoord}

Return the native `AstroCoord` type for a given propagator. Bridge between the
`propagate()` API (which uses propagator types) and the event system (which uses
`AstroCoord` types).
"""
coord_type(::CowellPropagator) = Cartesian
coord_type(::GaussVEPropagator) = Keplerian
coord_type(::MilankovichPropagator) = Milankovich
coord_type(::USM7Propagator) = USM7
coord_type(::USM6Propagator) = USM6
coord_type(::USMEMPropagator) = USMEM
coord_type(::EDromoPropagator) = EDromo
coord_type(::KSPropagator) = KustaanheimoStiefel
coord_type(::StiSchePropagator) = StiefelScheifele
coord_type(::GEqOEPropagator) = GEqOE

# ========================================================================================
# get_cartesian — convert integrator state to Cartesian
# ========================================================================================

"""
    get_cartesian(u, t, μ, ::Type{C}) where {C<:AstroCoord}
    get_cartesian(u, t, μ, ::Type{C}, config::RegularizedCoordinateConfig) where {C<:AstroCoord}

Convert a raw integrator state vector `u` to `Cartesian` via the specified
coordinate type `C` and AstroCoords.jl's universal conversion system. Regularized
coordinate sets require a `RegularizedCoordinateConfig`.
"""
get_cartesian(u, t, μ, ::Type{Cartesian}) = Cartesian(u)

function get_cartesian(u, t, μ, ::Type{C}) where {C<:AstroCoords.AstroCoord}
    return Cartesian(C(u), μ)
end

function get_cartesian(u, t, μ, ::Type{EDromo}, config::RegularizedCoordinateConfig)
    return Cartesian(EDromo(u), μ, t, config)
end

function get_cartesian(
    u, t, μ, ::Type{KustaanheimoStiefel}, config::RegularizedCoordinateConfig
)
    return Cartesian(KustaanheimoStiefel(u), μ, config)
end

function get_cartesian(
    u, t, μ, ::Type{StiefelScheifele}, config::RegularizedCoordinateConfig
)
    return Cartesian(StiefelScheifele(u), μ, t, config)
end

function get_cartesian(u, t, μ, ::Type{GEqOE}, config::RegularizedCoordinateConfig)
    return Cartesian(GEqOE(u), μ, config)
end

# ========================================================================================
# get_keplerian — convert integrator state to Keplerian
# ========================================================================================

"""
    get_keplerian(u, t, μ, ::Type{C}) where {C<:AstroCoord}
    get_keplerian(u, t, μ, ::Type{C}, config) where {C<:AstroCoord}

Convert integrator state to `Keplerian` via `get_cartesian` then AstroCoords
conversion. The returned `Keplerian` supports named field access (`.a`, `.e`,
`.i`, `.Ω`, `.ω`, `.f`, `.M`, `.E`).
"""
function get_keplerian(u, t, μ, ::Type{C}) where {C<:AstroCoords.AstroCoord}
    cart = get_cartesian(u, t, μ, C)
    return Keplerian(cart, μ)
end

function get_keplerian(
    u, t, μ, ::Type{C}, config::RegularizedCoordinateConfig
) where {C<:AstroCoords.AstroCoord}
    cart = get_cartesian(u, t, μ, C, config)
    return Keplerian(cart, μ)
end

# ========================================================================================
# get_physical_time — recover physical time from the independent variable
# ========================================================================================

"""
    get_physical_time(u, t, ::Type{C}) where {C<:AstroCoord}
    get_physical_time(u, t, ::Type{C}, config) where {C<:AstroCoord}

Recover physical time [s] from the integrator's independent variable. For
standard coordinate sets, returns `t` unchanged. For regularized sets (EDromo,
KS, StiefelScheifele), extracts physical time from the state vector via the
time element. GEqOE is regularized but uses physical time directly.
"""
get_physical_time(u, t, ::Type{C}) where {C<:AstroCoords.AstroCoord} = t

get_physical_time(u, t, ::Type{GEqOE}, ::RegularizedCoordinateConfig) = t

function get_physical_time(u, t, ::Type{EDromo}, config::RegularizedCoordinateConfig)
    return get_EDromo_time(u, t, config)
end

function get_physical_time(
    u, t, ::Type{KustaanheimoStiefel}, config::RegularizedCoordinateConfig
)
    return get_KS_time(u, config)
end

function get_physical_time(
    u, t, ::Type{StiefelScheifele}, config::RegularizedCoordinateConfig
)
    return get_stiefelscheifele_time(u, t, config)
end
