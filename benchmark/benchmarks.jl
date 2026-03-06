using AstroCoords
using AstroForceModels
using AstroPropagators
using BenchmarkTools
using ComponentArrays
using LinearAlgebra
using OrdinaryDiffEqAdamsBashforthMoulton
using SatelliteToolboxGravityModels
using SatelliteToolboxTransformations
using SpaceIndices
using StaticArraysCore

const SUITE = BenchmarkGroup()

SUITE["eom"] = BenchmarkGroup(["equations of motion"])
SUITE["eom"]["standard"] = BenchmarkGroup(["non-regularized"])
SUITE["eom"]["regularized"] = BenchmarkGroup(["regularized coordinates"])
SUITE["propagation"] = BenchmarkGroup(["full orbit propagation"])

# ---------------------
# Common state and parameters
# ---------------------
const _JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
const _t = 0.0

const _state_cart = [
    -1076.225324679696
    -6765.896364327722
    -332.3087833503755
    9.356857417032581
    -3.3123476319597557
    -1.1880157328553503
]

SpaceIndices.init()
const _eop_data = fetch_iers_eop()
const _grav_coeffs = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96))

# ---------------------
# Keplerian force model
# ---------------------
const _kep_model = KeplerianGravityAstroModel()
const _μ_kep = _kep_model.μ
const _p_kep = ComponentVector(; JD=_JD, μ=_μ_kep)
const _dynamics_kep = CentralBodyDynamicsModel(_kep_model)

# ---------------------
# Full force model (J2 36x36 + Sun + Moon + SRP + Drag)
# ---------------------
const _grav_model = GravityHarmonicsAstroModel(;
    gravity_model=_grav_coeffs, eop_data=_eop_data, order=36, degree=36
)
const _μ_full = GravityModels.gravity_constant(_grav_model.gravity_model) / 1E9
const _p_full = ComponentVector(; JD=_JD, μ=_μ_full)

const _sun_model = ThirdBodyModel(; body=SunBody(), eop_data=_eop_data)
const _moon_model = ThirdBodyModel(; body=MoonBody(), eop_data=_eop_data)
const _srp_model = SRPAstroModel(;
    satellite_srp_model=CannonballFixedSRP(0.2),
    sun_data=_sun_model,
    eop_data=_eop_data,
    shadow_model=Conical(),
)
const _drag_model = DragAstroModel(;
    satellite_drag_model=CannonballFixedDrag(0.2),
    atmosphere_model=JB2008(),
    eop_data=_eop_data,
)
const _dynamics_full = CentralBodyDynamicsModel(
    _grav_model, (_sun_model, _moon_model, _srp_model, _drag_model)
)

# ---------------------
# Standard propagator states (non-regularized)
# ---------------------
const _state_koe = Array(Keplerian(Cartesian(_state_cart), _μ_kep))
const _state_mil = Array(Milankovich(Cartesian(_state_cart), _μ_kep))
const _state_usm7 = Array(USM7(Cartesian(_state_cart), _μ_kep))
const _state_usm6 = Array(USM6(Cartesian(_state_cart), _μ_kep))
const _state_usmem = Array(USMEM(Cartesian(_state_cart), _μ_kep))

# ---------------------
# Regularized propagator states
# ---------------------
const _config_pt = RegularizedCoordinateConfig(
    _state_cart, _μ_kep; W=-1e-2, t₀=0.0, flag_time=PhysicalTime()
)
const _ϕ_pt = AstroCoords.compute_initial_phi(_state_cart, _μ_kep, _config_pt)

const _state_edromo = Array(EDromo(Cartesian(_state_cart), _μ_kep, _ϕ_pt, _config_pt))
const _state_ks = Array(KustaanheimoStiefel(Cartesian(_state_cart), _μ_kep, _config_pt))
const _state_stische = Array(
    StiefelScheifele(Cartesian(_state_cart), _μ_kep, _ϕ_pt, _config_pt)
)

const _config_geqoe = RegularizedCoordinateConfig(; W=-1e-2)
const _state_geqoe = Array(GEqOE(Cartesian(_state_cart), _μ_kep, _config_geqoe))

# =====================
# Standard EOM benchmarks
# =====================

# --- Cowell ---
SUITE["eom"]["standard"]["Cowell (Keplerian)"] = @benchmarkable Cowell_EOM(
    $_state_cart, $_p_kep, $_t, $_dynamics_kep
)
SUITE["eom"]["standard"]["Cowell (Full)"] = @benchmarkable Cowell_EOM(
    $_state_cart, $_p_full, $_t, $_dynamics_full
)

# --- GaussVE ---
SUITE["eom"]["standard"]["GaussVE (Keplerian)"] = @benchmarkable GaussVE_EOM(
    $_state_koe, $_p_kep, $_t, $_dynamics_kep
)
SUITE["eom"]["standard"]["GaussVE (Full)"] = @benchmarkable GaussVE_EOM(
    $_state_koe, $_p_full, $_t, $_dynamics_full
)

# --- Milankovich ---
SUITE["eom"]["standard"]["Milankovich (Keplerian)"] = @benchmarkable Milankovich_EOM(
    $_state_mil, $_p_kep, $_t, $_dynamics_kep
)
SUITE["eom"]["standard"]["Milankovich (Full)"] = @benchmarkable Milankovich_EOM(
    $_state_mil, $_p_full, $_t, $_dynamics_full
)

# --- USM7 ---
SUITE["eom"]["standard"]["USM7 (Keplerian)"] = @benchmarkable USM7_EOM(
    $_state_usm7, $_p_kep, $_t, $_dynamics_kep
)
SUITE["eom"]["standard"]["USM7 (Full)"] = @benchmarkable USM7_EOM(
    $_state_usm7, $_p_full, $_t, $_dynamics_full
)

# --- USM6 ---
SUITE["eom"]["standard"]["USM6 (Keplerian)"] = @benchmarkable USM6_EOM(
    $_state_usm6, $_p_kep, $_t, $_dynamics_kep
)
SUITE["eom"]["standard"]["USM6 (Full)"] = @benchmarkable USM6_EOM(
    $_state_usm6, $_p_full, $_t, $_dynamics_full
)

# --- USMEM ---
SUITE["eom"]["standard"]["USMEM (Keplerian)"] = @benchmarkable USMEM_EOM(
    $_state_usmem, $_p_kep, $_t, $_dynamics_kep
)
SUITE["eom"]["standard"]["USMEM (Full)"] = @benchmarkable USMEM_EOM(
    $_state_usmem, $_p_full, $_t, $_dynamics_full
)

# =====================
# Regularized EOM benchmarks
# =====================

# --- EDromo ---
SUITE["eom"]["regularized"]["EDromo (Keplerian)"] = @benchmarkable EDromo_EOM(
    $_state_edromo, $_p_kep, $_ϕ_pt, $_dynamics_kep, $_config_pt
)
SUITE["eom"]["regularized"]["EDromo (Full)"] = @benchmarkable EDromo_EOM(
    $_state_edromo, $_p_full, $_ϕ_pt, $_dynamics_full, $_config_pt
)

# --- Kustaanheimo-Stiefel ---
SUITE["eom"]["regularized"]["KS (Keplerian)"] = @benchmarkable KS_EOM(
    $_state_ks, $_p_kep, $_ϕ_pt, $_dynamics_kep, $_config_pt
)
SUITE["eom"]["regularized"]["KS (Full)"] = @benchmarkable KS_EOM(
    $_state_ks, $_p_full, $_ϕ_pt, $_dynamics_full, $_config_pt
)

# --- Stiefel-Scheifele ---
SUITE["eom"]["regularized"]["StiSche (Keplerian)"] = @benchmarkable StiSche_EOM(
    $_state_stische, $_p_kep, $_ϕ_pt, $_dynamics_kep, $_config_pt
)
SUITE["eom"]["regularized"]["StiSche (Full)"] = @benchmarkable StiSche_EOM(
    $_state_stische, $_p_full, $_ϕ_pt, $_dynamics_full, $_config_pt
)

# --- GEqOE ---
SUITE["eom"]["regularized"]["GEqOE (Keplerian)"] = @benchmarkable GEqOE_EOM(
    $_state_geqoe, $_p_kep, $_t, $_dynamics_kep, $_config_geqoe
)
SUITE["eom"]["regularized"]["GEqOE (Full)"] = @benchmarkable GEqOE_EOM(
    $_state_geqoe, $_p_full, $_t, $_dynamics_full, $_config_geqoe
)

# =====================
# Propagation benchmarks (2-day, full force model)
# =====================
const _tspan = (0.0, 2.0 * 86400.0)

const _W_full =
    potential(Cartesian(_state_cart), _p_full, 0.0, _grav_model) -
    potential(Cartesian(_state_cart), _p_full, 0.0, KeplerianGravityAstroModel(; μ=_μ_full))

const _config_prop = RegularizedCoordinateConfig(
    _state_cart, _μ_full; W=_W_full, t₀=0.0, flag_time=PhysicalTime()
)
const _config_geqoe_prop = RegularizedCoordinateConfig(; W=_W_full)

const _u0_cowell = copy(_state_cart)
const _u0_gaussve = Array(Keplerian(Cartesian(_state_cart), _μ_full))
const _u0_milankovich = Array(Milankovich(Cartesian(_state_cart), _μ_full))
const _u0_usm7 = Array(USM7(Cartesian(_state_cart), _μ_full))
const _u0_usm6 = Array(USM6(Cartesian(_state_cart), _μ_full))
const _u0_usmem = Array(USMEM(Cartesian(_state_cart), _μ_full))

const _ϕ_prop = AstroCoords.compute_initial_phi(_state_cart, _μ_full, _config_prop)
const _u0_edromo = Array(EDromo(Cartesian(_state_cart), _μ_full, _ϕ_prop, _config_prop))
const _u0_ks = Array(KustaanheimoStiefel(Cartesian(_state_cart), _μ_full, _config_prop))
const _u0_stische = Array(
    StiefelScheifele(Cartesian(_state_cart), _μ_full, _ϕ_prop, _config_prop)
)
const _u0_geqoe = Array(GEqOE(Cartesian(_state_cart), _μ_full, _config_geqoe_prop))

# --- Standard propagators ---
SUITE["propagation"]["Cowell"] = @benchmarkable propagate(
    $_u0_cowell, $_p_full, $_dynamics_full, $_tspan;
    prop_type=CowellPropagator(),
)
SUITE["propagation"]["GaussVE"] = @benchmarkable propagate(
    $_u0_gaussve, $_p_full, $_dynamics_full, $_tspan;
    prop_type=GaussVEPropagator(),
)
SUITE["propagation"]["Milankovich"] = @benchmarkable propagate(
    $_u0_milankovich, $_p_full, $_dynamics_full, $_tspan;
    prop_type=MilankovichPropagator(),
)
SUITE["propagation"]["USM7"] = @benchmarkable propagate(
    $_u0_usm7, $_p_full, $_dynamics_full, $_tspan;
    prop_type=USM7Propagator(),
)
SUITE["propagation"]["USM6"] = @benchmarkable propagate(
    $_u0_usm6, $_p_full, $_dynamics_full, $_tspan;
    prop_type=USM6Propagator(),
)
SUITE["propagation"]["USMEM"] = @benchmarkable propagate(
    $_u0_usmem, $_p_full, $_dynamics_full, $_tspan;
    prop_type=USMEMPropagator(),
)

# --- Regularized propagators ---
SUITE["propagation"]["EDromo"] = @benchmarkable propagate(
    $_u0_edromo, $_p_full, $_dynamics_full, $_tspan;
    prop_type=EDromoPropagator(), config=$_config_prop,
)
SUITE["propagation"]["KS"] = @benchmarkable propagate(
    $_u0_ks, $_p_full, $_dynamics_full, $_tspan;
    prop_type=KSPropagator(), config=$_config_prop,
)
SUITE["propagation"]["StiSche"] = @benchmarkable propagate(
    $_u0_stische, $_p_full, $_dynamics_full, $_tspan;
    prop_type=StiSchePropagator(), config=$_config_prop,
)
SUITE["propagation"]["GEqOE"] = @benchmarkable propagate(
    $_u0_geqoe, $_p_full, $_dynamics_full, $_tspan;
    prop_type=GEqOEPropagator(), config=$_config_geqoe_prop,
)

# ---------------------
# Tune and cache
# ---------------------
paramspath = joinpath(dirname(@__FILE__), "params.json")

if isfile(paramspath)
    loadparams!(SUITE, BenchmarkTools.load(paramspath)[1], :evals)
else
    tune!(SUITE)
    BenchmarkTools.save(paramspath, BenchmarkTools.params(SUITE))
end
