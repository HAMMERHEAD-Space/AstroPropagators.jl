const COMMON_JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)

mutable struct _MockIntegrator{UT,PT}
    u::UT
    t::Float64
    p::PT
end

const KEPLERIAN_U0 = [
    -1076.225324679696
    -6765.896364327722
    -332.3087833503755
    9.356857417032581
    -3.3123476319597557
    -1.1880157328553503
]

const HF2_U0 = [
    -1076.225324679696
    -6765.896364327722
    -332.3087833503755
    8.956857417032581
    -3.3123476319597557
    -1.1880157328553503
]

function negate_velocity(u)
    return [u[1], u[2], u[3], -u[4], -u[5], -u[6]]
end

# ── Force model setup ─────────────────────────────────────────────────

function setup_keplerian()
    grav_model = KeplerianGravityAstroModel()
    μ = grav_model.μ
    model_list = CentralBodyDynamicsModel(grav_model)
    p = ComponentVector(; JD=COMMON_JD, μ=μ)
    return (; model_list, p, μ, grav_model)
end

function setup_gravity_only()
    SpaceIndices.init()
    eop_data = fetch_iers_eop()
    grav_coeffs = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96))
    grav_model = GravityHarmonicsAstroModel(;
        gravity_model=grav_coeffs, eop_data=eop_data, order=36, degree=36
    )
    μ = GravityModels.gravity_constant(grav_model.gravity_model) / 1E9
    sun_third_body = ThirdBodyModel(; body=SunBody(), eop_data=eop_data)
    moon_third_body = ThirdBodyModel(; body=MoonBody(), eop_data=eop_data)
    model_list = CentralBodyDynamicsModel(grav_model, (sun_third_body, moon_third_body))
    p = ComponentVector(; JD=COMMON_JD, μ=μ)
    return (; model_list, p, μ, grav_model, eop_data)
end

function setup_full_hf(; srp_coeff=0.2)
    SpaceIndices.init()
    eop_data = fetch_iers_eop()
    grav_coeffs = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96))
    grav_model = GravityHarmonicsAstroModel(;
        gravity_model=grav_coeffs, eop_data=eop_data, order=36, degree=36
    )
    μ = GravityModels.gravity_constant(grav_model.gravity_model) / 1E9
    sun_third_body = ThirdBodyModel(; body=SunBody(), eop_data=eop_data)
    moon_third_body = ThirdBodyModel(; body=MoonBody(), eop_data=eop_data)
    satellite_srp_model = CannonballFixedSRP(srp_coeff)
    srp_model = SRPAstroModel(;
        satellite_srp_model=satellite_srp_model,
        sun_data=sun_third_body,
        eop_data=eop_data,
        shadow_model=Conical(),
    )
    satellite_drag_model = CannonballFixedDrag(0.2)
    drag_model = DragAstroModel(;
        satellite_drag_model=satellite_drag_model,
        atmosphere_model=JB2008(),
        eop_data=eop_data,
    )
    model_list = CentralBodyDynamicsModel(
        grav_model, (sun_third_body, moon_third_body, srp_model, drag_model)
    )
    p = ComponentVector(; JD=COMMON_JD, μ=μ)
    return (; model_list, p, μ, grav_model, eop_data)
end

# ── Standard propagator wrappers ──────────────────────────────────────

function run_cowell(u0, model_list, p, tspan; atol=1e-13, rtol=1e-13)
    EOM!(du, u, _p, t) = Cowell_EOM!(du, u, _p, t, model_list)
    prob = ODEProblem(EOM!, u0, tspan, p)
    sol = solve(prob, Vern9(); abstol=atol, reltol=rtol)
    return sol.u[end]
end

function run_gaussve(u0_cart, model_list, p, tspan; atol=1e-13, rtol=1e-13)
    μ = p.μ
    u0_koe = Array(Keplerian(Cartesian(u0_cart), μ))
    EOM!(du, u, _p, t) = GaussVE_EOM!(du, u, _p, t, model_list)
    prob = ODEProblem(EOM!, u0_koe, tspan, p)
    sol = solve(prob, Vern9(); abstol=atol, reltol=rtol)
    return Array(Cartesian(Keplerian(sol.u[end]), μ))
end

function run_milankovich(u0_cart, model_list, p, tspan; atol=1e-13, rtol=1e-13)
    μ = p.μ
    u0_mil = Array(Milankovich(Cartesian(u0_cart), μ))
    EOM!(du, u, _p, t) = Milankovich_EOM!(du, u, _p, t, model_list)
    prob = ODEProblem(EOM!, u0_mil, tspan, p)
    sol = solve(prob, Vern9(); abstol=atol, reltol=rtol)
    return Array(Cartesian(Milankovich(sol.u[end]), μ))
end

function run_usm7(u0_cart, model_list, p, tspan; atol=1e-13, rtol=1e-13)
    μ = p.μ
    u0_usm = Array(USM7(Cartesian(u0_cart), μ))
    EOM!(du, u, _p, t) = USM7_EOM!(du, u, _p, t, model_list)
    prob = ODEProblem(EOM!, u0_usm, tspan, p)
    sol = solve(prob, Vern9(); abstol=atol, reltol=rtol)
    return Array(Cartesian(USM7(sol.u[end]), μ))
end

function run_usm6(u0_cart, model_list, p, tspan; atol=1e-13, rtol=1e-13)
    μ = p.μ
    u0_usm = Array(USM6(Cartesian(u0_cart), μ))
    EOM!(du, u, _p, t) = USM6_EOM!(du, u, _p, t, model_list)
    prob = ODEProblem(EOM!, u0_usm, tspan, p)
    sol = solve(prob, Vern9(); abstol=atol, reltol=rtol)
    return Array(Cartesian(USM6(sol.u[end]), μ))
end

function run_usmem(u0_cart, model_list, p, tspan; atol=1e-13, rtol=1e-13)
    μ = p.μ
    u0_usm = Array(USMEM(Cartesian(u0_cart), μ))
    EOM!(du, u, _p, t) = USMEM_EOM!(du, u, _p, t, model_list)
    prob = ODEProblem(EOM!, u0_usm, tspan, p)
    sol = solve(prob, Vern9(); abstol=atol, reltol=rtol)
    return Array(Cartesian(USMEM(sol.u[end]), μ))
end

# ── Regularized propagator wrappers ───────────────────────────────────

function _regularized_W(u0_cart, grav_model, μ)
    p_base = ComponentVector(; JD=COMMON_JD)
    return (
        potential(Cartesian(u0_cart), p_base, 0.0, grav_model) -
        potential(Cartesian(u0_cart), p_base, 0.0, KeplerianGravityAstroModel(; μ=μ))
    )
end

function run_edromo(
    u0_cart,
    model_list,
    μ,
    grav_model,
    duration;
    flag_time=PhysicalTime(),
    atol=1e-15,
    rtol=1e-15,
    ϕ_range=6π,
)
    W = _regularized_W(u0_cart, grav_model, μ)
    config = RegularizedCoordinateConfig(u0_cart, μ; W=W, t₀=0.0, flag_time=flag_time)
    p_full = ComponentVector(; JD=COMMON_JD, μ=μ)
    ϕ₀ = compute_initial_phi(u0_cart, μ, config)
    u0 = Array(EDromo(Cartesian(u0_cart), μ, ϕ₀, config))
    tspan = (ϕ₀, ϕ₀ + ϕ_range)
    EOM!(du, u, _p, t) = EDromo_EOM!(du, u, _p, t, model_list, config)
    prob = ODEProblem(EOM!, u0, tspan, p_full)
    sol = solve(
        prob,
        Vern9();
        abstol=atol,
        reltol=rtol,
        callback=build_termination_callback(duration, EDromo, config),
    )
    return Array(Cartesian(EDromo(sol.u[end]), μ, sol.t[end], config))
end

function run_ks(
    u0_cart,
    model_list,
    μ,
    grav_model,
    duration;
    flag_time=PhysicalTime(),
    atol=1e-15,
    rtol=1e-15,
    ϕ_range=9π,
)
    W = _regularized_W(u0_cart, grav_model, μ)
    config = RegularizedCoordinateConfig(u0_cart, μ; W=W, t₀=0.0, flag_time=flag_time)
    p_full = ComponentVector(; JD=COMMON_JD, μ=μ)
    u0 = Array(KustaanheimoStiefel(Cartesian(u0_cart), μ, config))
    tspan = (0.0, ϕ_range)
    EOM!(du, u, _p, t) = KS_EOM!(du, u, _p, t, model_list, config)
    prob = ODEProblem(EOM!, u0, tspan, p_full)
    sol = solve(
        prob,
        Vern9();
        abstol=atol,
        reltol=rtol,
        callback=build_termination_callback(duration, KustaanheimoStiefel, config),
    )
    return Array(Cartesian(KustaanheimoStiefel(sol.u[end]), μ, config))
end

function run_stische(
    u0_cart,
    model_list,
    μ,
    grav_model,
    duration;
    flag_time=PhysicalTime(),
    atol=1e-15,
    rtol=1e-15,
    ϕ_range=6π,
)
    W = _regularized_W(u0_cart, grav_model, μ)
    config = RegularizedCoordinateConfig(u0_cart, μ; W=W, t₀=0.0, flag_time=flag_time)
    p_full = ComponentVector(; JD=COMMON_JD, μ=μ)
    ϕ₀ = compute_initial_phi(u0_cart, μ, config)
    u0 = Array(StiefelScheifele(Cartesian(u0_cart), μ, ϕ₀, config))
    tspan = (ϕ₀, ϕ₀ + ϕ_range)
    EOM!(du, u, _p, t) = StiSche_EOM!(du, u, _p, t, model_list, config)
    prob = ODEProblem(EOM!, u0, tspan, p_full)
    sol = solve(
        prob,
        Vern9();
        abstol=atol,
        reltol=rtol,
        callback=build_termination_callback(duration, StiefelScheifele, config),
    )
    return Array(Cartesian(StiefelScheifele(sol.u[end]), μ, sol.t[end], config))
end
