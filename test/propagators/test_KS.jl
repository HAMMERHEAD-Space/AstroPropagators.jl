@testset "Kustaanheimo-Stiefel Propagator Keplerian Physical Time" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
    p = ComponentVector(; JD=JD)

    grav_model = KeplerianGravityAstroModel()
    μ = grav_model.μ

    u0_cart = [
        -1076.225324679696
        -6765.896364327722
        -332.3087833503755
        9.356857417032581
        -3.3123476319597557
        -1.1880157328553503
    ] #km, km/s

    # For a pure Keplerian problem, the perturbing potential W₀ is 0.
    W = (
        potential(Cartesian(u0_cart), p, 0.0, grav_model) -
        potential(Cartesian(u0_cart), p, 0.0, KeplerianGravityAstroModel(μ=μ))
    )
    ks_config = RegularizedCoordinateConfig(
        u0_cart, μ; W=W, t₀=0.0, flag_time=PhysicalTime()
    )

    p_full = ComponentVector(; p..., μ=μ)

    u0_KS = Array(KustaanheimoStiefel(Cartesian(u0_cart), μ, ks_config))

    model_list = CentralBodyDynamicsModel(grav_model)
    # The independent variable is ϕ, so we integrate over one orbit (2π)
    tspan = (0.0, 9π)

    EOM!(du, u, p, t) = KS_EOM!(du, u, p, t, model_list, ks_config)

    prob = ODEProblem(EOM!, u0_KS, tspan, p_full)
    sol = solve(
        prob,
        VCABM();
        abstol=1e-15,
        reltol=1e-15,
        callback=end_KS_integration(86400.0, ks_config),
    )

    get_KS_time(sol.u[end], ks_config)

    NRG = orbitalNRG.(KustaanheimoStiefel.(sol.u), μ, [ks_config])
    @test NRG[1] ≈ NRG[end] rtol=1e-6

    final_state = Cartesian(KustaanheimoStiefel(sol.u[end]), μ, ks_config)

    expected_end = [
        29447.829229065504
        21027.31807433234
        -1675.1455650359862
        0.1548633780130051
        2.3814564036944668
        0.1401977642923555
    ]
    @test final_state ≈ expected_end rtol=1e-3
end

@testset "Kustaanheimo-Stiefel Propagator Keplerian Linear Time" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
    p = ComponentVector(; JD=JD)

    grav_model = KeplerianGravityAstroModel()
    μ = grav_model.μ

    u0_cart = [
        -1076.225324679696
        -6765.896364327722
        -332.3087833503755
        9.356857417032581
        -3.3123476319597557
        -1.1880157328553503
    ] #km, km/s

    # For a pure Keplerian problem, the perturbing potential W₀ is 0.
    W = (
        potential(Cartesian(u0_cart), p, 0.0, grav_model) -
        potential(Cartesian(u0_cart), p, 0.0, KeplerianGravityAstroModel(μ=μ))
    )
    ks_config = RegularizedCoordinateConfig(u0_cart, μ; W=W, t₀=0.0, flag_time=LinearTime())

    p_full = ComponentVector(; p..., μ=μ)

    u0_KS = Array(KustaanheimoStiefel(Cartesian(u0_cart), μ, ks_config))

    model_list = CentralBodyDynamicsModel(grav_model)
    # The independent variable is ϕ, so we integrate over one orbit (2π)
    tspan = (0.0, 9π)

    EOM!(du, u, p, t) = KS_EOM!(du, u, p, t, model_list, ks_config)

    prob = ODEProblem(EOM!, u0_KS, tspan, p_full)
    sol = solve(
        prob,
        VCABM();
        abstol=1e-15,
        reltol=1e-15,
        callback=end_KS_integration(86400.0, ks_config),
    )

    NRG = orbitalNRG.(KustaanheimoStiefel.(sol.u), μ, [ks_config])
    #TODO: IS THIS A BUG?
    @test NRG[1] ≈ NRG[end] rtol=1e-3

    final_state = Cartesian(KustaanheimoStiefel(sol.u[end]), μ, ks_config)

    expected_end = [
        29447.829229065504
        21027.31807433234
        -1675.1455650359862
        0.1548633780130051
        2.3814564036944668
        0.1401977642923555
    ]
    @test final_state ≈ expected_end rtol=1e-3
end

@testset "Kustaanheimo-Stiefel Propagator High-Fidelity Physical Time" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
    p = ComponentVector(; JD=JD)

    SpaceIndices.init()
    eop_data = fetch_iers_eop()
    grav_coeffs = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96))

    grav_model = GravityHarmonicsAstroModel(;
        gravity_model=grav_coeffs, eop_data=eop_data, order=36, degree=36
    )
    sun_third_body = ThirdBodyModel(; body=SunBody(), eop_data=eop_data)
    moon_third_body = ThirdBodyModel(; body=MoonBody(), eop_data=eop_data)

    satellite_srp_model = CannonballFixedSRP(0.2)
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

    u0 = [
        -1076.225324679696
        -6765.896364327722
        -332.3087833503755
        9.356857417032581
        -3.3123476319597557
        -1.1880157328553503
    ] #km, km/s

    model_list = CentralBodyDynamicsModel(
        grav_model, (sun_third_body, moon_third_body, srp_model, drag_model)
    )

    μ = GravityModels.gravity_constant(grav_coeffs) / 1E9

    W =
        potential(Cartesian(u0), p, 0.0, grav_model) -
        potential(Cartesian(u0), p, 0.0, KeplerianGravityAstroModel(μ=μ))
    ks_config = RegularizedCoordinateConfig(u0, μ; W=W, t₀=0.0, flag_time=PhysicalTime())

    tspan = (0.0, 9π)

    p_full = ComponentVector(; p..., μ=μ)

    u0_ks = Array(KustaanheimoStiefel(Cartesian(u0), μ, ks_config))

    EOM!(du, u, p, t) = KS_EOM!(du, u, p, t, model_list, ks_config)

    prob = ODEProblem(EOM!, u0_ks, tspan, p_full)
    sol = solve(
        prob,
        VCABM();
        abstol=1e-13,
        reltol=1e-13,
        callback=end_KS_integration(86400.0, ks_config),
    )

    final_state = Cartesian(KustaanheimoStiefel(sol.u[end]), μ, ks_config)

    # Regression Test
    expected_end = [
        29209.156599133385
        22221.211933611725
        -1539.730469049243
        0.020136024262479676
        2.3045202516241976
        0.15104860229599937
    ]
    @test final_state ≈ expected_end rtol=1e-3
end

@testset "Kustaanheimo-Stiefel Propagator High-Fidelity Linear Time" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
    p = ComponentVector(; JD=JD)

    SpaceIndices.init()
    eop_data = fetch_iers_eop()
    grav_coeffs = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96))

    grav_model = GravityHarmonicsAstroModel(;
        gravity_model=grav_coeffs, eop_data=eop_data, order=36, degree=36
    )
    sun_third_body = ThirdBodyModel(; body=SunBody(), eop_data=eop_data)
    moon_third_body = ThirdBodyModel(; body=MoonBody(), eop_data=eop_data)

    satellite_srp_model = CannonballFixedSRP(0.2)
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

    u0 = [
        -1076.225324679696
        -6765.896364327722
        -332.3087833503755
        9.356857417032581
        -3.3123476319597557
        -1.1880157328553503
    ] #km, km/s

    model_list = CentralBodyDynamicsModel(
        grav_model, (sun_third_body, moon_third_body, srp_model, drag_model)
    )

    μ = GravityModels.gravity_constant(grav_coeffs) / 1E9

    W =
        potential(Cartesian(u0), p, 0.0, grav_model) -
        potential(Cartesian(u0), p, 0.0, KeplerianGravityAstroModel(μ=μ))
    ks_config = RegularizedCoordinateConfig(u0, μ; W=W, t₀=0.0, flag_time=LinearTime())

    tspan = (0.0, 9π)

    p_full = ComponentVector(; p..., μ=μ)

    u0_ks = Array(KustaanheimoStiefel(Cartesian(u0), μ, ks_config))

    EOM!(du, u, p, t) = KS_EOM!(du, u, p, t, model_list, ks_config)

    prob = ODEProblem(EOM!, u0_ks, tspan, p_full)
    sol = solve(
        prob,
        VCABM();
        abstol=1e-13,
        reltol=1e-13,
        callback=end_KS_integration(86400.0, ks_config),
    )

    final_state = Cartesian(KustaanheimoStiefel(sol.u[end]), μ, ks_config)

    # Regression Test
    expected_end = [
        29209.15685582137
        22221.213061986924
        -1539.7303991049012
        0.020135930608283184
        2.304520171352468
        0.15104860706652987
    ]
    @test final_state ≈ expected_end rtol=1e-3
end

@testset "Kustaanheimo-Stiefel Propagator High-Fidelity Physical Time 2" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
    p = ComponentVector(; JD=JD)

    SpaceIndices.init()
    eop_data = fetch_iers_eop()
    grav_coeffs = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96))

    grav_model = GravityHarmonicsAstroModel(;
        gravity_model=grav_coeffs, eop_data=eop_data, order=36, degree=36
    )
    sun_third_body = ThirdBodyModel(; body=SunBody(), eop_data=eop_data)
    moon_third_body = ThirdBodyModel(; body=MoonBody(), eop_data=eop_data)

    satellite_srp_model = CannonballFixedSRP(0.5)
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

    u0 = [
        -1076.225324679696
        -6765.896364327722
        -332.3087833503755
        8.956857417032581
        -3.3123476319597557
        -1.1880157328553503
    ] #km, km/s

    model_list = CentralBodyDynamicsModel(
        grav_model, (sun_third_body, moon_third_body, srp_model, drag_model)
    )

    μ = GravityModels.gravity_constant(grav_coeffs) / 1E9

    ks_config = RegularizedCoordinateConfig(
        u0,
        μ;
        W=potential(Cartesian(u0), p, 0.0, grav_model) -
          potential(Cartesian(u0), p, 0.0, KeplerianGravityAstroModel(μ=μ)),
        flag_time=PhysicalTime(),
    )

    tspan = (0.0, 60π)

    p_full = ComponentVector(; p..., μ=μ)

    u0_ks = Array(KustaanheimoStiefel(Cartesian(u0), μ, ks_config))

    EOM!(du, u, p, t) = KS_EOM!(du, u, p, t, model_list, ks_config)

    prob = ODEProblem(EOM!, u0_ks, tspan, p_full)
    sol = solve(
        prob,
        VCABM();
        abstol=1e-13,
        reltol=1e-13,
        callback=end_KS_integration(3.0 * 86400.0, ks_config),
    )

    @assert sol.retcode == ReturnCode.Terminated

    final_state = Cartesian(KustaanheimoStiefel(sol.u[end]), μ, ks_config)

    # Regression Test
    expected_end = [
        -6464.325397344593
        -2366.920974325633
        503.4543800242743
        4.7894749700537345
        -7.898985518614695
        -1.0683968399621873
    ]
    @test final_state ≈ expected_end rtol=1e-3
end

@testset "Kustaanheimo-Stiefel Propagator High-Fidelity Linear Time 2" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
    p = ComponentVector(; JD=JD)

    SpaceIndices.init()
    eop_data = fetch_iers_eop()
    grav_coeffs = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96))

    grav_model = GravityHarmonicsAstroModel(;
        gravity_model=grav_coeffs, eop_data=eop_data, order=36, degree=36
    )
    sun_third_body = ThirdBodyModel(; body=SunBody(), eop_data=eop_data)
    moon_third_body = ThirdBodyModel(; body=MoonBody(), eop_data=eop_data)

    satellite_srp_model = CannonballFixedSRP(0.5)
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

    u0 = [
        -1076.225324679696
        -6765.896364327722
        -332.3087833503755
        8.956857417032581
        -3.3123476319597557
        -1.1880157328553503
    ] #km, km/s

    model_list = CentralBodyDynamicsModel(
        grav_model, (sun_third_body, moon_third_body, srp_model, drag_model)
    )

    μ = GravityModels.gravity_constant(grav_coeffs) / 1E9

    ks_config = RegularizedCoordinateConfig(
        u0,
        μ;
        W=potential(Cartesian(u0), p, 0.0, grav_model) -
          potential(Cartesian(u0), p, 0.0, KeplerianGravityAstroModel(μ=μ)),
        flag_time=LinearTime(),
    )

    tspan = (0.0, 60π)

    p_full = ComponentVector(; p..., μ=μ)

    u0_ks = Array(KustaanheimoStiefel(Cartesian(u0), μ, ks_config))

    EOM!(du, u, p, t) = KS_EOM!(du, u, p, t, model_list, ks_config)

    prob = ODEProblem(EOM!, u0_ks, tspan, p_full)
    sol = solve(
        prob,
        VCABM();
        abstol=1e-13,
        reltol=1e-13,
        callback=end_KS_integration(3.0 * 86400.0, ks_config),
    )

    @assert sol.retcode == ReturnCode.Terminated

    final_state = Cartesian(KustaanheimoStiefel(sol.u[end]), μ, ks_config)

    # Regression Test
    expected_end = [
        -6464.311859841494
        -2366.9431006298246
        503.4513556099948
        4.7894966732659325
        -7.89897738711462
        -1.0683985082163687
    ]
    @test final_state ≈ expected_end rtol=1e-3
end
