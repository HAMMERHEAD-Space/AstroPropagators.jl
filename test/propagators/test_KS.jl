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
        Vern9();
        abstol=1e-15,
        reltol=1e-15,
        callback=end_KS_integration(86400.0, ks_config),
    )

    get_KS_time(sol.u[end], ks_config)

    NRG = orbitalNRG.(KustaanheimoStiefel.(sol.u), μ, [ks_config])
    @test NRG[1] ≈ NRG[end] rtol=1e-6
    h = norm.(angularMomentumVector.(KustaanheimoStiefel.(sol.u), μ, [ks_config]))
    @test h[1] ≈ h[end]

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
        Vern9();
        abstol=1e-15,
        reltol=1e-15,
        callback=end_KS_integration(86400.0, ks_config),
    )

    NRG = orbitalNRG.(KustaanheimoStiefel.(sol.u), μ, [ks_config])
    @test NRG[1] ≈ NRG[end] rtol=1e-6
    h = norm.(angularMomentumVector.(KustaanheimoStiefel.(sol.u), μ, [ks_config]))
    @test h[1] ≈ h[end]

    final_state = Cartesian(KustaanheimoStiefel(sol.u[end]), μ, ks_config)

    expected_end = [
        29447.829229065504
        21027.31807433234
        -1675.1455650359862
        0.1548633780130051
        2.3814564036944668
        0.1401977642923555
    ]
    @test final_state ≈ expected_end rtol=1e-10
end

@testset "Kustaanheimo-Stiefel Propagator High-Fidelity Regression Physical Time" begin
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
        Vern9();
        abstol=1e-13,
        reltol=1e-13,
        callback=end_KS_integration(3.0 * 86400.0, ks_config),
    )

    @assert sol.retcode == ReturnCode.Terminated

    final_state = Cartesian(KustaanheimoStiefel(sol.u[end]), μ, ks_config)

    # Regression Test
    expected_end = [
        -6450.191420482322
        -2390.1543192175777
        500.30668770666614
        4.812541976308283
        -7.890469898354672
        -1.0701929804778323
    ]
    @test final_state ≈ expected_end rtol=1e-3
end

@testset "Kustaanheimo-Stiefel Propagator High-Fidelity Regression Linear Time" begin
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
        Vern9();
        abstol=1e-13,
        reltol=1e-13,
        callback=end_KS_integration(3.0 * 86400.0, ks_config),
    )

    @assert sol.retcode == ReturnCode.Terminated

    final_state = Cartesian(KustaanheimoStiefel(sol.u[end]), μ, ks_config)

    # Regression Test
    expected_end = [
        -6450.189926464047
        -2390.1567676103496
        500.3063554871579
        4.812544408893667
        -7.890468995144382
        -1.070193169485886
    ]
    @test final_state ≈ expected_end rtol=1e-3
end
