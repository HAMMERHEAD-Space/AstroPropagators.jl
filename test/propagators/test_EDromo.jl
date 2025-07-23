@testset "EDromo Propagator Keplerian Physical Time" begin
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
    edromo_config = RegularizedCoordinateConfig(
        u0_cart, μ; W=W, t₀=0.0, flag_time=PhysicalTime()
    )

    p_full = ComponentVector(; p..., μ=μ)

    ϕ₀ = compute_initial_phi(u0_cart, μ, edromo_config)

    u0_EDromo = Array(EDromo(Cartesian(u0_cart), μ, ϕ₀, edromo_config))

    model_list = CentralBodyDynamicsModel(grav_model)
    # The independent variable is ϕ, so we integrate over one orbit (2π)
    tspan = (ϕ₀, ϕ₀ + 6π)

    EOM!(du, u, p, t) = EDromo_EOM!(du, u, p, t, model_list, edromo_config)

    prob = ODEProblem(EOM!, u0_EDromo, tspan, p_full)
    sol = solve(
        prob,
        VCABM();
        abstol=1e-15,
        reltol=1e-15,
        callback=end_EDromo_integration(86400.0, edromo_config),
    )

    NRG = orbitalNRG.(EDromo.(sol.u), μ, sol.t, [edromo_config])
    @test NRG[1] ≈ NRG[end]

    final_state = Cartesian(EDromo(sol.u[end]), μ, sol.t[end], edromo_config)

    expected_end = [
        29447.829229065504
        21027.31807433234
        -1675.1455650359862
        0.1548633780130051
        2.3814564036944668
        0.1401977642923555
    ]
    @test final_state ≈ expected_end rtol=1e-6
end

@testset "EDromo Propagator Keplerian Constant Time" begin
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
    W =
        potential(Cartesian(u0_cart), p, 0.0, grav_model) -
        potential(Cartesian(u0_cart), p, 0.0, KeplerianGravityAstroModel(μ=μ))
    edromo_config = RegularizedCoordinateConfig(
        u0_cart, μ; W=W, t₀=0.0, flag_time=ConstantTime()
    )

    p_full = ComponentVector(; p..., μ=μ)

    ϕ₀ = compute_initial_phi(u0_cart, μ, edromo_config)

    u0_EDromo = Array(EDromo(Cartesian(u0_cart), μ, ϕ₀, edromo_config))

    model_list = CentralBodyDynamicsModel(grav_model)
    # The independent variable is ϕ, so we integrate over one orbit (2π)
    tspan = (ϕ₀, ϕ₀ + 6π)

    EOM!(du, u, p, t) = EDromo_EOM!(du, u, p, t, model_list, edromo_config)

    prob = ODEProblem(EOM!, u0_EDromo, tspan, p_full)
    sol = solve(
        prob,
        VCABM();
        abstol=1e-15,
        reltol=1e-15,
        callback=end_EDromo_integration(86400.0, edromo_config),
    )

    NRG = zeros(length(sol.u))
    for i in 1:length(sol.u)
        NRG[i] = orbitalNRG(EDromo(sol.u[i]), grav_model.μ, sol.t[i], edromo_config)
    end

    states = zeros(6, length(sol.u))
    for i in 1:length(sol.u)
        states[:, i] = Array(Cartesian(EDromo(sol.u[i]), μ, sol.t[i], edromo_config))
    end

    @test NRG[1] ≈ NRG[end]

    final_state = Cartesian(EDromo(sol.u[end]), μ, sol.t[end], edromo_config)

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

@testset "EDromo Propagator Keplerian Linear Time" begin
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
    W =
        potential(Cartesian(u0_cart), p, 0.0, grav_model) -
        potential(Cartesian(u0_cart), p, 0.0, KeplerianGravityAstroModel(μ=μ))
    edromo_config = RegularizedCoordinateConfig(
        u0_cart, μ; W=W, t₀=0.0, flag_time=LinearTime()
    )

    p_full = ComponentVector(; p..., μ=μ)

    ϕ₀ = compute_initial_phi(u0_cart, μ, edromo_config)

    u0_EDromo = Array(EDromo(Cartesian(u0_cart), μ, ϕ₀, edromo_config))

    model_list = CentralBodyDynamicsModel(grav_model)
    # The independent variable is ϕ, so we integrate over one orbit (2π)
    tspan = (ϕ₀, ϕ₀ + 6π)

    EOM!(du, u, p, t) = EDromo_EOM!(du, u, p, t, model_list, edromo_config)

    prob = ODEProblem(EOM!, u0_EDromo, tspan, p_full)
    sol = solve(
        prob,
        VCABM();
        abstol=1e-15,
        reltol=1e-15,
        callback=end_EDromo_integration(86400.0, edromo_config),
    )

    NRG = zeros(length(sol.u))
    for i in 1:length(sol.u)
        NRG[i] = orbitalNRG(EDromo(sol.u[i]), grav_model.μ, sol.t[i], edromo_config)
    end

    states = zeros(6, length(sol.u))
    for i in 1:length(sol.u)
        states[:, i] = Array(Cartesian(EDromo(sol.u[i]), μ, sol.t[i], edromo_config))
    end

    @test NRG[1] ≈ NRG[end]

    final_state = Cartesian(EDromo(sol.u[end]), μ, sol.t[end], edromo_config)

    expected_end = [
        29447.829229065504
        21027.31807433234
        -1675.1455650359862
        0.1548633780130051
        2.3814564036944668
        0.1401977642923555
    ]
    #TODO: IS THIS A BUG?
    @test final_state ≈ expected_end rtol=2e-1
end

@testset "EDromo Propagator High-Fidelity Physical Time" begin
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
    edromo_config = RegularizedCoordinateConfig(
        u0, μ; W=W, t₀=0.0, flag_time=PhysicalTime()
    )

    ϕ₀ = compute_initial_phi(u0, μ, edromo_config)

    tspan = (ϕ₀, ϕ₀ + 6π)

    p_full = ComponentVector(; p..., μ=μ)

    u0_EDromo = Array(EDromo(Cartesian(u0), μ, ϕ₀, edromo_config))

    EOM!(du, u, p, t) = EDromo_EOM!(du, u, p, t, model_list, edromo_config)

    prob = ODEProblem(EOM!, u0_EDromo, tspan, p_full)
    sol = solve(
        prob,
        VCABM();
        abstol=1e-13,
        reltol=1e-13,
        callback=end_EDromo_integration(86400.0, edromo_config),
    )

    sol.u

    states = zeros(6, length(sol.u))
    for i in 1:length(sol.u)
        states[:, i] = Array(Cartesian(EDromo(sol.u[i]), μ, sol.t[i], edromo_config))
    end

    final_state = Cartesian(EDromo(sol.u[end]), μ, sol.t[end], edromo_config)

    # Regression Test
    expected_end = [
        29245.74497253034
        22127.193058906043
        -1549.7695621180876
        0.03440754290088714
        2.3126076757200003
        0.1501127574349222
    ]
    @test final_state ≈ expected_end rtol=1e-6
end

@testset "EDromo Propagator High-Fidelity Constant Time" begin
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

    edromo_config = RegularizedCoordinateConfig(
        u0,
        μ;
        W=potential(Cartesian(u0), p, 0.0, grav_model) -
          potential(Cartesian(u0), p, 0.0, KeplerianGravityAstroModel(μ=μ)),
        flag_time=ConstantTime(),
    )

    ϕ₀ = compute_initial_phi(u0, μ, edromo_config)

    tspan = (ϕ₀, ϕ₀ + 6π)

    p_full = ComponentVector(; p..., μ=μ)

    u0_EDromo = Array(EDromo(Cartesian(u0), μ, ϕ₀, edromo_config))

    EOM!(du, u, p, t) = EDromo_EOM!(du, u, p, t, model_list, edromo_config)

    prob = ODEProblem(EOM!, u0_EDromo, tspan, p_full)
    sol = solve(
        prob,
        VCABM();
        abstol=1e-13,
        reltol=1e-13,
        callback=end_EDromo_integration(86400.0, edromo_config),
    )

    sol.u

    states = zeros(6, length(sol.u))
    for i in 1:length(sol.u)
        states[:, i] = Array(Cartesian(EDromo(sol.u[i]), μ, sol.t[i], edromo_config))
    end

    final_state = Cartesian(EDromo(sol.u[end]), μ, sol.t[end], edromo_config)

    # Regression Test
    expected_end = [
        29245.74497253034
        22127.193058906043
        -1549.7695621180876
        0.03440754290088714
        2.3126076757200003
        0.1501127574349222
    ]
    @test final_state ≈ expected_end rtol=1e-6
end

@testset "EDromo Propagator High-Fidelity Linear Time" begin
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

    edromo_config = RegularizedCoordinateConfig(
        u0,
        μ;
        W=potential(Cartesian(u0), p, 0.0, grav_model) -
          potential(Cartesian(u0), p, 0.0, KeplerianGravityAstroModel(μ=μ)),
        flag_time=LinearTime(),
    )

    ϕ₀ = compute_initial_phi(u0, μ, edromo_config)

    tspan = (ϕ₀, ϕ₀ + 6π)

    p_full = ComponentVector(; p..., μ=μ)

    u0_EDromo = Array(EDromo(Cartesian(u0), μ, ϕ₀, edromo_config))

    EOM!(du, u, p, t) = EDromo_EOM!(du, u, p, t, model_list, edromo_config)

    prob = ODEProblem(EOM!, u0_EDromo, tspan, p_full)
    sol = solve(
        prob,
        VCABM();
        abstol=1e-13,
        reltol=1e-13,
        callback=end_EDromo_integration(86400.0, edromo_config),
    )

    states = zeros(6, length(sol.u))
    for i in 1:length(sol.u)
        states[:, i] = Array(Cartesian(EDromo(sol.u[i]), μ, sol.t[i], edromo_config))
    end

    final_state = Cartesian(EDromo(sol.u[end]), μ, sol.t[end], edromo_config)

    # Regression Test
    expected_end = [
        29245.74497253034
        22127.193058906043
        -1549.7695621180876
        0.03440754290088714
        2.3126076757200003
        0.1501127574349222
    ]
    @test final_state ≈ expected_end rtol=1e-6
end

@testset "EDromo Propagator High-Fidelity Physical Time 2" begin
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

    edromo_config = RegularizedCoordinateConfig(
        u0,
        μ;
        W=potential(Cartesian(u0), p, 0.0, grav_model) -
          potential(Cartesian(u0), p, 0.0, KeplerianGravityAstroModel(μ=μ)),
        flag_time=PhysicalTime(),
    )

    ϕ₀ = compute_initial_phi(u0, μ, edromo_config)

    tspan = (ϕ₀, ϕ₀ + 50π)

    p_full = ComponentVector(; p..., μ=μ)

    u0_EDromo = Array(EDromo(Cartesian(u0), μ, ϕ₀, edromo_config))

    EOM!(du, u, p, t) = EDromo_EOM!(du, u, p, t, model_list, edromo_config)

    prob = ODEProblem(EOM!, u0_EDromo, tspan, p_full)
    sol = solve(
        prob,
        VCABM();
        abstol=1e-13,
        reltol=1e-13,
        callback=end_EDromo_integration(3.0 * 86400.0, edromo_config),
    )

    @assert sol.retcode == ReturnCode.Terminated

    states = zeros(6, length(sol.u))
    for i in 1:length(sol.u)
        states[:, i] = Array(Cartesian(EDromo(sol.u[i]), μ, sol.t[i], edromo_config))
    end

    final_state = Cartesian(EDromo(sol.u[end]), μ, sol.t[end], edromo_config)

    # Regression Test
    expected_end = [
        -8372.895299659556
        5579.539458768319
        1276.7499246213724
        -0.41955615768671234
        -7.385062119946734
        -0.4797259172503372
    ]
    @test final_state ≈ expected_end rtol=1e-6
end

@testset "EDromo Propagator High-Fidelity Constant Time 2" begin
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

    edromo_config = RegularizedCoordinateConfig(
        u0,
        μ;
        W=potential(Cartesian(u0), p, 0.0, grav_model) -
          potential(Cartesian(u0), p, 0.0, KeplerianGravityAstroModel(μ=μ)),
        flag_time=ConstantTime(),
    )

    ϕ₀ = compute_initial_phi(u0, μ, edromo_config)

    tspan = (ϕ₀, ϕ₀ + 50π)

    p_full = ComponentVector(; p..., μ=μ)

    u0_EDromo = Array(EDromo(Cartesian(u0), μ, ϕ₀, edromo_config))

    EOM!(du, u, p, t) = EDromo_EOM!(du, u, p, t, model_list, edromo_config)

    prob = ODEProblem(EOM!, u0_EDromo, tspan, p_full)
    sol = solve(
        prob,
        VCABM();
        abstol=1e-13,
        reltol=1e-13,
        callback=end_EDromo_integration(3.0 * 86400.0, edromo_config),
    )

    @assert sol.retcode == ReturnCode.Terminated

    states = zeros(6, length(sol.u))
    for i in 1:length(sol.u)
        states[:, i] = Array(Cartesian(EDromo(sol.u[i]), μ, sol.t[i], edromo_config))
    end

    final_state = Cartesian(EDromo(sol.u[end]), μ, sol.t[end], edromo_config)

    # Regression Test
    expected_end = [
        -8372.895299659556
        5579.539458768319
        1276.7499246213724
        -0.41955615768671234
        -7.385062119946734
        -0.4797259172503372
    ]
    @test final_state ≈ expected_end rtol=1e-6
end

@testset "EDromo Propagator High-Fidelity Linear Time 2" begin
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

    edromo_config = RegularizedCoordinateConfig(
        u0,
        μ;
        W=potential(Cartesian(u0), p, 0.0, grav_model) -
          potential(Cartesian(u0), p, 0.0, KeplerianGravityAstroModel(μ=μ)),
        flag_time=LinearTime(),
    )

    ϕ₀ = compute_initial_phi(u0, μ, edromo_config)

    tspan = (ϕ₀, ϕ₀ + 50π)

    p_full = ComponentVector(; p..., μ=μ)

    u0_EDromo = Array(EDromo(Cartesian(u0), μ, ϕ₀, edromo_config))

    EOM!(du, u, p, t) = EDromo_EOM!(du, u, p, t, model_list, edromo_config)

    prob = ODEProblem(EOM!, u0_EDromo, tspan, p_full)
    sol = solve(
        prob,
        VCABM();
        abstol=1e-13,
        reltol=1e-13,
        callback=end_EDromo_integration(3.0 * 86400.0, edromo_config),
    )

    @assert sol.retcode == ReturnCode.Terminated

    states = zeros(6, length(sol.u))
    for i in 1:length(sol.u)
        states[:, i] = Array(Cartesian(EDromo(sol.u[i]), μ, sol.t[i], edromo_config))
    end

    final_state = Cartesian(EDromo(sol.u[end]), μ, sol.t[end], edromo_config)

    # Regression Test
    expected_end = [
        -8372.895299659556
        5579.539458768319
        1276.7499246213724
        -0.41955615768671234
        -7.385062119946734
        -0.4797259172503372
    ]
    @test final_state ≈ expected_end rtol=1e-6
end
