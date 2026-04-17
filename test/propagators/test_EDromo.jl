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
        Vern9();
        abstol=1e-15,
        reltol=1e-15,
        callback=build_termination_callback(86400.0, EDromo, edromo_config),
    )

    NRG = orbitalNRG.(EDromo.(sol.u), μ, sol.t, [edromo_config])
    @test NRG[1] ≈ NRG[end]
    h = norm.(angularMomentumVector.(EDromo.(sol.u), μ, sol.t, [edromo_config]))
    @test h[1] ≈ h[end]

    final_state = Cartesian(EDromo(sol.u[end]), μ, sol.t[end], edromo_config)

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
        Vern9();
        abstol=1e-15,
        reltol=1e-15,
        callback=build_termination_callback(86400.0, EDromo, edromo_config),
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
    h = norm.(angularMomentumVector.(EDromo.(sol.u), μ, sol.t, [edromo_config]))
    @test h[1] ≈ h[end]

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
        Vern9();
        abstol=1e-15,
        reltol=1e-15,
        callback=build_termination_callback(86400.0, EDromo, edromo_config),
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
    h = norm.(angularMomentumVector.(EDromo.(sol.u), μ, sol.t, [edromo_config]))
    @test h[1] ≈ h[end]

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
