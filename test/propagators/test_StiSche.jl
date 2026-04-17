@testset "Stiefel-Scheifele Propagator Keplerian Physical Time" begin
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
    stische_config = RegularizedCoordinateConfig(
        u0_cart, μ; W=W, t₀=0.0, flag_time=PhysicalTime()
    )

    p_full = ComponentVector(; p..., μ=μ)

    ϕ₀ = compute_initial_phi(u0_cart, μ, stische_config)

    u0_stische = Array(StiefelScheifele(Cartesian(u0_cart), μ, ϕ₀, stische_config))

    model_list = CentralBodyDynamicsModel(grav_model)
    # The independent variable is ϕ, so we integrate over one orbit (2π)
    tspan = (ϕ₀, ϕ₀ + 6π)

    EOM!(du, u, p, t) = StiSche_EOM!(du, u, p, t, model_list, stische_config)

    prob = ODEProblem(EOM!, u0_stische, tspan, p_full)
    sol = solve(
        prob,
        Vern9();
        abstol=1e-15,
        reltol=1e-15,
        callback=build_termination_callback(86400.0, StiefelScheifele, stische_config),
    )

    NRG = orbitalNRG.(StiefelScheifele.(sol.u), μ, sol.t, [stische_config])
    @test NRG[1] ≈ NRG[end]
    h = norm.(angularMomentumVector.(StiefelScheifele.(sol.u), μ, sol.t, [stische_config]))
    @test h[1] ≈ h[end]

    final_state = Cartesian(StiefelScheifele(sol.u[end]), μ, sol.t[end], stische_config)

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

@testset "Stiefel-Scheifele Propagator Keplerian Linear Time" begin
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
    stische_config = RegularizedCoordinateConfig(
        u0_cart, μ; W=W, t₀=0.0, flag_time=LinearTime()
    )

    p_full = ComponentVector(; p..., μ=μ)

    ϕ₀ = compute_initial_phi(u0_cart, μ, stische_config)

    u0_stische = Array(StiefelScheifele(Cartesian(u0_cart), μ, ϕ₀, stische_config))

    model_list = CentralBodyDynamicsModel(grav_model)
    # The independent variable is ϕ, so we integrate over one orbit (2π)
    tspan = (ϕ₀, ϕ₀ + 6π)

    EOM!(du, u, p, t) = StiSche_EOM!(du, u, p, t, model_list, stische_config)

    prob = ODEProblem(EOM!, u0_stische, tspan, p_full)
    sol = solve(
        prob,
        Vern9();
        abstol=1e-15,
        reltol=1e-15,
        callback=build_termination_callback(86400.0, StiefelScheifele, stische_config),
    )

    NRG = orbitalNRG.(StiefelScheifele.(sol.u), μ, sol.t, [stische_config])
    @test NRG[1] ≈ NRG[end]
    h = norm.(angularMomentumVector.(StiefelScheifele.(sol.u), μ, sol.t, [stische_config]))
    @test h[1] ≈ h[end]

    final_state = Cartesian(StiefelScheifele(sol.u[end]), μ, sol.t[end], stische_config)

    expected_end = [
        29447.82922962221
        21027.31807281826
        -1675.1455651936633
        0.15486337817570853
        2.381456403763686
        0.1401977642799883
    ]
    @test final_state ≈ expected_end rtol=1e-3
end
