@testset "GEqOE Propagator Keplerian" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)

    grav_model = KeplerianGravityAstroModel()
    μ = grav_model.μ
    p = ComponentVector(; JD=JD, μ=μ)

    u0 = [
        -1076.225324679696
        -6765.896364327722
        -332.3087833503755
        9.356857417032581
        -3.3123476319597557
        -1.1880157328553503
    ] #km, km/s

    W =
        potential(Cartesian(u0), p, 0.0, grav_model) -
        potential(Cartesian(u0), p, 0.0, KeplerianGravityAstroModel(; μ=μ))

    config = RegularizedCoordinateConfig(; W=W)
    u0_geqoe = Array(GEqOE(Cartesian(u0), μ, config))

    model_list = CentralBodyDynamicsModel(grav_model)
    tspan = (0.0, 86400.0)

    EOM!(du, u, p, t) = GEqOE_EOM!(du, u, p, t, model_list, config)

    prob = ODEProblem(EOM!, u0_geqoe, tspan, p)
    sol = solve(prob, Vern9(); abstol=1e-13, reltol=1e-13)

    cart_first = Array(Cartesian(GEqOE(sol.u[1]), μ, config))
    cart_last = Array(Cartesian(GEqOE(sol.u[end]), μ, config))
    NRG_first = 0.5 * sum(abs2, cart_first[4:6]) - μ / sqrt(sum(abs2, cart_first[1:3]))
    NRG_last = 0.5 * sum(abs2, cart_last[4:6]) - μ / sqrt(sum(abs2, cart_last[1:3]))
    @test NRG_first ≈ NRG_last

    expected_end = [
        29447.829229065504
        21027.31807433234
        -1675.1455650359862
        0.1548633780130051
        2.3814564036944668
        0.1401977642923555
    ]
    @test Cartesian(GEqOE(sol.u[end]), μ, config) ≈ expected_end
end
