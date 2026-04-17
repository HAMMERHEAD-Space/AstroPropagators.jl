@testset "Cowell Propagator Keplerian" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
    p = ComponentVector(; JD=JD)

    grav_model = KeplerianGravityAstroModel()

    u0 = [
        -1076.225324679696
        -6765.896364327722
        -332.3087833503755
        9.356857417032581
        -3.3123476319597557
        -1.1880157328553503
    ] #km, km/s

    model_list = CentralBodyDynamicsModel(grav_model)
    tspan = (0.0, 86400.0)

    EOM!(du, u, p, t) = Cowell_EOM!(du, u, p, t, model_list)

    prob = ODEProblem(EOM!, u0, tspan, p)
    sol = solve(prob, Vern9(); abstol=1e-13, reltol=1e-13)

    NRG = orbitalNRG.(Cartesian.(sol.u), grav_model.μ)

    @test NRG[1] ≈ NRG[end]
    h = angularMomentumQuantity.(sol.u)
    @test h[1] ≈ h[end]
    sol.u[end]
    # Regression Test
    expected_end = [
        29447.829229065504
        21027.31807433234
        -1675.1455650359862
        0.1548633780130051
        2.3814564036944668
        0.1401977642923555
    ]
    @test sol.u[end] ≈ expected_end
end
