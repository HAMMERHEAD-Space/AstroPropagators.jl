@testset "Modified Equinoctial Propagator Keplerian" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)

    grav_model = KeplerianGravityAstroModel()
    p = ComponentVector(; JD=JD, μ=grav_model.μ)

    u0 = [
        -1076.225324679696
        -6765.896364327722
        -332.3087833503755
        9.356857417032581
        -3.3123476319597557
        -1.1880157328553503
    ] #km, km/s

    u0_meq = Array(ModEq(Cartesian(u0), p.μ))

    model_list = CentralBodyDynamicsModel(grav_model)
    tspan = (0.0, 86400.0)

    EOM!(du, u, p, t) = ModEq_EOM!(du, u, p, t, model_list)

    prob = ODEProblem(EOM!, u0_meq, tspan, p)
    sol = solve(prob, Vern9(); abstol=1e-13, reltol=1e-13)

    NRG = orbitalNRG.(ModEq.(sol.u), p.μ)

    @test NRG[1] ≈ NRG[end]
    h = norm.(angularMomentumVector.(ModEq.(sol.u), p.μ))
    @test h[1] ≈ h[end]

    # Comparison Against Cowell
    expected_end = [
        29447.829229065504
        21027.31807433234
        -1675.1455650359862
        0.1548633780130051
        2.3814564036944668
        0.1401977642923555
    ]
    @test Cartesian(ModEq(sol.u[end]), p.μ) ≈ expected_end
end
