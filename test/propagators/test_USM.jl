@testset "USM7 Propagator Keplerian" begin
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

    model_list = CentralBodyDynamicsModel(grav_model)
    tspan = (0.0, 86400.0)

    EOM!(du, u, p, t) = USM7_EOM!(du, u, p, t, model_list)

    u0_USM7 = Array(USM7(Cartesian(u0), p.μ))

    prob = ODEProblem(EOM!, u0_USM7, tspan, p)
    sol = solve(prob, Vern9(); abstol=1e-13, reltol=1e-13)

    NRG = orbitalNRG.(USM7.(sol.u), p.μ)

    @test NRG[1] ≈ NRG[end]

    h = norm.(angularMomentumVector.(USM7.(sol.u), p.μ))
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
    @test Cartesian(USM7(sol.u[end]), p.μ) ≈ expected_end
end

@testset "USM6 Propagator Keplerian" begin
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

    model_list = CentralBodyDynamicsModel(grav_model)
    tspan = (0.0, 86400.0)

    EOM!(du, u, p, t) = USM6_EOM!(du, u, p, t, model_list)

    u0_USM6 = Array(USM6(Cartesian(u0), p.μ))

    prob = ODEProblem(EOM!, u0_USM6, tspan, p)
    sol = solve(prob, Vern9(); abstol=1e-13, reltol=1e-13)

    NRG = orbitalNRG.(USM6.(sol.u), p.μ)

    @test NRG[1] ≈ NRG[end]

    h = norm.(angularMomentumVector.(USM6.(sol.u), p.μ))
    @test h[1] ≈ h[end]

    expected_end = [
        29447.829229065504
        21027.31807433234
        -1675.1455650359862
        0.1548633780130051
        2.3814564036944668
        0.1401977642923555
    ]
    @test Cartesian(USM6(sol.u[end]), p.μ) ≈ expected_end rtol=1e-3
end

@testset "USMEM Propagator Keplerian" begin
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

    model_list = CentralBodyDynamicsModel(grav_model)
    tspan = (0.0, 86400.0)

    EOM!(du, u, p, t) = USMEM_EOM!(du, u, p, t, model_list)

    u0_USMEM = Array(USMEM(Cartesian(u0), p.μ))

    prob = ODEProblem(EOM!, u0_USMEM, tspan, p)
    sol = solve(prob, Vern9(); abstol=1e-13, reltol=1e-13)

    NRG = orbitalNRG.(USMEM.(sol.u), p.μ)

    @test NRG[1] ≈ NRG[end]

    h = norm.(angularMomentumVector.(USMEM.(sol.u), p.μ))
    @test h[1] ≈ h[end]

    #Comparison Against Cowell
    expected_end = [
        29447.829229065504
        21027.31807433234
        -1675.1455650359862
        0.1548633780130051
        2.3814564036944668
        0.1401977642923555
    ]
    @test Cartesian(USMEM(sol.u[end]), p.μ) ≈ expected_end rtol=1e-3
end
