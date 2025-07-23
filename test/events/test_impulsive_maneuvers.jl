@testset "Cowell Propagator Keplerian with Maneuver" begin
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

    model_list = (grav_model,)
    tspan = (0.0, 86400.0)

    deltaV = [0.05; 0.01; 0.01]
    function condition(
        u::AbstractVector, t::Number, integrator::T
    ) where {T<:SciMLBase.DEIntegrator}
        return t - 43200.0
    end

    burn = ContinuousCallback(
        condition, (integrator) -> impulsive_burn!(integrator, deltaV)
    )

    EOM!(du, u, p, t) = Cowell_EOM!(du, u, p, t, model_list)

    prob = ODEProblem(EOM!, u0, tspan, p)
    sol = solve(prob, VCABM(); abstol=1e-13, reltol=1e-13, callback=burn)

    NRG = orbitalNRG.(Cartesian.(sol.u), grav_model.μ)

    @test NRG[1] != NRG[end]

    # Regression Test
    expected_end = [
        29390.280395821836
        18637.945967159154
        -1768.361355756133
        0.47323343997331674
        2.572684107343496
        0.13273831002165992
    ]
    @test sol.u[end] ≈ expected_end
end

@testset "Gauss Variational Propagator Keplerian with Maneuver" begin
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

    model_list = (grav_model,)
    tspan = (0.0, 86400.0)

    deltaV = [0.05; 0.01; 0.01]
    function condition(
        u::AbstractVector, t::Number, integrator::T
    ) where {T<:SciMLBase.DEIntegrator}
        return t - 43200.0
    end

    burn = ContinuousCallback(
        condition,
        (integrator) -> impulsive_burn!(integrator, deltaV; coordinate_set=Keplerian),
    )

    EOM!(du, u, p, t) = GaussVE_EOM!(du, u, p, t, model_list)

    u0_koe = Array(Keplerian(Cartesian(u0), p.μ))

    prob = ODEProblem(EOM!, u0_koe, tspan, p)
    sol = solve(prob, VCABM(); abstol=1e-13, reltol=1e-13, callback=burn)

    NRG = orbitalNRG.(Keplerian.(sol.u), grav_model.μ)

    @test NRG[1] != NRG[end]

    # Comparison with Cowell
    expected_end = [
        29390.280395821836
        18637.945967159154
        -1768.361355756133
        0.47323343997331674
        2.572684107343496
        0.13273831002165992
    ]
    @test Cartesian(Keplerian(sol.u[end]), p.μ) ≈ expected_end rtol = 1e-4
end

@testset "Milankovich Propagator Keplerian with Maneuver" begin
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

    model_list = (grav_model,)
    tspan = (0.0, 86400.0)

    deltaV = [0.05; 0.01; 0.01]
    function condition(
        u::AbstractVector, t::Number, integrator::T
    ) where {T<:SciMLBase.DEIntegrator}
        return t - 43200.0
    end

    burn = ContinuousCallback(
        condition,
        (integrator) -> impulsive_burn!(integrator, deltaV; coordinate_set=Milankovich),
    )

    EOM!(du, u, p, t) = Milankovich_EOM!(du, u, p, t, model_list)

    u0_mil = Array(Milankovich(Cartesian(u0), p.μ))

    prob = ODEProblem(EOM!, u0_mil, tspan, p)
    sol = solve(prob, VCABM(); abstol=1e-13, reltol=1e-13, callback=burn)

    NRG = orbitalNRG.(Milankovich.(sol.u), grav_model.μ)

    @test NRG[1] != NRG[end]

    # Comparison with Cowell
    expected_end = [
        29390.280395821836
        18637.945967159154
        -1768.361355756133
        0.47323343997331674
        2.572684107343496
        0.13273831002165992
    ]
    @test Cartesian(Milankovich(sol.u[end]), p.μ) ≈ expected_end rtol = 1e-2
end

@testset "USM7 Propagator Keplerian with Maneuver" begin
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

    model_list = (grav_model,)
    tspan = (0.0, 86400.0)

    deltaV = [0.05; 0.01; 0.01]
    function condition(
        u::AbstractVector, t::Number, integrator::T
    ) where {T<:SciMLBase.DEIntegrator}
        return t - 43200.0
    end

    burn = ContinuousCallback(
        condition, (integrator) -> impulsive_burn!(integrator, deltaV; coordinate_set=USM7)
    )

    EOM!(du, u, p, t) = USM7_EOM!(du, u, p, t, model_list)

    u0_mil = Array(USM7(Cartesian(u0), p.μ))

    prob = ODEProblem(EOM!, u0_mil, tspan, p)
    sol = solve(prob, VCABM(); abstol=1e-13, reltol=1e-13, callback=burn)

    NRG = orbitalNRG.(USM7.(sol.u), grav_model.μ)

    @test NRG[1] != NRG[end]

    # Comparison with Cowell
    expected_end = [
        29390.280395821836
        18637.945967159154
        -1768.361355756133
        0.47323343997331674
        2.572684107343496
        0.13273831002165992
    ]
    @test Cartesian(USM7(sol.u[end]), p.μ) ≈ expected_end rtol = 1e-4
end

@testset "USM6 Propagator Keplerian with Maneuver" begin
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

    model_list = (grav_model,)
    tspan = (0.0, 86400.0)

    deltaV = [0.05; 0.01; 0.01]
    function condition(
        u::AbstractVector, t::Number, integrator::T
    ) where {T<:SciMLBase.DEIntegrator}
        return t - 43200.0
    end

    burn = ContinuousCallback(
        condition, (integrator) -> impulsive_burn!(integrator, deltaV; coordinate_set=USM6)
    )

    EOM!(du, u, p, t) = USM6_EOM!(du, u, p, t, model_list)

    u0_mil = Array(USM6(Cartesian(u0), p.μ))

    prob = ODEProblem(EOM!, u0_mil, tspan, p)
    sol = solve(prob, VCABM(); abstol=1e-13, reltol=1e-13, callback=burn)

    NRG = orbitalNRG.(USM6.(sol.u), grav_model.μ)

    @test NRG[1] != NRG[end]

    # Comparison with Cowell
    expected_end = [
        29390.280395821836
        18637.945967159154
        -1768.361355756133
        0.47323343997331674
        2.572684107343496
        0.13273831002165992
    ]
    @test Cartesian(USM6(sol.u[end]), p.μ) ≈ expected_end rtol = 1e-4
end

@testset "USMEM Propagator Keplerian with Maneuver" begin
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

    model_list = (grav_model,)
    tspan = (0.0, 86400.0)

    deltaV = [0.05; 0.01; 0.01]
    function condition(
        u::AbstractVector, t::Number, integrator::T
    ) where {T<:SciMLBase.DEIntegrator}
        return t - 43200.0
    end

    burn = ContinuousCallback(
        condition, (integrator) -> impulsive_burn!(integrator, deltaV; coordinate_set=USMEM)
    )

    EOM!(du, u, p, t) = USMEM_EOM!(du, u, p, t, model_list)

    u0_mil = Array(USMEM(Cartesian(u0), p.μ))

    prob = ODEProblem(EOM!, u0_mil, tspan, p)
    sol = solve(prob, VCABM(); abstol=1e-13, reltol=1e-13, callback=burn)

    NRG = orbitalNRG.(USMEM.(sol.u), grav_model.μ)

    @test NRG[1] != NRG[end]

    # Comparison with Cowell
    expected_end = [
        29390.280395821836
        18637.945967159154
        -1768.361355756133
        0.47323343997331674
        2.572684107343496
        0.13273831002165992
    ]
    @test Cartesian(USMEM(sol.u[end]), p.μ) ≈ expected_end rtol = 1e-4
end
