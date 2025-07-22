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

    model_list = CentralBodyDynamicsModel(grav_model)
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

    model_list = CentralBodyDynamicsModel(grav_model)
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

    model_list = CentralBodyDynamicsModel(grav_model)
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

    model_list = CentralBodyDynamicsModel(grav_model)
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

    model_list = CentralBodyDynamicsModel(grav_model)
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

    model_list = CentralBodyDynamicsModel(grav_model)
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

@testset "EDromo Propagator Keplerian with Maneuver (Physical Time)" begin
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

    W = (
        potential(Cartesian(u0), p, 0.0, grav_model) -
        potential(Cartesian(u0), p, 0.0, KeplerianGravityAstroModel(μ=p.μ))
    )
    edromo_config = RegularizedCoordinateConfig(
        u0, p.μ; W=W, t₀=0.0, flag_time=PhysicalTime()
    )

    ϕ₀ = compute_initial_phi(u0, p.μ, edromo_config)

    model_list = CentralBodyDynamicsModel(grav_model)
    tspan = (ϕ₀, ϕ₀ + 200π)

    deltaV = [0.05; 0.01; 0.01]

    burn_callback = EDromo_burn(43200.0, deltaV, edromo_config)
    end_callback = end_EDromo_integration(86400.0, edromo_config)

    callback_set = CallbackSet(burn_callback, end_callback)

    EOM!(du, u, p, t) = EDromo_EOM!(du, u, p, t, model_list, edromo_config)

    u0_edromo = Array(EDromo(Cartesian(u0), p.μ, ϕ₀, edromo_config))

    prob = ODEProblem(EOM!, u0_edromo, tspan, p)
    sol = solve(prob, VCABM(); abstol=1e-13, reltol=1e-13, callback=callback_set)

    NRG = zeros(length(sol.u))
    for i in 1:length(sol.u)
        NRG[i] = orbitalNRG(EDromo(sol.u[i]), grav_model.μ, sol.t[i], edromo_config)
    end

    states = zeros(6, length(sol.u))
    for i in 1:length(sol.u)
        states[:, i] = Array(Cartesian(EDromo(sol.u[i]), p.μ, sol.t[i], edromo_config))
    end

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

    final_state = Cartesian(EDromo(sol.u[end]), p.μ, sol.t[end], edromo_config)

    @test final_state ≈ expected_end rtol = 1e-4
end

@testset "EDromo Propagator Keplerian with Maneuver (Constant Time)" begin
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

    W = (
        potential(Cartesian(u0), p, 0.0, grav_model) -
        potential(Cartesian(u0), p, 0.0, KeplerianGravityAstroModel(μ=p.μ))
    )
    edromo_config = RegularizedCoordinateConfig(
        u0, p.μ; W=W, t₀=0.0, flag_time=ConstantTime()
    )

    ϕ₀ = compute_initial_phi(u0, p.μ, edromo_config)

    model_list = CentralBodyDynamicsModel(grav_model)
    tspan = (ϕ₀, ϕ₀ + 200π)

    deltaV = [0.05; 0.01; 0.01]

    burn_callback = EDromo_burn(43200.0, deltaV, edromo_config)
    end_callback = end_EDromo_integration(86400.0, edromo_config)

    callback_set = CallbackSet(burn_callback, end_callback)

    EOM!(du, u, p, t) = EDromo_EOM!(du, u, p, t, model_list, edromo_config)

    u0_edromo = Array(EDromo(Cartesian(u0), p.μ, ϕ₀, edromo_config))

    prob = ODEProblem(EOM!, u0_edromo, tspan, p)
    sol = solve(prob, VCABM(); abstol=1e-13, reltol=1e-13, callback=callback_set)

    NRG = zeros(length(sol.u))
    for i in 1:length(sol.u)
        NRG[i] = orbitalNRG(EDromo(sol.u[i]), grav_model.μ, sol.t[i], edromo_config)
    end

    states = zeros(6, length(sol.u))
    for i in 1:length(sol.u)
        states[:, i] = Array(Cartesian(EDromo(sol.u[i]), p.μ, sol.t[i], edromo_config))
    end

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

    final_state = Cartesian(EDromo(sol.u[end]), p.μ, sol.t[end], edromo_config)

    @test final_state ≈ expected_end rtol = 1e-4
end

@testset "EDromo Propagator Keplerian with Maneuver (Linear Time)" begin
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

    W = (
        potential(Cartesian(u0), p, 0.0, grav_model) -
        potential(Cartesian(u0), p, 0.0, KeplerianGravityAstroModel(μ=p.μ))
    )
    edromo_config = RegularizedCoordinateConfig(
        u0, p.μ; W=W, t₀=0.0, flag_time=LinearTime()
    )

    ϕ₀ = compute_initial_phi(u0, p.μ, edromo_config)

    model_list = CentralBodyDynamicsModel(grav_model)
    tspan = (ϕ₀, ϕ₀ + 200π)

    deltaV = [0.05; 0.01; 0.01]

    burn_callback = EDromo_burn(43200.0, deltaV, edromo_config)
    end_callback = end_EDromo_integration(86400.0, edromo_config)

    callback_set = CallbackSet(burn_callback, end_callback)

    EOM!(du, u, p, t) = EDromo_EOM!(du, u, p, t, model_list, edromo_config)

    u0_edromo = Array(EDromo(Cartesian(u0), p.μ, ϕ₀, edromo_config))

    prob = ODEProblem(EOM!, u0_edromo, tspan, p)
    sol = solve(prob, VCABM(); abstol=1e-13, reltol=1e-13, callback=callback_set)

    NRG = zeros(length(sol.u))
    for i in 1:length(sol.u)
        NRG[i] = orbitalNRG(EDromo(sol.u[i]), grav_model.μ, sol.t[i], edromo_config)
    end

    states = zeros(6, length(sol.u))
    for i in 1:length(sol.u)
        states[:, i] = Array(Cartesian(EDromo(sol.u[i]), p.μ, sol.t[i], edromo_config))
    end

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

    final_state = Cartesian(EDromo(sol.u[end]), p.μ, sol.t[end], edromo_config)

    @test final_state ≈ expected_end rtol = 3e-1
end

@testset "Kustaanheimo-Stiefel Propagator Keplerian with Maneuver (Physical Time)" begin
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

    W = (
        potential(Cartesian(u0), p, 0.0, grav_model) -
        potential(Cartesian(u0), p, 0.0, KeplerianGravityAstroModel(μ=p.μ))
    )
    ks_config = RegularizedCoordinateConfig(u0, p.μ; W=W, t₀=0.0, flag_time=PhysicalTime())

    model_list = CentralBodyDynamicsModel(grav_model)
    tspan = (0.0, 9π)

    deltaV = [0.05; 0.01; 0.01]

    burn_callback = KS_burn(43200.0, deltaV, ks_config)
    end_callback = end_KS_integration(86400.0, ks_config)

    callback_set = CallbackSet(burn_callback, end_callback)

    EOM!(du, u, p, t) = KS_EOM!(du, u, p, t, model_list, ks_config)

    u0_ks = Array(KustaanheimoStiefel(Cartesian(u0), p.μ, ks_config))

    prob = ODEProblem(EOM!, u0_ks, tspan, p)
    sol = solve(prob, VCABM(); abstol=1e-13, reltol=1e-13, callback=callback_set)

    NRG = orbitalNRG.(KustaanheimoStiefel.(sol.u), grav_model.μ, [ks_config])

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

    final_state = Cartesian(KustaanheimoStiefel(sol.u[end]), p.μ, ks_config)

    @test final_state ≈ expected_end rtol = 1e-4
end

@testset "Kustaanheimo-Stiefel Propagator Keplerian with Maneuver (Linear Time)" begin
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

    W = (
        potential(Cartesian(u0), p, 0.0, grav_model) -
        potential(Cartesian(u0), p, 0.0, KeplerianGravityAstroModel(μ=p.μ))
    )
    ks_config = RegularizedCoordinateConfig(u0, p.μ; W=W, t₀=0.0, flag_time=LinearTime())

    model_list = CentralBodyDynamicsModel(grav_model)
    tspan = (0.0, 9π)

    deltaV = [0.05; 0.01; 0.01]

    burn_callback = KS_burn(43200.0, deltaV, ks_config)
    end_callback = end_KS_integration(86400.0, ks_config)

    callback_set = CallbackSet(burn_callback, end_callback)

    EOM!(du, u, p, t) = KS_EOM!(du, u, p, t, model_list, ks_config)

    u0_ks = Array(KustaanheimoStiefel(Cartesian(u0), p.μ, ks_config))

    prob = ODEProblem(EOM!, u0_ks, tspan, p)
    sol = solve(prob, VCABM(); abstol=1e-13, reltol=1e-13, callback=callback_set)

    NRG = orbitalNRG.(KustaanheimoStiefel.(sol.u), grav_model.μ, [ks_config])

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

    final_state = Cartesian(KustaanheimoStiefel(sol.u[end]), p.μ, ks_config)

    @test final_state ≈ expected_end rtol = 1e-4
end
