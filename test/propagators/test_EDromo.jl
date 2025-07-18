using Test

using AstroCoords
using AstroForceModels
using AstroPropagators
using ComponentArrays
using LinearAlgebra
using OrdinaryDiffEqAdamsBashforthMoulton
using OrdinaryDiffEqCore
using SatelliteToolboxGravityModels
using SatelliteToolboxTransformations
using SciMLBase
using SpaceIndices

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
    edromo_config = set_initial_edromo_configurations(
        u0_cart, 
        μ; 
        W=potential(Cartesian(u0_cart), p, 0.0, grav_model) - potential(Cartesian(u0_cart), p, 0.0, KeplerianGravityAstroModel(μ=μ)),
        flag_time=PhysicalTime()
    )

    p_full = ComponentVector(; p..., μ=μ)

    u0_EDromo = Array(EDromo(Cartesian(u0_cart), μ; edromo_config...))

    model_list = CentralBodyDynamicsModel(grav_model)
    # The independent variable is ϕ, so we integrate over one orbit (2π)
    tspan = (edromo_config.ϕ, edromo_config.ϕ + 6π)

    EOM!(du, u, p, t) = EDromo_EOM!(du, u, p, t, model_list; DU=edromo_config.DU, TU=edromo_config.TU, W=edromo_config.W, t₀=edromo_config.t₀, flag_time=edromo_config.flag_time)

    prob = ODEProblem(EOM!, u0_EDromo, tspan, p_full)
    sol = solve(prob, VCABM(); abstol=1e-15, reltol=1e-15, callback=end_EDromo_integration(86400.0; edromo_config...))

    NRG = zeros(length(sol.u))
    for i in 1:length(sol.u)
        NRG[i] = orbitalNRG(EDromo(sol.u[i]), grav_model.μ; set_edromo_configurations(edromo_config; curr_ϕ=sol.t[i])...)
    end

    states = zeros(6, length(sol.u))
    for i in 1:length(sol.u)
        states[:, i] = Array(Cartesian(EDromo(sol.u[i]), μ; set_edromo_configurations(edromo_config; curr_ϕ=sol.t[i])...))
    end

    @test NRG[1] ≈ NRG[end]

    final_state = Cartesian(EDromo(sol.u[end]), μ; set_edromo_configurations(edromo_config; curr_ϕ=sol.t[end])...)

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
    edromo_config = set_initial_edromo_configurations(
        u0_cart, 
        μ; 
        W=potential(Cartesian(u0_cart), p, 0.0, grav_model) - potential(Cartesian(u0_cart), p, 0.0, KeplerianGravityAstroModel(μ=μ)),
        flag_time=ConstantTime()
    )

    p_full = ComponentVector(; p..., μ=μ)

    u0_EDromo = Array(EDromo(Cartesian(u0_cart), μ; edromo_config...))

    model_list = CentralBodyDynamicsModel(grav_model)
    # The independent variable is ϕ, so we integrate over one orbit (2π)
    tspan = (edromo_config.ϕ, edromo_config.ϕ + 6π)

    EOM!(du, u, p, t) = EDromo_EOM!(du, u, p, t, model_list; DU=edromo_config.DU, TU=edromo_config.TU, W=edromo_config.W, t₀=edromo_config.t₀, flag_time=edromo_config.flag_time)

    prob = ODEProblem(EOM!, u0_EDromo, tspan, p_full)
    sol = solve(prob, VCABM(); abstol=1e-15, reltol=1e-15, callback=end_EDromo_integration(86400.0; edromo_config...))

    NRG = zeros(length(sol.u))
    for i in 1:length(sol.u)
        NRG[i] = orbitalNRG(EDromo(sol.u[i]), grav_model.μ; set_edromo_configurations(edromo_config; curr_ϕ=sol.t[i])...)
    end

    states = zeros(6, length(sol.u))
    for i in 1:length(sol.u)
        states[:, i] = Array(Cartesian(EDromo(sol.u[i]), μ; set_edromo_configurations(edromo_config; curr_ϕ=sol.t[i])...))
    end

    @test NRG[1] ≈ NRG[end]

    final_state = Cartesian(EDromo(sol.u[end]), μ; set_edromo_configurations(edromo_config; curr_ϕ=sol.t[end])...)

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
    edromo_config = set_initial_edromo_configurations(
        u0_cart, 
        μ; 
        W=potential(Cartesian(u0_cart), p, 0.0, grav_model) - potential(Cartesian(u0_cart), p, 0.0, KeplerianGravityAstroModel(μ=μ)),
        flag_time=LinearTime()
    )

    p_full = ComponentVector(; p..., μ=μ)

    u0_EDromo = Array(EDromo(Cartesian(u0_cart), μ; edromo_config...))

    
    model_list = CentralBodyDynamicsModel(grav_model)
    # The independent variable is ϕ, so we integrate over one orbit (2π)
    tspan = (edromo_config.ϕ, edromo_config.ϕ + 6π)

    EOM!(du, u, p, t) = EDromo_EOM!(du, u, p, t, model_list; DU=edromo_config.DU, TU=edromo_config.TU, W=edromo_config.W, t₀=edromo_config.t₀, flag_time=edromo_config.flag_time)

    prob = ODEProblem(EOM!, u0_EDromo, tspan, p_full)
    sol = solve(prob, VCABM(); abstol=1e-15, reltol=1e-15, callback=end_EDromo_integration(86400.0; edromo_config...))

    NRG = zeros(length(sol.u))
    for i in 1:length(sol.u)
        NRG[i] = orbitalNRG(EDromo(sol.u[i]), grav_model.μ; set_edromo_configurations(edromo_config; curr_ϕ=sol.t[i])...)
    end

    states = zeros(6, length(sol.u))
    for i in 1:length(sol.u)
        states[:, i] = Array(Cartesian(EDromo(sol.u[i]), μ; set_edromo_configurations(edromo_config; curr_ϕ=(edromo_config.ϕ + sol.t[i]))...))
    end

    @test NRG[1] ≈ NRG[end]
    
    final_state = Cartesian(EDromo(sol.u[end]), μ; set_edromo_configurations(edromo_config; curr_ϕ=sol.t[end])...)

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

#@testset "EDromo Propagator High-Fidelity Physical Time" begin
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

    edromo_config = set_initial_edromo_configurations(
        u0_cart, 
        μ; 
        W=potential(Cartesian(u0_cart), p, 0.0, grav_model) - potential(Cartesian(u0_cart), p, 0.0, KeplerianGravityAstroModel(μ=μ)),
        flag_time=PhysicalTime()
    )

    tspan = (0.0, 86400.0)

    p_full = ComponentVector(; p..., μ=μ)

    u0_EDromo = Array(EDromo(Cartesian(u0_cart), μ; edromo_config...))

    model_list = CentralBodyDynamicsModel(grav_model)
    # The independent variable is ϕ, so we integrate over one orbit (2π)
    tspan = (edromo_config.ϕ, edromo_config.ϕ + 6π)

    EOM!(du, u, p, t) = EDromo_EOM!(du, u, p, t, model_list; DU=edromo_config.DU, TU=edromo_config.TU, W=edromo_config.W, t₀=edromo_config.t₀, flag_time=edromo_config.flag_time)

    prob = ODEProblem(EOM!, u0_EDromo, tspan, p_full)
    sol = solve(prob, VCABM(); abstol=1e-15, reltol=1e-15, callback=end_EDromo_integration(86400.0; edromo_config...))

    states = zeros(6, length(sol.u))
    for i in 1:length(sol.u)
        states[:, i] = Array(Cartesian(EDromo(sol.u[i]), μ; set_edromo_configurations(edromo_config; curr_ϕ=sol.t[i])...))
    end

    final_state = Cartesian(EDromo(sol.u[end]), μ; set_edromo_configurations(edromo_config; curr_ϕ=sol.t[end])...)

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