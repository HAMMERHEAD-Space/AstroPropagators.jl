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

@testset "GEqOE Propagator High-Fidelity" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)

    SpaceIndices.init()
    eop_data = fetch_iers_eop()
    grav_coeffs = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96))

    grav_model = GravityHarmonicsAstroModel(;
        gravity_model=grav_coeffs, eop_data=eop_data, order=36, degree=36
    )
    μ = GravityModels.gravity_constant(grav_model.gravity_model) / 1E9
    p = ComponentVector(; JD=JD, μ=μ)

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

    W =
        potential(Cartesian(u0), p, 0.0, grav_model) -
        potential(Cartesian(u0), p, 0.0, KeplerianGravityAstroModel(; μ=μ))

    config = RegularizedCoordinateConfig(; W=W)
    u0_geqoe = Array(GEqOE(Cartesian(u0), μ, config))

    tspan = (0.0, 86400.0)

    EOM!(du, u, p, t) = GEqOE_EOM!(du, u, p, t, model_list, config)

    prob = ODEProblem(EOM!, u0_geqoe, tspan, p)
    sol = solve(prob, Vern9(); abstol=1e-13, reltol=1e-13)

    expected_end = [
        29209.16404907953
        22221.199335560723
        -1539.7320425979071
        0.020138201496128487
        2.3045214269873857
        0.15104845625911167
    ]
    @test Cartesian(GEqOE(sol.u[end]), μ, config) ≈ expected_end rtol = 1e-2
end

@testset "GEqOE Propagator High-Fidelity Regression" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)

    SpaceIndices.init()
    eop_data = fetch_iers_eop()
    grav_coeffs = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96))

    grav_model = GravityHarmonicsAstroModel(;
        gravity_model=grav_coeffs, eop_data=eop_data, order=36, degree=36
    )
    μ = GravityModels.gravity_constant(grav_model.gravity_model) / 1E9
    p = ComponentVector(; JD=JD, μ=μ)

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
        8.956857417032581
        -3.3123476319597557
        -1.1880157328553503
    ] #km, km/s

    model_list = CentralBodyDynamicsModel(
        grav_model, (sun_third_body, moon_third_body, srp_model, drag_model)
    )

    W =
        potential(Cartesian(u0), p, 0.0, grav_model) -
        potential(Cartesian(u0), p, 0.0, KeplerianGravityAstroModel(; μ=μ))

    config = RegularizedCoordinateConfig(; W=W)
    u0_geqoe = Array(GEqOE(Cartesian(u0), μ, config))

    tspan = (0.0, 3 * 86400.0)

    EOM!(du, u, p, t) = GEqOE_EOM!(du, u, p, t, model_list, config)

    prob = ODEProblem(EOM!, u0_geqoe, tspan, p)
    sol = solve(prob, Vern9(); abstol=1e-13, reltol=1e-13)

    expected_end = [
        -6462.555199025645
        -2369.8382849120076
        503.0595947121262
        4.792369569779785
        -7.897922283599572
        -1.06862260690453
    ]
    @test Cartesian(GEqOE(sol.u[end]), μ, config) ≈ expected_end rtol = 2e0
end
