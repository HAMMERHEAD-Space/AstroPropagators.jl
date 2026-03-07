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

@testset "Cowell Propagator High-Fidelity Regression" begin
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

    satellite_srp_model = CannonballFixedSRP(0.5)
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

    tspan = (0.0, 3 * 86400.0)

    EOM!(du, u, p, t) = Cowell_EOM!(du, u, p, t, model_list)

    prob = ODEProblem(EOM!, u0, tspan, p)
    sol = solve(prob, Vern9(); abstol=1e-13, reltol=1e-13)

    expected_end = [
        -6450.188981590596
        -2390.1583151356413
        500.3061454352241
        4.812545945042311
        -7.890468423940125
        -1.070193288763016
    ]
    @test sol.u[end] ≈ expected_end rtol=1e-3
end
