using AstroMeasurements

@testset "AstroPropagatorsMeasurementsExt" begin
    setup = setup_keplerian()
    (; model_list, p, μ) = setup

    u0_cart = [
        -1076.225324679696
        -6765.896364327722
        -332.3087833503755
        9.356857417032581
        -3.3123476319597557
        -1.1880157328553503
    ]

    eop_data = fetch_iers_eop()

    # ================================================================
    # elevation_condition — basic construction and evaluation
    # ================================================================
    @testset "elevation_condition" begin
        gs = RadarGroundStation("TestGS", deg2rad(38.9), deg2rad(-77.0), 0.0)

        cond = elevation_condition(Cartesian, gs, eop_data)
        @test cond isa Function

        mock_integrator = (p=p,)
        g = cond(u0_cart, 0.0, mock_integrator)
        @test g isa Number

        el_ref = AstroMeasurements.compute_elevation(
            SVector{6}(u0_cart...), gs, p.JD, eop_data
        )
        @test g ≈ el_ref - deg2rad(10.0) rtol = 1e-12
    end

    @testset "elevation_condition (custom el_min)" begin
        gs = RadarGroundStation("TestGS", deg2rad(38.9), deg2rad(-77.0), 0.0)

        cond5 = elevation_condition(Cartesian, gs, eop_data; el_min=deg2rad(5.0))
        cond20 = elevation_condition(Cartesian, gs, eop_data; el_min=deg2rad(20.0))

        mock_integrator = (p=p,)
        g5 = cond5(u0_cart, 0.0, mock_integrator)
        g20 = cond20(u0_cart, 0.0, mock_integrator)

        @test g5 - g20 ≈ deg2rad(20.0 - 5.0) rtol = 1e-12
    end

    @testset "elevation_condition (integration)" begin
        gs = RadarGroundStation("TestGS", deg2rad(38.9), deg2rad(-77.0), 0.0)

        cond = elevation_condition(Cartesian, gs, eop_data; el_min=deg2rad(5.0))

        tspan = (0.0, 86400.0)
        EOM!(du, u, _p, t) = Cowell_EOM!(du, u, _p, t, model_list)

        passes = []
        cb = ContinuousCallback(
            cond,
            (integrator) -> push!(passes, (t=integrator.t, event=:aos)),
            (integrator) -> push!(passes, (t=integrator.t, event=:los)),
        )

        prob = ODEProblem(EOM!, u0_cart, tspan, p)
        sol = solve(prob, Vern9(); abstol=1e-13, reltol=1e-13, callback=cb)

        @test sol.retcode == ReturnCode.Success

        if length(passes) >= 2
            aos_events = filter(e -> e.event == :aos, passes)
            for ev in aos_events
                state_at_event = sol(ev.t)
                el = AstroMeasurements.compute_elevation(
                    SVector{6}(state_at_event...), gs, p.JD + ev.t / 86400.0, eop_data
                )
                @test abs(el - deg2rad(5.0)) < 1e-4
            end
        end
    end

    # ================================================================
    # relative_distance_condition
    # ================================================================
    @testset "relative_distance_condition" begin
        target_pos = SVector{3}(7000.0, 0.0, 0.0)
        target_fn = (JD) -> target_pos

        d_target = 10000.0
        cond = relative_distance_condition(Cartesian, target_fn, d_target)
        @test cond isa Function

        mock_integrator = (p=p,)
        g = cond(u0_cart, 0.0, mock_integrator)

        r_sc = SVector{3}(u0_cart[1], u0_cart[2], u0_cart[3])
        expected = norm(r_sc - target_pos) - d_target
        @test g ≈ expected rtol = 1e-12
    end

    @testset "relative_distance_condition (integration)" begin
        target_pos = SVector{3}(7000.0, 0.0, 0.0)
        target_fn = (JD) -> target_pos

        cond = relative_distance_condition(Cartesian, target_fn, 8000.0)

        tspan = (0.0, 86400.0)
        EOM!(du, u, _p, t) = Cowell_EOM!(du, u, _p, t, model_list)

        crossing_log = []
        cb = ContinuousCallback(
            cond,
            (integrator) -> push!(crossing_log, integrator.t),
            (integrator) -> push!(crossing_log, integrator.t),
        )

        prob = ODEProblem(EOM!, u0_cart, tspan, p)
        sol = solve(prob, Vern9(); abstol=1e-13, reltol=1e-13, callback=cb)

        @test sol.retcode == ReturnCode.Success

        for t_ev in crossing_log
            state = sol(t_ev)
            r_sc = SVector{3}(state[1], state[2], state[3])
            d = norm(r_sc - target_pos)
            @test abs(d - 8000.0) < 1e-4
        end
    end

    # ================================================================
    # angular_separation_condition
    # ================================================================
    @testset "angular_separation_condition" begin
        beacon_fn = (JD) -> SVector{3}(1.496e8, 0.0, 0.0)
        observer_fn = (JD) -> SVector{3}(6378.0, 0.0, 0.0)
        threshold = deg2rad(30.0)

        cond = angular_separation_condition(Cartesian, beacon_fn, observer_fn, threshold)
        @test cond isa Function

        mock_integrator = (p=p,)
        g = cond(u0_cart, 0.0, mock_integrator)
        @test g isa Number

        r_sc = SVector{3}(u0_cart[1], u0_cart[2], u0_cart[3])
        r_beacon = SVector{3}(1.496e8, 0.0, 0.0)
        r_observer = SVector{3}(6378.0, 0.0, 0.0)
        d_sc = r_sc - r_observer
        d_beacon = r_beacon - r_observer
        cosθ = dot(d_sc, d_beacon) / (norm(d_sc) * norm(d_beacon))
        cosθ = clamp(cosθ, -1.0, 1.0)
        expected = acos(cosθ) - threshold
        @test g ≈ expected rtol = 1e-10
    end

    @testset "angular_separation_condition (integration)" begin
        beacon_fn = (JD) -> SVector{3}(1.496e8, 0.0, 0.0)
        observer_fn = (JD) -> SVector{3}(6378.0, 0.0, 0.0)

        cond = angular_separation_condition(
            Cartesian, beacon_fn, observer_fn, deg2rad(90.0)
        )

        tspan = (0.0, 86400.0)
        EOM!(du, u, _p, t) = Cowell_EOM!(du, u, _p, t, model_list)

        crossing_log = []
        cb = ContinuousCallback(
            cond,
            (integrator) -> push!(crossing_log, integrator.t),
            (integrator) -> push!(crossing_log, integrator.t),
        )

        prob = ODEProblem(EOM!, u0_cart, tspan, p)
        sol = solve(prob, Vern9(); abstol=1e-13, reltol=1e-13, callback=cb)

        @test sol.retcode == ReturnCode.Success
    end

    # ================================================================
    # Extension with non-Cartesian coordinates (Keplerian)
    # ================================================================
    @testset "Extension with Keplerian coords" begin
        target_fn = (JD) -> SVector{3}(7000.0, 0.0, 0.0)

        cond_cart = relative_distance_condition(Cartesian, target_fn, 8000.0)
        cond_kep = relative_distance_condition(Keplerian, target_fn, 8000.0)

        mock_integrator = (p=p,)
        g_cart = cond_cart(u0_cart, 0.0, mock_integrator)

        kep_u0 = Array(Keplerian(Cartesian(u0_cart), μ))
        g_kep = cond_kep(kep_u0, 0.0, mock_integrator)

        @test g_cart ≈ g_kep rtol = 1e-8
    end

    # ================================================================
    # Extension with utility combinators
    # ================================================================
    @testset "Extension with and_condition" begin
        gs = RadarGroundStation("TestGS", deg2rad(38.9), deg2rad(-77.0), 0.0)

        el_cond = elevation_condition(Cartesian, gs, eop_data)
        alt_cond = altitude_condition(Cartesian, 300.0)

        combined = and_condition(el_cond, alt_cond)

        mock_integrator = (p=p,)
        g = combined(u0_cart, 0.0, mock_integrator)

        g_el = el_cond(u0_cart, 0.0, mock_integrator)
        g_alt = alt_cond(u0_cart, 0.0, mock_integrator)
        @test g ≈ min(g_el, g_alt)
    end
end
