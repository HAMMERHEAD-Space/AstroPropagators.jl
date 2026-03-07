using StaticArraysCore

@testset "Event Detectors" begin
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

    kep0 = Keplerian(Cartesian(u0_cart), μ)
    T_orbit = 2π * sqrt(kep0.a^3 / μ)
    r_peri = kep0.a * (1 - kep0.e)
    r_apo = kep0.a * (1 + kep0.e)
    n_orbits_in_day = 86400.0 / T_orbit

    # ================================================================
    # State Adapter
    # ================================================================
    @testset "State Adapter" begin
        @test coord_type(CowellPropagator()) == Cartesian
        @test coord_type(GaussVEPropagator()) == Keplerian
        @test coord_type(EDromoPropagator()) == EDromo
        @test coord_type(KSPropagator()) == KustaanheimoStiefel
        @test coord_type(StiSchePropagator()) == StiefelScheifele
        @test coord_type(GEqOEPropagator()) == GEqOE
        @test coord_type(MilankovichPropagator()) == Milankovich
        @test coord_type(USM7Propagator()) == USM7
        @test coord_type(USM6Propagator()) == USM6
        @test coord_type(USMEMPropagator()) == USMEM

        cart = get_cartesian(u0_cart, 0.0, μ, Cartesian)
        @test cart ≈ u0_cart

        kep_u0 = Array(Keplerian(Cartesian(u0_cart), μ))
        cart2 = get_cartesian(kep_u0, 0.0, μ, Keplerian)
        @test Array(cart2) ≈ u0_cart rtol = 1e-12

        kep_rt = get_keplerian(u0_cart, 0.0, μ, Cartesian)
        @test kep_rt.a ≈ kep0.a rtol = 1e-12
        @test kep_rt.e ≈ kep0.e rtol = 1e-12

        @test get_physical_time(u0_cart, 100.0, Cartesian) == 100.0
        @test get_physical_time(u0_cart, 100.0, Keplerian) == 100.0
    end

    # ================================================================
    # Orbital Detectors — Apside
    # ================================================================
    @testset "Apside Detector (Cowell)" begin
        cond = apside_condition(Cartesian)

        tspan = (0.0, 86400.0)
        EOM!(du, u, _p, t) = Cowell_EOM!(du, u, _p, t, model_list)

        peri_log = []
        apo_log = []
        peri_affect! =
            (integrator) -> begin
                cart = Cartesian(integrator.u)
                r = sqrt(cart[1]^2 + cart[2]^2 + cart[3]^2)
                rdot = (cart[1]*cart[4] + cart[2]*cart[5] + cart[3]*cart[6])
                push!(peri_log, (t=integrator.t, r=r, rdot=rdot))
            end
        apo_affect! =
            (integrator) -> begin
                cart = Cartesian(integrator.u)
                r = sqrt(cart[1]^2 + cart[2]^2 + cart[3]^2)
                rdot = (cart[1]*cart[4] + cart[2]*cart[5] + cart[3]*cart[6])
                push!(apo_log, (t=integrator.t, r=r, rdot=rdot))
            end

        cb = ContinuousCallback(cond, peri_affect!, apo_affect!)

        prob = ODEProblem(EOM!, u0_cart, tspan, p)
        sol = solve(prob, Vern9(); abstol=1e-13, reltol=1e-13, callback=cb)

        n_expected = floor(Int, n_orbits_in_day)
        @test length(peri_log) >= n_expected
        @test length(apo_log) >= n_expected

        for ev in peri_log
            @test abs(ev.rdot) < 1e-8
            @test ev.r ≈ r_peri rtol = 1e-6
        end

        for ev in apo_log
            @test abs(ev.rdot) < 1e-8
            @test ev.r ≈ r_apo rtol = 1e-6
        end

        peri_times = [e.t for e in peri_log]
        for i in 2:length(peri_times)
            @test peri_times[i] - peri_times[i - 1] ≈ T_orbit rtol = 1e-6
        end
    end

    # ================================================================
    # Orbital Detectors — Node
    # ================================================================
    @testset "Node Detector (Cowell)" begin
        cond = node_condition(Cartesian)

        tspan = (0.0, 86400.0)
        EOM!(du, u, _p, t) = Cowell_EOM!(du, u, _p, t, model_list)

        asc_log = []
        desc_log = []
        cb = ContinuousCallback(
            cond,
            (integrator) ->
                push!(asc_log, (t=integrator.t, z=integrator.u[3], vz=integrator.u[6])),
            (integrator) ->
                push!(desc_log, (t=integrator.t, z=integrator.u[3], vz=integrator.u[6])),
        )

        prob = ODEProblem(EOM!, u0_cart, tspan, p)
        sol = solve(prob, Vern9(); abstol=1e-13, reltol=1e-13, callback=cb)

        @test length(asc_log) + length(desc_log) >= 2 * floor(Int, n_orbits_in_day)

        for ev in asc_log
            @test abs(ev.z) < 1e-9
            @test ev.vz > 0
        end
        for ev in desc_log
            @test abs(ev.z) < 1e-9
            @test ev.vz < 0
        end
    end

    # ================================================================
    # Orbital Detectors — True Anomaly Crossing
    # ================================================================
    @testset "True Anomaly Crossing (Cowell)" begin
        f_target = π
        cond = true_anomaly_condition(Cartesian, f_target)

        tspan = (0.0, 86400.0)
        EOM!(du, u, _p, t) = Cowell_EOM!(du, u, _p, t, model_list)

        event_log = []
        affect_fn! =
            (integrator) -> begin
                kep = Keplerian(Cartesian(integrator.u), integrator.p.μ)
                push!(event_log, (t=integrator.t, f=kep.f))
            end

        cb = ContinuousCallback(cond, affect_fn!, affect_fn!)

        prob = ODEProblem(EOM!, u0_cart, tspan, p)
        sol = solve(prob, Vern9(); abstol=1e-13, reltol=1e-13, callback=cb)

        @test length(event_log) >= floor(Int, n_orbits_in_day)

        for ev in event_log
            @test abs(sin(ev.f - f_target)) < 1e-10
        end
    end

    # ================================================================
    # Orbital Detectors — Argument of Latitude
    # ================================================================
    @testset "Argument of Latitude Crossing (Cowell)" begin
        u_target = π / 2
        cond = argument_of_latitude_condition(Cartesian, u_target)

        tspan = (0.0, 86400.0)
        EOM!(du, u, _p, t) = Cowell_EOM!(du, u, _p, t, model_list)

        event_log = []
        affect_fn! =
            (integrator) -> begin
                kep = Keplerian(Cartesian(integrator.u), integrator.p.μ)
                push!(event_log, (t=integrator.t, aol=(kep.ω + kep.f)))
            end

        cb = ContinuousCallback(cond, affect_fn!, affect_fn!)

        prob = ODEProblem(EOM!, u0_cart, tspan, p)
        sol = solve(prob, Vern9(); abstol=1e-13, reltol=1e-13, callback=cb)

        @test length(event_log) >= 1

        for ev in event_log
            @test abs(sin(ev.aol - u_target)) < 1e-10
        end
    end

    # ================================================================
    # Orbital Detectors — Mean Anomaly
    # ================================================================
    @testset "Mean Anomaly Crossing (Cowell)" begin
        M_target = 0.0
        cond = mean_anomaly_condition(Cartesian, M_target)

        tspan = (0.0, 86400.0)
        EOM!(du, u, _p, t) = Cowell_EOM!(du, u, _p, t, model_list)

        event_log = []
        affect_fn! =
            (integrator) -> begin
                kep = Keplerian(Cartesian(integrator.u), integrator.p.μ)
                push!(event_log, (t=integrator.t, M=kep.M))
            end

        cb = ContinuousCallback(cond, affect_fn!, affect_fn!)

        prob = ODEProblem(EOM!, u0_cart, tspan, p)
        sol = solve(prob, Vern9(); abstol=1e-13, reltol=1e-13, callback=cb)

        @test length(event_log) >= floor(Int, n_orbits_in_day)

        for ev in event_log
            @test abs(sin(ev.M - M_target)) < 1e-10
        end
    end

    # ================================================================
    # Orbital Detectors — RAAN (sanity check, slow drift in Kepler)
    # ================================================================
    @testset "RAAN Crossing (Cowell)" begin
        Ω0 = kep0.Ω
        Ω_target = Ω0 + 0.001
        cond = raan_condition(Cartesian, Ω_target)

        mock_integrator = (p=p,)
        g_val = cond(u0_cart, 0.0, mock_integrator)
        @test g_val isa Number
        @test abs(g_val - sin(Ω0 - Ω_target)) < 1e-10
    end

    # ================================================================
    # Geometric Detectors — Altitude
    # ================================================================
    @testset "Altitude Detector (Cowell)" begin
        h_target = 400.0
        cond = altitude_condition(Cartesian, h_target)

        tspan = (0.0, 86400.0)
        EOM!(du, u, _p, t) = Cowell_EOM!(du, u, _p, t, model_list)

        event_log = []
        affect_fn! =
            (integrator) -> begin
                r = sqrt(integrator.u[1]^2 + integrator.u[2]^2 + integrator.u[3]^2)
                push!(event_log, (t=integrator.t, alt=r - 6378.137))
            end

        cb = ContinuousCallback(cond, affect_fn!, affect_fn!)

        prob = ODEProblem(EOM!, u0_cart, tspan, p)
        sol = solve(prob, Vern9(); abstol=1e-13, reltol=1e-13, callback=cb)

        @test length(event_log) >= 2
        for ev in event_log
            @test abs(ev.alt - h_target) < 1e-6
        end
    end

    # ================================================================
    # Geometric Detectors — Altitude (Non-Earth body)
    # ================================================================
    @testset "Altitude Detector (custom R_body)" begin
        R_moon = 1737.4
        cond = altitude_condition(Cartesian, 100.0; R_body=R_moon)
        mock_integrator = (p=p,)
        g = cond(u0_cart, 0.0, mock_integrator)
        r0 = sqrt(u0_cart[1]^2 + u0_cart[2]^2 + u0_cart[3]^2)
        @test g ≈ r0 - (R_moon + 100.0) rtol = 1e-12
    end

    # ================================================================
    # Geometric Detectors — Latitude
    # ================================================================
    @testset "Latitude Detector (Cowell)" begin
        lat_target = 0.0
        cond = latitude_condition(Cartesian, lat_target)

        tspan = (0.0, 86400.0)
        EOM!(du, u, _p, t) = Cowell_EOM!(du, u, _p, t, model_list)

        event_log = []
        affect_fn! =
            (integrator) -> begin
                r = sqrt(integrator.u[1]^2 + integrator.u[2]^2 + integrator.u[3]^2)
                lat = asin(integrator.u[3] / r)
                push!(event_log, (t=integrator.t, lat=lat))
            end

        cb = ContinuousCallback(cond, affect_fn!, affect_fn!)

        prob = ODEProblem(EOM!, u0_cart, tspan, p)
        sol = solve(prob, Vern9(); abstol=1e-13, reltol=1e-13, callback=cb)

        @test length(event_log) >= 2 * floor(Int, n_orbits_in_day)
        for ev in event_log
            @test abs(ev.lat) < 1e-9
        end
    end

    # ================================================================
    # Utility Detectors — Date
    # ================================================================
    @testset "Date Detector (Cowell)" begin
        t_target = 43200.0
        cond = date_condition(Cartesian, t_target)

        tspan = (0.0, 86400.0)
        EOM!(du, u, _p, t) = Cowell_EOM!(du, u, _p, t, model_list)

        event_time = Ref(0.0)
        cb = ContinuousCallback(cond, (integrator) -> begin
            event_time[] = integrator.t
        end)

        prob = ODEProblem(EOM!, u0_cart, tspan, p)
        sol = solve(prob, Vern9(); abstol=1e-13, reltol=1e-13, callback=cb)

        @test abs(event_time[] - t_target) < 1e-9
    end

    # ================================================================
    # Utility Detectors — Negate
    # ================================================================
    @testset "Negate Detector" begin
        cond = node_condition(Cartesian)
        neg_cond = negate_condition(cond)

        mock_integrator = (p=p,)
        g_val = cond(u0_cart, 0.0, mock_integrator)
        neg_g_val = neg_cond(u0_cart, 0.0, mock_integrator)
        @test g_val ≈ -neg_g_val
        @test g_val != 0.0
    end

    # ================================================================
    # Utility Detectors — Boolean Composition
    # ================================================================
    @testset "Boolean Detectors" begin
        g1 = node_condition(Cartesian)
        g2 = altitude_condition(Cartesian, 400.0)

        and_cond = and_condition(g1, g2)
        or_cond = or_condition(g1, g2)

        mock_integrator = (p=p,)
        v1 = g1(u0_cart, 0.0, mock_integrator)
        v2 = g2(u0_cart, 0.0, mock_integrator)

        @test and_cond(u0_cart, 0.0, mock_integrator) ≈ min(v1, v2)
        @test or_cond(u0_cart, 0.0, mock_integrator) ≈ max(v1, v2)

        same_g = and_condition(g1, g1)
        @test same_g(u0_cart, 0.0, mock_integrator) ≈ v1
    end

    # ================================================================
    # Utility Detectors — Shift
    # ================================================================
    @testset "Shift Condition" begin
        cond = date_condition(Cartesian, 1000.0)
        shifted = shift_condition(cond, 100.0, Cartesian)

        mock_integrator = (p=p,)
        @test cond(u0_cart, 1000.0, mock_integrator) ≈ 0.0
        @test shifted(u0_cart, 900.0, mock_integrator) ≈ 0.0
    end

    # ================================================================
    # Maneuver System — build_maneuver_callback (Cowell, time trigger)
    # ================================================================
    @testset "Build Maneuver Callback (Cowell, TimeTrigger)" begin
        deltaV = [0.05, 0.01, 0.01]

        cb = build_maneuver_callback(TimeTrigger(43200.0), FixedDeltaV(deltaV), Cartesian)

        tspan = (0.0, 86400.0)
        EOM!(du, u, _p, t) = Cowell_EOM!(du, u, _p, t, model_list)
        prob = ODEProblem(EOM!, u0_cart, tspan, p)
        sol = solve(prob, Vern9(); abstol=1e-13, reltol=1e-13, callback=cb)

        NRG_before = orbitalNRG(Cartesian(sol(43200.0 - 1.0)), μ)
        NRG_after = orbitalNRG(Cartesian(sol(43200.0 + 1.0)), μ)
        @test abs(NRG_after - NRG_before) > 1e-4

        expected_end = [
            29390.280395821836
            18637.945967159154
            -1768.361355756133
            0.47323343997331674
            2.572684107343496
            0.13273831002165992
        ]
        @test sol.u[end] ≈ expected_end rtol = 1e-8
    end

    # ================================================================
    # Maneuver System — build_maneuver_callback with RTN frame
    # ================================================================
    @testset "Build Maneuver Callback (Cowell, RTN Frame)" begin
        deltaV_rtn = [0.0, 0.05, 0.0]

        cb = build_maneuver_callback(
            TimeTrigger(43200.0), FixedDeltaV(deltaV_rtn, RTNFrame()), Cartesian
        )

        tspan = (0.0, 86400.0)
        EOM!(du, u, _p, t) = Cowell_EOM!(du, u, _p, t, model_list)
        prob = ODEProblem(EOM!, u0_cart, tspan, p)
        sol = solve(prob, Vern9(); abstol=1e-13, reltol=1e-13, callback=cb)

        NRG_before = orbitalNRG(Cartesian(sol(43200.0 - 1.0)), μ)
        NRG_after = orbitalNRG(Cartesian(sol(43200.0 + 1.0)), μ)
        @test abs(NRG_after - NRG_before) > 1e-4

        sol_noburn = solve(prob, Vern9(); abstol=1e-13, reltol=1e-13)
        @test !(sol.u[end] ≈ sol_noburn.u[end])
    end

    # ================================================================
    # Maneuver System — Event Triggered Burn
    # ================================================================
    @testset "Event Triggered Burn (Cowell, apoapsis)" begin
        deltaV = [0.01, 0.0, 0.0]

        apoapsis_cond = apside_condition(Cartesian)

        cb = build_maneuver_callback(
            EventTrigger(apoapsis_cond), FixedDeltaV(deltaV), Cartesian
        )

        tspan = (0.0, 86400.0)
        EOM!(du, u, _p, t) = Cowell_EOM!(du, u, _p, t, model_list)
        prob = ODEProblem(EOM!, u0_cart, tspan, p)
        sol = solve(prob, Vern9(); abstol=1e-13, reltol=1e-13, callback=cb)

        NRG = orbitalNRG.(Cartesian.(sol.u), μ)
        NRG_init = NRG[1]
        changes = findall(i -> abs(NRG[i] - NRG_init) > 1e-6, 2:length(NRG))
        @test length(changes) >= 1

        sol_noburn = solve(prob, Vern9(); abstol=1e-13, reltol=1e-13)
        @test !(sol.u[end] ≈ sol_noburn.u[end])
    end

    # ================================================================
    # Maneuver System — Computed DeltaV
    # ================================================================
    @testset "Computed DeltaV (Cowell)" begin
        compute_dv = (cart, t, p) -> SVector{3}(0.01, 0.0, 0.0)

        cb_computed = build_maneuver_callback(
            TimeTrigger(43200.0), ComputedDeltaV(compute_dv), Cartesian
        )
        cb_fixed = build_maneuver_callback(
            TimeTrigger(43200.0), FixedDeltaV([0.01, 0.0, 0.0]), Cartesian
        )

        tspan = (0.0, 86400.0)
        EOM!(du, u, _p, t) = Cowell_EOM!(du, u, _p, t, model_list)
        prob = ODEProblem(EOM!, u0_cart, tspan, p)

        sol_computed = solve(
            prob, Vern9(); abstol=1e-13, reltol=1e-13, callback=cb_computed
        )
        sol_fixed = solve(prob, Vern9(); abstol=1e-13, reltol=1e-13, callback=cb_fixed)

        @test sol_computed.u[end] ≈ sol_fixed.u[end] rtol = 1e-10
    end

    # ================================================================
    # Maneuver System — Schedule (verify BOTH burns execute)
    # ================================================================
    @testset "Maneuver Schedule (Cowell)" begin
        burns = [
            (TimeTrigger(21600.0), FixedDeltaV([0.01, 0.0, 0.0])),
            (TimeTrigger(43200.0), FixedDeltaV([0.0, 0.01, 0.0])),
        ]
        cbs = build_maneuver_schedule(burns, Cartesian)

        tspan = (0.0, 86400.0)
        EOM!(du, u, _p, t) = Cowell_EOM!(du, u, _p, t, model_list)
        prob = ODEProblem(EOM!, u0_cart, tspan, p)
        sol = solve(prob, Vern9(); abstol=1e-13, reltol=1e-13, callback=cbs)

        NRG_pre1 = orbitalNRG(Cartesian(sol(21600.0 - 1.0)), μ)
        NRG_post1 = orbitalNRG(Cartesian(sol(21600.0 + 1.0)), μ)
        NRG_pre2 = orbitalNRG(Cartesian(sol(43200.0 - 1.0)), μ)
        NRG_post2 = orbitalNRG(Cartesian(sol(43200.0 + 1.0)), μ)

        @test abs(NRG_post1 - NRG_pre1) > 1e-6
        @test abs(NRG_post2 - NRG_pre2) > 1e-6

        @test NRG_post1 ≈ NRG_pre2 rtol = 1e-6

        cb_single = build_maneuver_callback(
            TimeTrigger(21600.0), FixedDeltaV([0.01, 0.0, 0.0]), Cartesian
        )
        sol_single = solve(prob, Vern9(); abstol=1e-13, reltol=1e-13, callback=cb_single)
        @test !(sol.u[end] ≈ sol_single.u[end])
    end

    # ================================================================
    # Event Actions — Terminate
    # ================================================================
    @testset "TerminateAction" begin
        cond = apside_condition(Cartesian)
        cb = build_event_callback(cond, TerminateAction(), Cartesian)

        tspan = (0.0, 86400.0)
        EOM!(du, u, _p, t) = Cowell_EOM!(du, u, _p, t, model_list)
        prob = ODEProblem(EOM!, u0_cart, tspan, p)
        sol = solve(prob, Vern9(); abstol=1e-13, reltol=1e-13, callback=cb)

        @test sol.retcode == ReturnCode.Terminated
        @test sol.t[end] < T_orbit

        cart_end = Cartesian(sol.u[end])
        r = SVector{3}(cart_end[1], cart_end[2], cart_end[3])
        v = SVector{3}(cart_end[4], cart_end[5], cart_end[6])
        @test abs(dot(r, v)) < 1e-6
    end

    # ================================================================
    # Event Actions — Log
    # ================================================================
    @testset "LogAction" begin
        cond = apside_condition(Cartesian)
        log = []
        cb = build_event_callback(cond, LogAction(log), Cartesian)

        tspan = (0.0, 86400.0)
        EOM!(du, u, _p, t) = Cowell_EOM!(du, u, _p, t, model_list)
        prob = ODEProblem(EOM!, u0_cart, tspan, p)
        sol = solve(prob, Vern9(); abstol=1e-13, reltol=1e-13, callback=cb)

        @test length(log) >= 2 * floor(Int, n_orbits_in_day)
        @test all(e -> haskey(e, :t) && haskey(e, :state), log)

        for ev in log
            cart = ev.state
            r = SVector{3}(cart[1], cart[2], cart[3])
            v = SVector{3}(cart[4], cart[5], cart[6])
            @test abs(dot(r, v)) < 1e-6
        end
    end

    # ================================================================
    # Event Actions — ManeuverAction
    # ================================================================
    @testset "ManeuverAction" begin
        cond = apside_condition(Cartesian)
        action = ManeuverAction(FixedDeltaV([0.01, 0.0, 0.0]))
        cb = build_event_callback(cond, action, Cartesian)

        tspan = (0.0, 86400.0)
        EOM!(du, u, _p, t) = Cowell_EOM!(du, u, _p, t, model_list)
        prob = ODEProblem(EOM!, u0_cart, tspan, p)
        sol = solve(prob, Vern9(); abstol=1e-13, reltol=1e-13, callback=cb)

        NRG = orbitalNRG.(Cartesian.(sol.u), μ)
        @test NRG[1] != NRG[end]

        sol_noburn = solve(prob, Vern9(); abstol=1e-13, reltol=1e-13)
        @test !(sol.u[end] ≈ sol_noburn.u[end])
    end

    # ================================================================
    # Event Actions — ContinueAction
    # ================================================================
    @testset "ContinueAction" begin
        cond = apside_condition(Cartesian)
        cb = build_event_callback(cond, ContinueAction(), Cartesian)

        tspan = (0.0, 86400.0)
        EOM!(du, u, _p, t) = Cowell_EOM!(du, u, _p, t, model_list)
        prob = ODEProblem(EOM!, u0_cart, tspan, p)

        sol = solve(prob, Vern9(); abstol=1e-13, reltol=1e-13, callback=cb)
        sol_noburn = solve(prob, Vern9(); abstol=1e-13, reltol=1e-13)

        @test sol.u[end] ≈ sol_noburn.u[end] rtol = 1e-10
    end

    # ================================================================
    # Termination Callback
    # ================================================================
    @testset "build_termination_callback (Cowell)" begin
        cb = build_termination_callback(50000.0, Cartesian)

        tspan = (0.0, 86400.0)
        EOM!(du, u, _p, t) = Cowell_EOM!(du, u, _p, t, model_list)
        prob = ODEProblem(EOM!, u0_cart, tspan, p)
        sol = solve(prob, Vern9(); abstol=1e-13, reltol=1e-13, callback=cb)

        @test sol.retcode == ReturnCode.Terminated
        @test abs(sol.t[end] - 50000.0) < 1e-8
    end

    # ================================================================
    # Maneuver System — GaussVE propagator
    # ================================================================
    @testset "Build Maneuver Callback (GaussVE, TimeTrigger)" begin
        deltaV = [0.05, 0.01, 0.01]

        cb = build_maneuver_callback(TimeTrigger(43200.0), FixedDeltaV(deltaV), Keplerian)

        u0_koe = Array(Keplerian(Cartesian(u0_cart), μ))
        tspan = (0.0, 86400.0)
        EOM!(du, u, _p, t) = GaussVE_EOM!(du, u, _p, t, model_list)
        prob = ODEProblem(EOM!, u0_koe, tspan, p)
        sol = solve(prob, Vern9(); abstol=1e-13, reltol=1e-13, callback=cb)

        expected_end = [
            29390.280395821836
            18637.945967159154
            -1768.361355756133
            0.47323343997331674
            2.572684107343496
            0.13273831002165992
        ]
        @test Cartesian(Keplerian(sol.u[end]), μ) ≈ expected_end rtol = 1e-4
    end

    # ================================================================
    # Maneuver System — EDromo propagator
    # ================================================================
    @testset "Build Maneuver Callback (EDromo, TimeTrigger)" begin
        grav_model = setup.grav_model

        W = (
            potential(Cartesian(u0_cart), p, 0.0, grav_model) -
            potential(Cartesian(u0_cart), p, 0.0, KeplerianGravityAstroModel(μ=μ))
        )
        config = RegularizedCoordinateConfig(
            u0_cart, μ; W=W, t₀=0.0, flag_time=PhysicalTime()
        )

        ϕ₀ = compute_initial_phi(u0_cart, μ, config)
        u0_edromo = Array(EDromo(Cartesian(u0_cart), μ, ϕ₀, config))
        tspan = (ϕ₀, ϕ₀ + 200π)

        deltaV = [0.05, 0.01, 0.01]

        burn_cb = build_maneuver_callback(
            TimeTrigger(43200.0), FixedDeltaV(deltaV), EDromo, config
        )
        end_cb = build_termination_callback(86400.0, EDromo, config)
        cbs = CallbackSet(burn_cb, end_cb)

        EOM!(du, u, _p, t) = EDromo_EOM!(du, u, _p, t, model_list, config)
        prob = ODEProblem(EOM!, u0_edromo, tspan, p)
        sol = solve(prob, Vern9(); abstol=1e-13, reltol=1e-13, callback=cbs)

        final_state = Cartesian(EDromo(sol.u[end]), μ, sol.t[end], config)

        expected_end = [
            29390.280395821836
            18637.945967159154
            -1768.361355756133
            0.47323343997331674
            2.572684107343496
            0.13273831002165992
        ]
        @test final_state ≈ expected_end rtol = 1e-4
    end

    # ================================================================
    # Maneuver System — KS propagator
    # ================================================================
    @testset "Build Maneuver Callback (KS, TimeTrigger)" begin
        grav_model = setup.grav_model

        W = (
            potential(Cartesian(u0_cart), p, 0.0, grav_model) -
            potential(Cartesian(u0_cart), p, 0.0, KeplerianGravityAstroModel(μ=μ))
        )
        config = RegularizedCoordinateConfig(
            u0_cart, μ; W=W, t₀=0.0, flag_time=PhysicalTime()
        )

        u0_ks = Array(KustaanheimoStiefel(Cartesian(u0_cart), μ, config))
        tspan = (0.0, 9π)

        deltaV = [0.05, 0.01, 0.01]

        burn_cb = build_maneuver_callback(
            TimeTrigger(43200.0), FixedDeltaV(deltaV), KustaanheimoStiefel, config
        )
        end_cb = build_termination_callback(86400.0, KustaanheimoStiefel, config)
        cbs = CallbackSet(burn_cb, end_cb)

        EOM!(du, u, _p, t) = KS_EOM!(du, u, _p, t, model_list, config)
        prob = ODEProblem(EOM!, u0_ks, tspan, p)
        sol = solve(prob, Vern9(); abstol=1e-13, reltol=1e-13, callback=cbs)

        final_state = Cartesian(KustaanheimoStiefel(sol.u[end]), μ, config)

        expected_end = [
            29390.280395821836
            18637.945967159154
            -1768.361355756133
            0.47323343997331674
            2.572684107343496
            0.13273831002165992
        ]
        @test final_state ≈ expected_end rtol = 1e-4
    end

    # ================================================================
    # Orbital Detectors — Regularized propagator (EDromo)
    # ================================================================
    @testset "Apside Detector (EDromo)" begin
        grav_model = setup.grav_model

        W = (
            potential(Cartesian(u0_cart), p, 0.0, grav_model) -
            potential(Cartesian(u0_cart), p, 0.0, KeplerianGravityAstroModel(μ=μ))
        )
        config = RegularizedCoordinateConfig(
            u0_cart, μ; W=W, t₀=0.0, flag_time=PhysicalTime()
        )

        ϕ₀ = compute_initial_phi(u0_cart, μ, config)
        u0_edromo = Array(EDromo(Cartesian(u0_cart), μ, ϕ₀, config))
        tspan = (ϕ₀, ϕ₀ + 200π)

        cond = apside_condition(EDromo, config)
        end_cb = build_termination_callback(86400.0, EDromo, config)

        peri_log = []
        apo_log = []
        peri_fn! =
            (integrator) -> begin
                cart = Cartesian(EDromo(integrator.u), μ, integrator.t, config)
                push!(peri_log, sqrt(cart[1]^2 + cart[2]^2 + cart[3]^2))
            end
        apo_fn! =
            (integrator) -> begin
                cart = Cartesian(EDromo(integrator.u), μ, integrator.t, config)
                push!(apo_log, sqrt(cart[1]^2 + cart[2]^2 + cart[3]^2))
            end

        apside_cb = ContinuousCallback(cond, peri_fn!, apo_fn!)
        cbs = CallbackSet(apside_cb, end_cb)

        EOM!(du, u, _p, t) = EDromo_EOM!(du, u, _p, t, model_list, config)
        prob = ODEProblem(EOM!, u0_edromo, tspan, p)
        sol = solve(prob, Vern9(); abstol=1e-13, reltol=1e-13, callback=cbs)

        n_expected = floor(Int, n_orbits_in_day)
        @test length(peri_log) >= n_expected
        @test length(apo_log) >= n_expected

        for r in peri_log
            @test r ≈ r_peri rtol = 1e-4
        end
        for r in apo_log
            @test r ≈ r_apo rtol = 1e-4
        end
    end

    # ================================================================
    # Date Detector — EDromo (regularized)
    # ================================================================
    @testset "Date Detector (EDromo)" begin
        grav_model = setup.grav_model

        W = (
            potential(Cartesian(u0_cart), p, 0.0, grav_model) -
            potential(Cartesian(u0_cart), p, 0.0, KeplerianGravityAstroModel(μ=μ))
        )
        config = RegularizedCoordinateConfig(
            u0_cart, μ; W=W, t₀=0.0, flag_time=PhysicalTime()
        )

        ϕ₀ = compute_initial_phi(u0_cart, μ, config)
        u0_edromo = Array(EDromo(Cartesian(u0_cart), μ, ϕ₀, config))
        tspan = (ϕ₀, ϕ₀ + 200π)

        t_target = 43200.0
        cond = date_condition(EDromo, config, t_target)
        end_cb = build_termination_callback(86400.0, EDromo, config)

        event_phys_time = Ref(0.0)
        detect_cb = ContinuousCallback(
            cond,
            (integrator) -> begin
                event_phys_time[] = get_physical_time(
                    integrator.u, integrator.t, EDromo, config
                )
            end,
        )

        EOM!(du, u, _p, t) = EDromo_EOM!(du, u, _p, t, model_list, config)
        prob = ODEProblem(EOM!, u0_edromo, tspan, p)
        sol = solve(
            prob,
            Vern9();
            abstol=1e-13,
            reltol=1e-13,
            callback=CallbackSet(detect_cb, end_cb),
        )

        @test abs(event_phys_time[] - t_target) < 1e-6
    end
end
