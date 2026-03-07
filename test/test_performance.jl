@testset "Aqua.jl" begin
    Aqua.test_all(
        AstroPropagators;
        ambiguities=(recursive = false),
        deps_compat=(check_extras = false),
    )
end

@testset "JET Testing" begin
    rep = JET.test_package(
        AstroPropagators; toplevel_logger=nothing, target_modules=(@__MODULE__,)
    )
end

@testset "EOM Allocations" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)

    SpaceIndices.init()
    eop_data = fetch_iers_eop()
    grav_coeffs = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96))

    grav_model = GravityHarmonicsAstroModel(;
        gravity_model=grav_coeffs,
        eop_data=eop_data,
        order=36,
        degree=36,
        P=zeros(37, 37),
        dP=zeros(37, 37),
    )
    p = ComponentVector(;
        JD=JD, μ=GravityModels.gravity_constant(grav_model.gravity_model) / 1E9
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

    model_list = CentralBodyDynamicsModel(
        grav_model, (sun_third_body, moon_third_body, srp_model, drag_model)
    )

    @test length(
        check_allocs(GaussVE_EOM, (Vector{Float64}, typeof(p), Float64, typeof(model_list)))
    ) == 0
    @test length(
        check_allocs(
            GaussVE_EOM!,
            (Vector{Float64}, Vector{Float64}, typeof(p), Float64, typeof(model_list)),
        ),
    ) == 0
    @test length(
        check_allocs(Cowell_EOM, (Vector{Float64}, typeof(p), Float64, typeof(model_list)))
    ) == 0
    @test length(
        check_allocs(
            Cowell_EOM!,
            (Vector{Float64}, Vector{Float64}, typeof(p), Float64, typeof(model_list)),
        ),
    ) == 0
    @test length(
        check_allocs(USM7_EOM, (Vector{Float64}, typeof(p), Float64, typeof(model_list)))
    ) == 0
    @test length(
        check_allocs(
            USM7_EOM!,
            (Vector{Float64}, Vector{Float64}, typeof(p), Float64, typeof(model_list)),
        ),
    ) == 0
    @test length(
        check_allocs(USM6_EOM, (Vector{Float64}, typeof(p), Float64, typeof(model_list)))
    ) == 0
    @test length(
        check_allocs(
            USM6_EOM!,
            (Vector{Float64}, Vector{Float64}, typeof(p), Float64, typeof(model_list)),
        ),
    ) == 0
    @test length(
        check_allocs(USMEM_EOM, (Vector{Float64}, typeof(p), Float64, typeof(model_list)))
    ) == 0
    @test length(
        check_allocs(
            USMEM_EOM!,
            (Vector{Float64}, Vector{Float64}, typeof(p), Float64, typeof(model_list)),
        ),
    ) == 0
    @test length(
        check_allocs(
            Milankovich_EOM, (Vector{Float64}, typeof(p), Float64, typeof(model_list))
        ),
    ) == 0
    @test length(
        check_allocs(
            Milankovich_EOM!,
            (Vector{Float64}, Vector{Float64}, typeof(p), Float64, typeof(model_list)),
        ),
    ) == 0
end

@testset "Regularized Coordinate EOM Allocations" begin
    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)

    SpaceIndices.init()
    eop_data = fetch_iers_eop()
    grav_coeffs = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96))

    grav_model = GravityHarmonicsAstroModel(;
        gravity_model=grav_coeffs,
        eop_data=eop_data,
        order=36,
        degree=36,
        P=zeros(37, 37),
        dP=zeros(37, 37),
    )
    p = ComponentVector(;
        JD=JD, μ=GravityModels.gravity_constant(grav_model.gravity_model) / 1E9
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

    model_list = CentralBodyDynamicsModel(
        grav_model, (sun_third_body, moon_third_body, srp_model, drag_model)
    )

    u0_cart = [
        29390.280395821836
        18637.945967159154
        -1768.361355756133
        0.47323343997331674
        2.572684107343496
        0.13273831002165992
    ]

    μ = p.μ

    W = (
        potential(Cartesian(u0_cart), p, 0.0, grav_model) -
        potential(Cartesian(u0_cart), p, 0.0, KeplerianGravityAstroModel(μ=μ))
    )

    config = RegularizedCoordinateConfig(u0_cart, μ; W=W, t₀=0.0, flag_time=PhysicalTime())

    @test length(
        check_allocs(
            (u, p, ϕ, models) -> EDromo_EOM(u, p, ϕ, models, config),
            (Vector{Float64}, typeof(p), Float64, typeof(model_list)),
        ),
    ) == 0
    @test length(
        check_allocs(
            (du, u, p, ϕ, models) -> EDromo_EOM!(du, u, p, ϕ, models, config),
            (Vector{Float64}, Vector{Float64}, typeof(p), Float64, typeof(model_list)),
        ),
    ) == 0

    @test length(
        check_allocs(
            (u, p, ϕ, models) -> KS_EOM(u, p, ϕ, models, config),
            (Vector{Float64}, typeof(p), Float64, typeof(model_list)),
        ),
    ) == 0
    @test length(
        check_allocs(
            (du, u, p, ϕ, models) -> KS_EOM!(du, u, p, ϕ, models, config),
            (Vector{Float64}, Vector{Float64}, typeof(p), Float64, typeof(model_list)),
        ),
    ) == 0

    @test length(
        check_allocs(
            (u, p, ϕ, models) -> StiSche_EOM(u, p, ϕ, models, config),
            (Vector{Float64}, typeof(p), Float64, typeof(model_list)),
        ),
    ) == 0
    @test length(
        check_allocs(
            (du, u, p, ϕ, models) -> StiSche_EOM!(du, u, p, ϕ, models, config),
            (Vector{Float64}, Vector{Float64}, typeof(p), Float64, typeof(model_list)),
        ),
    ) == 0

    config_geqoe = RegularizedCoordinateConfig(; W=W)

    @test length(
        check_allocs(
            (u, p, t, models) -> GEqOE_EOM(u, p, t, models, config_geqoe),
            (Vector{Float64}, typeof(p), Float64, typeof(model_list)),
        ),
    ) == 0
    @test length(
        check_allocs(
            (du, u, p, t, models) -> GEqOE_EOM!(du, u, p, t, models, config_geqoe),
            (Vector{Float64}, Vector{Float64}, typeof(p), Float64, typeof(model_list)),
        ),
    ) == 0
end

@testset "Event Condition Allocations" begin
    u0_cart = [
        -1076.225324679696
        -6765.896364327722
        -332.3087833503755
        9.356857417032581
        -3.3123476319597557
        -1.1880157328553503
    ]

    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
    grav_model = KeplerianGravityAstroModel()
    μ = grav_model.μ
    p = ComponentVector(; JD=JD, μ=μ)

    mock_integrator = (p=p,)
    MI = typeof(mock_integrator)

    # -- State adapter: get_cartesian / get_keplerian --
    @testset "get_cartesian (Cartesian)" begin
        @test length(
            check_allocs(
                get_cartesian, (Vector{Float64}, Float64, Float64, Type{Cartesian})
            ),
        ) == 0
    end

    @testset "get_cartesian (Keplerian)" begin
        @test length(
            check_allocs(
                get_cartesian, (Vector{Float64}, Float64, Float64, Type{Keplerian})
            ),
        ) == 0
    end

    @testset "get_keplerian (Cartesian)" begin
        @test length(
            check_allocs(
                get_keplerian, (Vector{Float64}, Float64, Float64, Type{Cartesian})
            ),
        ) == 0
    end

    # -- Orbital detector conditions --
    @testset "apside_condition (Cartesian)" begin
        cond = apside_condition(Cartesian)
        @test length(check_allocs(cond, (Vector{Float64}, Float64, MI))) == 0
    end

    @testset "node_condition (Cartesian)" begin
        cond = node_condition(Cartesian)
        @test length(check_allocs(cond, (Vector{Float64}, Float64, MI))) == 0
    end

    @testset "true_anomaly_condition (Cartesian)" begin
        cond = true_anomaly_condition(Cartesian, π)
        @test length(check_allocs(cond, (Vector{Float64}, Float64, MI))) == 0
    end

    @testset "argument_of_latitude_condition (Cartesian)" begin
        cond = argument_of_latitude_condition(Cartesian, π / 2)
        @test length(check_allocs(cond, (Vector{Float64}, Float64, MI))) == 0
    end

    @testset "mean_anomaly_condition (Cartesian)" begin
        cond = mean_anomaly_condition(Cartesian, 0.0)
        @test length(check_allocs(cond, (Vector{Float64}, Float64, MI))) == 0
    end

    @testset "raan_condition (Cartesian)" begin
        cond = raan_condition(Cartesian, 0.0)
        @test length(check_allocs(cond, (Vector{Float64}, Float64, MI))) == 0
    end

    # -- Geometric detector conditions --
    @testset "altitude_condition (Cartesian)" begin
        cond = altitude_condition(Cartesian, 400.0)
        @test length(check_allocs(cond, (Vector{Float64}, Float64, MI))) == 0
    end

    @testset "latitude_condition (Cartesian)" begin
        cond = latitude_condition(Cartesian, 0.0)
        @test length(check_allocs(cond, (Vector{Float64}, Float64, MI))) == 0
    end

    @testset "longitude_condition (Cartesian)" begin
        cond = longitude_condition(Cartesian, 0.0)
        @test length(check_allocs(cond, (Vector{Float64}, Float64, MI))) == 0
    end

    # -- Utility detector conditions --
    @testset "date_condition (Cartesian)" begin
        cond = date_condition(Cartesian, 43200.0)
        @test length(check_allocs(cond, (Vector{Float64}, Float64, MI))) == 0
    end

    @testset "negate_condition" begin
        cond = negate_condition(node_condition(Cartesian))
        @test length(check_allocs(cond, (Vector{Float64}, Float64, MI))) == 0
    end

    @testset "and_condition" begin
        cond = and_condition(
            node_condition(Cartesian), altitude_condition(Cartesian, 400.0)
        )
        @test length(check_allocs(cond, (Vector{Float64}, Float64, MI))) == 0
    end

    @testset "or_condition" begin
        cond = or_condition(node_condition(Cartesian), altitude_condition(Cartesian, 400.0))
        @test length(check_allocs(cond, (Vector{Float64}, Float64, MI))) == 0
    end

    @testset "shift_condition" begin
        cond = shift_condition(node_condition(Cartesian), 60.0, Cartesian)
        @test length(check_allocs(cond, (Vector{Float64}, Float64, MI))) == 0
    end

    # -- Eclipse and beta_angle conditions (require sun data) --
    eop_data = fetch_iers_eop()
    sun_third_body = ThirdBodyModel(; body=SunBody(), eop_data=eop_data)

    @testset "eclipse_condition (Cartesian)" begin
        cond = eclipse_condition(Cartesian, Conical(), sun_third_body)
        @test length(check_allocs(cond, (Vector{Float64}, Float64, MI))) == 0
    end

    @testset "beta_angle_condition (Cartesian)" begin
        cond = beta_angle_condition(Cartesian, sun_third_body, deg2rad(60.0))
        @test length(check_allocs(cond, (Vector{Float64}, Float64, MI))) == 0
    end

    # -- Keplerian-based conditions (with coordinate conversion) --
    @testset "apside_condition (Keplerian)" begin
        cond = apside_condition(Keplerian)
        @test length(check_allocs(cond, (Vector{Float64}, Float64, MI))) == 0
    end

    @testset "node_condition (Keplerian)" begin
        cond = node_condition(Keplerian)
        @test length(check_allocs(cond, (Vector{Float64}, Float64, MI))) == 0
    end

    # -- Regularized coordinate conditions --
    @testset "Regularized condition allocations" begin
        W = (
            potential(Cartesian(u0_cart), p, 0.0, grav_model) -
            potential(Cartesian(u0_cart), p, 0.0, KeplerianGravityAstroModel(μ=μ))
        )
        config = RegularizedCoordinateConfig(
            u0_cart, μ; W=W, t₀=0.0, flag_time=PhysicalTime()
        )

        @testset "apside_condition (EDromo)" begin
            cond = apside_condition(EDromo, config)
            @test length(check_allocs(cond, (Vector{Float64}, Float64, MI))) == 0
        end

        @testset "node_condition (EDromo)" begin
            cond = node_condition(EDromo, config)
            @test length(check_allocs(cond, (Vector{Float64}, Float64, MI))) == 0
        end

        @testset "date_condition (EDromo)" begin
            cond = date_condition(EDromo, config, 43200.0)
            @test length(check_allocs(cond, (Vector{Float64}, Float64, MI))) == 0
        end

        @testset "get_physical_time (EDromo)" begin
            @test length(
                check_allocs(
                    (u, t) -> get_physical_time(u, t, EDromo, config),
                    (Vector{Float64}, Float64),
                ),
            ) == 0
        end

        @testset "apside_condition (KustaanheimoStiefel)" begin
            cond = apside_condition(KustaanheimoStiefel, config)
            @test length(check_allocs(cond, (Vector{Float64}, Float64, MI))) == 0
        end

        @testset "apside_condition (StiefelScheifele)" begin
            cond = apside_condition(StiefelScheifele, config)
            @test length(check_allocs(cond, (Vector{Float64}, Float64, MI))) == 0
        end

        config_geqoe = RegularizedCoordinateConfig(; W=W)

        @testset "apside_condition (GEqOE)" begin
            cond = apside_condition(GEqOE, config_geqoe)
            @test length(check_allocs(cond, (Vector{Float64}, Float64, MI))) == 0
        end
    end

    # -- Maneuver builder conditions --
    @testset "TimeTrigger condition (Cartesian)" begin
        trigger = TimeTrigger(43200.0)
        cond = AstroPropagators._build_condition(trigger, Cartesian)
        @test length(check_allocs(cond, (Vector{Float64}, Float64, MI))) == 0
    end

    @testset "TimeTrigger condition (EDromo)" begin
        W = (
            potential(Cartesian(u0_cart), p, 0.0, grav_model) -
            potential(Cartesian(u0_cart), p, 0.0, KeplerianGravityAstroModel(μ=μ))
        )
        config = RegularizedCoordinateConfig(
            u0_cart, μ; W=W, t₀=0.0, flag_time=PhysicalTime()
        )
        trigger = TimeTrigger(43200.0)
        cond = AstroPropagators._build_condition(trigger, EDromo, config)
        @test length(check_allocs(cond, (Vector{Float64}, Float64, MI))) == 0
    end
end

@testset "Maneuver Affect Allocations" begin
    u0_cart = SVector{6,Float64}(
        -1076.225324679696,
        -6765.896364327722,
        -332.3087833503755,
        9.356857417032581,
        -3.3123476319597557,
        -1.1880157328553503,
    )

    JD = date_to_jd(2024, 1, 5, 12, 0, 0.0)
    grav_model = KeplerianGravityAstroModel()
    μ = grav_model.μ
    p = ComponentVector(; JD=JD, μ=μ)

    mock_cart = _MockIntegrator(u0_cart, 0.0, p)
    MIC = typeof(mock_cart)

    # -- _apply_burn! --
    @testset "_apply_burn! (Cartesian)" begin
        @test length(
            check_allocs(
                AstroPropagators._apply_burn!, (MIC, SVector{3,Float64}, Type{Cartesian})
            ),
        ) == 0
    end

    @testset "_apply_burn! (Keplerian)" begin
        u0_kep = params(Keplerian(Cartesian(u0_cart), μ))
        mock_kep = _MockIntegrator(u0_kep, 0.0, p)
        MIK = typeof(mock_kep)
        @test length(
            check_allocs(
                AstroPropagators._apply_burn!, (MIK, SVector{3,Float64}, Type{Keplerian})
            ),
        ) == 0
    end

    # -- _build_affect: FixedDeltaV with all thrust frames --
    @testset "FixedDeltaV affect (Cartesian, InertialFrame)" begin
        effect = FixedDeltaV(SVector{3}(0.1, 0.0, 0.0), InertialFrame())
        affect! = AstroPropagators._build_affect(effect, Cartesian)
        @test length(check_allocs(affect!, (MIC,))) == 0
    end

    @testset "FixedDeltaV affect (Cartesian, RTNFrame)" begin
        effect = FixedDeltaV(SVector{3}(0.1, 0.0, 0.0), RTNFrame())
        affect! = AstroPropagators._build_affect(effect, Cartesian)
        @test length(check_allocs(affect!, (MIC,))) == 0
    end

    @testset "FixedDeltaV affect (Cartesian, VNBFrame)" begin
        effect = FixedDeltaV(SVector{3}(0.1, 0.0, 0.0), VNBFrame())
        affect! = AstroPropagators._build_affect(effect, Cartesian)
        @test length(check_allocs(affect!, (MIC,))) == 0
    end

    # -- _build_affect: ComputedDeltaV --
    @testset "ComputedDeltaV affect (Cartesian, InertialFrame)" begin
        compute = (cart, t, p) -> SVector{3,Float64}(0.1, 0.0, 0.0)
        effect = ComputedDeltaV(compute, InertialFrame())
        affect! = AstroPropagators._build_affect(effect, Cartesian)
        @test length(check_allocs(affect!, (MIC,))) == 0
    end

    @testset "ComputedDeltaV affect (Cartesian, RTNFrame)" begin
        compute = (cart, t, p) -> SVector{3,Float64}(0.0, 0.1, 0.0)
        effect = ComputedDeltaV(compute, RTNFrame())
        affect! = AstroPropagators._build_affect(effect, Cartesian)
        @test length(check_allocs(affect!, (MIC,))) == 0
    end

    # -- _build_affect: FixedDeltaV with Keplerian coord type --
    @testset "FixedDeltaV affect (Keplerian, InertialFrame)" begin
        u0_kep = params(Keplerian(Cartesian(u0_cart), μ))
        mock_kep = _MockIntegrator(u0_kep, 0.0, p)
        MIK = typeof(mock_kep)
        effect = FixedDeltaV(SVector{3}(0.1, 0.0, 0.0), InertialFrame())
        affect! = AstroPropagators._build_affect(effect, Keplerian)
        @test length(check_allocs(affect!, (MIK,))) == 0
    end

    # -- Regularized coordinate affects --
    @testset "Regularized affect allocations" begin
        W = (
            potential(Cartesian(u0_cart), p, 0.0, grav_model) -
            potential(Cartesian(u0_cart), p, 0.0, KeplerianGravityAstroModel(μ=μ))
        )
        config = RegularizedCoordinateConfig(
            Vector(u0_cart), μ; W=W, t₀=0.0, flag_time=PhysicalTime()
        )

        @testset "_apply_burn! (EDromo)" begin
            u0_ed = params(EDromo(Cartesian(u0_cart), μ, 0.0, config))
            mock_ed = _MockIntegrator(u0_ed, 0.0, p)
            MIE = typeof(mock_ed)
            @test length(
                check_allocs(
                    (integ, dv) -> AstroPropagators._apply_burn!(integ, dv, EDromo, config),
                    (MIE, SVector{3,Float64}),
                ),
            ) == 0
        end

        @testset "FixedDeltaV affect (EDromo)" begin
            u0_ed = params(EDromo(Cartesian(u0_cart), μ, 0.0, config))
            mock_ed = _MockIntegrator(u0_ed, 0.0, p)
            MIE = typeof(mock_ed)
            effect = FixedDeltaV(SVector{3}(0.1, 0.0, 0.0), InertialFrame())
            affect! = AstroPropagators._build_affect(effect, EDromo, config)
            @test length(check_allocs(affect!, (MIE,))) == 0
        end

        @testset "_apply_burn! (KustaanheimoStiefel)" begin
            u0_ks = params(KustaanheimoStiefel(Cartesian(u0_cart), μ, config))
            mock_ks = _MockIntegrator(u0_ks, 0.0, p)
            MIKS = typeof(mock_ks)
            @test length(
                check_allocs(
                    (integ, dv) -> AstroPropagators._apply_burn!(
                        integ, dv, KustaanheimoStiefel, config
                    ),
                    (MIKS, SVector{3,Float64}),
                ),
            ) == 0
        end

        @testset "_apply_burn! (StiefelScheifele)" begin
            ϕ₀ = compute_initial_phi(Vector(u0_cart), μ, config)
            u0_ss = params(StiefelScheifele(Cartesian(u0_cart), μ, ϕ₀, config))
            mock_ss = _MockIntegrator(u0_ss, ϕ₀, p)
            MISS = typeof(mock_ss)
            @test length(
                check_allocs(
                    (integ, dv) ->
                        AstroPropagators._apply_burn!(integ, dv, StiefelScheifele, config),
                    (MISS, SVector{3,Float64}),
                ),
            ) == 0
        end

        config_geqoe = RegularizedCoordinateConfig(; W=W)

        @testset "_apply_burn! (GEqOE)" begin
            u0_ge = params(GEqOE(Cartesian(u0_cart), μ, config_geqoe))
            mock_ge = _MockIntegrator(u0_ge, 0.0, p)
            MIGE = typeof(mock_ge)
            @test length(
                check_allocs(
                    (integ, dv) ->
                        AstroPropagators._apply_burn!(integ, dv, GEqOE, config_geqoe),
                    (MIGE, SVector{3,Float64}),
                ),
            ) == 0
        end
    end

    # -- Nth-occurrence wrappers --
    @testset "_wrap_condition_with_count" begin
        cond = apside_condition(Cartesian)
        wrapped = AstroPropagators._wrap_condition_with_count(cond, 3)
        MI_cond = typeof((p=p,))
        @test length(check_allocs(wrapped, (Vector{Float64}, Float64, MI_cond))) == 0
    end

    @testset "_wrap_affect_with_count" begin
        effect = FixedDeltaV(SVector{3}(0.1, 0.0, 0.0), InertialFrame())
        affect! = AstroPropagators._build_affect(effect, Cartesian)
        wrapped = AstroPropagators._wrap_affect_with_count(affect!, 3)
        @test length(check_allocs(wrapped, (MIC,))) == 0
    end
end
