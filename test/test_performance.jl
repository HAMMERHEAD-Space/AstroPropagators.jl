@testset "Aqua.jl" begin
    Aqua.test_all(
        AstroPropagators;
        ambiguities=(recursive = false),
        deps_compat=(check_extras = false),
    )
end

@testset "JET Testing" begin
    rep = JET.test_package(AstroPropagators; toplevel_logger=nothing)
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
end
