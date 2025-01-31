@testset "Aqua.jl" begin
    Aqua.test_all(AstroPropagators; ambiguities=(recursive = false))
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
        gravity_model=grav_coeffs, eop_data=eop_data, order=36, degree=36, P=zeros(37, 37), dP=zeros(37, 37)
    )
    p = ComponentVector(;
        JD=JD, μ=GravityModels.gravity_constant(grav_model.gravity_model) / 1E9
    )

    sun_third_body = ThirdBodyModel(; body=SunBody(), eop_data=eop_data)
    moon_third_body = ThirdBodyModel(; body=MoonBody(), eop_data=eop_data)

    satellite_srp_model = CannonballFixedSRP(0.2)
    srp_model = SRPAstroModel(; satellite_srp_model=satellite_srp_model, sun_data=sun_third_body, eop_data=eop_data, shadow_model=Conical())

    satellite_drag_model = CannonballFixedDrag(0.2)
    drag_model = DragAstroModel(;satellite_drag_model=satellite_drag_model, atmosphere_model=JB2008(), eop_data=eop_data)

    model_list = (grav_model, sun_third_body, moon_third_body, srp_model, drag_model)

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
