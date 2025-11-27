@testset "USM7 State Differentiability" begin
    SpaceIndices.init()
    for backend in _BACKENDS
        testname = "USM7 Differentiability " * backend[1]
        @testset "$testname" begin
            f_fd, df_fd = value_and_jacobian(
                (x) -> USM7_EOM(x, _p2, _t, _model_list),
                AutoFiniteDiff(),
                Array(_state_usm7),
            )

            f_ad, df_ad = value_and_jacobian(
                (x) -> Array(USM7_EOM(x, _p2, _t, _model_list)),
                backend[2],
                Array(_state_usm7),
            )

            @test f_fd ≈ f_ad
            @test df_fd ≈ df_ad atol = 1e-5
        end
    end
    SpaceIndices.destroy()
end

@testset "USM7 Time Differentiability" begin
    SpaceIndices.init()
    for backend in _BACKENDS
        testname = "USM7 Differentiability " * backend[1]
        @testset "$testname" begin
            f_fd, df_fd = value_and_derivative(
                (x) -> USM7_EOM(Array(_state_usm7), _p2, x, _model_list),
                AutoFiniteDiff(),
                _t,
            )

            f_ad, df_ad = value_and_derivative(
                (x) -> Array(USM7_EOM(Array(_state_usm7), _p2, x, _model_list)),
                backend[2],
                _t,
            )

            @test f_fd ≈ f_ad
            @test df_fd ≈ df_ad atol = 1e-5
        end
    end
    SpaceIndices.destroy()
end

@testset "USM7 Parameter Differentiability" begin
    SpaceIndices.init()

    function dynamics_params(x::AbstractArray{T}) where {T<:Number}
        satellite_srp_model = CannonballFixedSRP(x[1])
        srp_model = SRPAstroModel(;
            satellite_srp_model=satellite_srp_model,
            sun_data=_sun_model,
            eop_data=_eop_data,
            shadow_model=Conical(),
        )

        satellite_drag_model = CannonballFixedDrag(x[2])
        drag_model = DragAstroModel(;
            satellite_drag_model=satellite_drag_model,
            atmosphere_model=JB2008(),
            eop_data=_eop_data,
        )

        models = (_sun_model, _moon_model, srp_model, drag_model)
        model_list = CentralBodyDynamicsModel(_grav_model, models)

        return USM7_EOM(_state_usm7, _p2, _t, model_list)
    end

    for backend in _BACKENDS
        if backend[1] == "Enzyme"
            backend = (
                "Enzyme", AutoEnzyme(; mode=Enzyme.set_runtime_activity(Enzyme.Forward))
            )
        end
        testname = "USM7 Differentiability " * backend[1]
        @testset "$testname" begin
            f_fd, df_fd = value_and_jacobian(
                (x) -> Array(dynamics_params(x)), AutoFiniteDiff(), [0.2; 0.2]
            )

            f_ad, df_ad = value_and_jacobian(
                (x) -> Array(dynamics_params(x)), backend[2], [0.2; 0.2]
            )

            @test f_fd ≈ f_ad
            @test df_fd ≈ df_ad atol = 1e-5
        end
    end
    SpaceIndices.destroy()
end

@testset "USM6 State Differentiability" begin
    SpaceIndices.init()
    for backend in _BACKENDS
        testname = "USM6 Differentiability " * backend[1]
        @testset "$testname" begin
            f_fd, df_fd = value_and_jacobian(
                (x) -> USM6_EOM(x, _p2, _t, _model_list),
                AutoFiniteDiff(),
                Array(_state_usm6),
            )

            f_ad, df_ad = value_and_jacobian(
                (x) -> Array(USM6_EOM(x, _p2, _t, _model_list)),
                backend[2],
                Array(_state_usm6),
            )

            @test f_fd ≈ f_ad
            @test df_fd ≈ df_ad atol = 1e-5
        end
    end
    SpaceIndices.destroy()
end

@testset "USM6 Time Differentiability" begin
    SpaceIndices.init()
    for backend in _BACKENDS
        testname = "USM6 Differentiability " * backend[1]
        @testset "$testname" begin
            f_fd, df_fd = value_and_derivative(
                (x) -> USM6_EOM(Array(_state_usm6), _p2, x, _model_list),
                AutoFiniteDiff(),
                _t,
            )

            f_ad, df_ad = value_and_derivative(
                (x) -> Array(USM6_EOM(Array(_state_usm6), _p2, x, _model_list)),
                backend[2],
                _t,
            )

            @test f_fd ≈ f_ad
            @test df_fd ≈ df_ad atol = 1e-5
        end
    end
    SpaceIndices.destroy()
end

@testset "USM6 Parameter Differentiability" begin
    SpaceIndices.init()

    function dynamics_params(x::AbstractArray{T}) where {T<:Number}
        satellite_srp_model = CannonballFixedSRP(x[1])
        srp_model = SRPAstroModel(;
            satellite_srp_model=satellite_srp_model,
            sun_data=_sun_model,
            eop_data=_eop_data,
            shadow_model=Conical(),
        )

        satellite_drag_model = CannonballFixedDrag(x[2])
        drag_model = DragAstroModel(;
            satellite_drag_model=satellite_drag_model,
            atmosphere_model=JB2008(),
            eop_data=_eop_data,
        )

        models = (_sun_model, _moon_model, srp_model, drag_model)
        model_list = CentralBodyDynamicsModel(_grav_model, models)

        return USM6_EOM(_state_usm6, _p2, _t, model_list)
    end

    for backend in _BACKENDS
        if backend[1] == "Enzyme"
            backend = (
                "Enzyme", AutoEnzyme(; mode=Enzyme.set_runtime_activity(Enzyme.Forward))
            )
        end
        testname = "USM6 Differentiability " * backend[1]
        @testset "$testname" begin
            f_fd, df_fd = value_and_jacobian(
                (x) -> Array(dynamics_params(x)), AutoFiniteDiff(), [0.2; 0.2]
            )

            f_ad, df_ad = value_and_jacobian(
                (x) -> Array(dynamics_params(x)), backend[2], [0.2; 0.2]
            )

            @test f_fd ≈ f_ad
            @test df_fd ≈ df_ad atol = 1e-5
        end
    end
    SpaceIndices.destroy()
end

@testset "USMEM State Differentiability" begin
    SpaceIndices.init()
    for backend in _BACKENDS
        testname = "USMEM Differentiability " * backend[1]
        @testset "$testname" begin
            f_fd, df_fd = value_and_jacobian(
                (x) -> USMEM_EOM(x, _p2, _t, _model_list),
                AutoFiniteDiff(),
                Array(_state_usmem),
            )

            f_ad, df_ad = value_and_jacobian(
                (x) -> Array(USMEM_EOM(x, _p2, _t, _model_list)),
                backend[2],
                Array(_state_usmem),
            )

            @test f_fd ≈ f_ad
            @test df_fd ≈ df_ad atol = 1e-5
        end
    end
    SpaceIndices.destroy()
end

@testset "USMEM Time Differentiability" begin
    SpaceIndices.init()
    for backend in _BACKENDS
        testname = "USMEM Differentiability " * backend[1]
        @testset "$testname" begin
            f_fd, df_fd = value_and_derivative(
                (x) -> USMEM_EOM(Array(_state_usmem), _p2, x, _model_list),
                AutoFiniteDiff(),
                _t,
            )

            f_ad, df_ad = value_and_derivative(
                (x) -> Array(USMEM_EOM(Array(_state_usmem), _p2, x, _model_list)),
                backend[2],
                _t,
            )

            @test f_fd ≈ f_ad
            @test df_fd ≈ df_ad atol = 1e-5
        end
    end
    SpaceIndices.destroy()
end

@testset "USMEM Parameter Differentiability" begin
    SpaceIndices.init()

    function dynamics_params(x::AbstractArray{T}) where {T<:Number}
        satellite_srp_model = CannonballFixedSRP(x[1])
        srp_model = SRPAstroModel(;
            satellite_srp_model=satellite_srp_model,
            sun_data=_sun_model,
            eop_data=_eop_data,
            shadow_model=Conical(),
        )

        satellite_drag_model = CannonballFixedDrag(x[2])
        drag_model = DragAstroModel(;
            satellite_drag_model=satellite_drag_model,
            atmosphere_model=JB2008(),
            eop_data=_eop_data,
        )

        models = (_sun_model, _moon_model, srp_model, drag_model)
        model_list = CentralBodyDynamicsModel(_grav_model, models)

        return USMEM_EOM(_state_usmem, _p2, _t, model_list)
    end

    for backend in _BACKENDS
        if backend[1] == "Enzyme"
            backend = (
                "Enzyme", AutoEnzyme(; mode=Enzyme.set_runtime_activity(Enzyme.Forward))
            )
        end
        testname = "USMEM Differentiability " * backend[1]
        @testset "$testname" begin
            f_fd, df_fd = value_and_jacobian(
                (x) -> Array(dynamics_params(x)), AutoFiniteDiff(), [0.2; 0.2]
            )

            f_ad, df_ad = value_and_jacobian(
                (x) -> Array(dynamics_params(x)), backend[2], [0.2; 0.2]
            )

            @test f_fd ≈ f_ad
            @test df_fd ≈ df_ad atol = 1e-5
        end
    end
    SpaceIndices.destroy()
end
