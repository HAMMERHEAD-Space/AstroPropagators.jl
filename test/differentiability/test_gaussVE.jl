@testset "Gauss VE State Differentiability" begin
    SpaceIndices.init()
    for backend in _BACKENDS
        testname = "Gauss VE Differentiability " * backend[1]
        @testset "$testname" begin
            f_fd, df_fd = value_and_jacobian(
                (x) -> GaussVE_EOM(x, _p2, _t, _model_list),
                AutoFiniteDiff(),
                Array(_state_koe),
            )

            f_ad, df_ad = value_and_jacobian(
                (x) -> Array(GaussVE_EOM(x, _p2, _t, _model_list)),
                backend[2],
                Array(_state_koe),
            )

            @test f_fd ≈ f_ad
            @test df_fd ≈ df_ad atol = 1e-5
        end
    end
    SpaceIndices.destroy()
end

@testset "Gauss VE Time Differentiability" begin
    SpaceIndices.init()
    for backend in _BACKENDS
        testname = "Gauss VE Differentiability " * backend[1]
        @testset "$testname" begin
            f_fd, df_fd = value_and_derivative(
                (x) -> GaussVE_EOM(Array(_state_koe), _p2, x, _model_list),
                AutoFiniteDiff(),
                _t,
            )

            f_ad, df_ad = value_and_derivative(
                (x) -> Array(GaussVE_EOM(Array(_state_koe), _p2, x, _model_list)),
                backend[2],
                _t,
            )

            @test f_fd ≈ f_ad
            @test df_fd ≈ df_ad atol = 1e-5
        end
    end
    SpaceIndices.destroy()
end

@testset "Gauss VE Parameter Differentiability" begin
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

        return GaussVE_EOM(_state_koe, _p2, _t, model_list)
    end

    for backend in _BACKENDS
        if backend[1] == "Enzyme"
            backend = (
                "Enzyme", AutoEnzyme(; mode=Enzyme.set_runtime_activity(Enzyme.Forward))
            )
        end
        testname = "Gauss VE Differentiability " * backend[1]
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
