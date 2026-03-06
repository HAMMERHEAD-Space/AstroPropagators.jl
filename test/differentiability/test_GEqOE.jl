@testset "GEqOE State Differentiability" begin
    SpaceIndices.init()
    for backend in _BACKENDS
        testname = "GEqOE Differentiability " * backend[1]
        @testset "$testname" begin
            f_fd, df_fd = value_and_jacobian(
                (x) -> GEqOE_EOM(x, _p2, _t, _model_list, _config_geqoe),
                AutoFiniteDiff(),
                Array(_state_geqoe),
            )

            f_ad, df_ad = value_and_jacobian(
                (x) -> Array(GEqOE_EOM(x, _p2, _t, _model_list, _config_geqoe)),
                backend[2],
                Array(_state_geqoe),
            )

            @test f_fd ≈ f_ad
            @test df_fd ≈ df_ad atol = 1e-5
        end
    end
    SpaceIndices.destroy()
end

@testset "GEqOE Time Differentiability" begin
    SpaceIndices.init()
    for backend in _BACKENDS
        testname = "GEqOE Differentiability " * backend[1]
        @testset "$testname" begin
            f_fd, df_fd = value_and_derivative(
                (x) -> GEqOE_EOM(Array(_state_geqoe), _p2, x, _model_list, _config_geqoe),
                AutoFiniteDiff(),
                _t,
            )

            f_ad, df_ad = value_and_derivative(
                (x) -> Array(
                    GEqOE_EOM(Array(_state_geqoe), _p2, x, _model_list, _config_geqoe)
                ),
                backend[2],
                _t,
            )

            @test f_fd ≈ f_ad
            @test df_fd ≈ df_ad atol = 1e-5
        end
    end
    SpaceIndices.destroy()
end

@testset "GEqOE Parameter Differentiability" begin
    SpaceIndices.init()

    function dynamics_params_geqoe(x::AbstractArray{T}) where {T<:Number}
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

        return GEqOE_EOM(_state_geqoe, _p2, _t, model_list, _config_geqoe)
    end

    for backend in _BACKENDS
        if backend[1] == "Enzyme"
            backend = (
                "Enzyme", AutoEnzyme(; mode=Enzyme.set_runtime_activity(Enzyme.Forward))
            )
        end
        testname = "GEqOE Differentiability " * backend[1]
        @testset "$testname" begin
            f_fd, df_fd = value_and_jacobian(
                (x) -> Array(dynamics_params_geqoe(x)), AutoFiniteDiff(), [0.2; 0.2]
            )

            f_ad, df_ad = value_and_jacobian(
                (x) -> Array(dynamics_params_geqoe(x)), backend[2], [0.2; 0.2]
            )

            @test f_fd ≈ f_ad
            @test df_fd ≈ df_ad atol = 1e-5
        end
    end
    SpaceIndices.destroy()
end
