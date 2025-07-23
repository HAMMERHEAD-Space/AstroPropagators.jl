@testset "Stiefel-Scheifele State Differentiability" begin
    SpaceIndices.init()

    time_types = (
        ("PhysicalTime", _state_stische_pt, _config_pt, _ϕ_pt),
        ("LinearTime", _state_stische_lt, _config_lt, _ϕ_lt),
    )

    for backend in _BACKENDS
        if backend[1] == "Zygote"
            continue
        end
        for (time_type, state, config, ϕ) in time_types
            testname = "Stiefel-Scheifele Differentiability " * backend[1] * " " * time_type
            @testset "$testname" begin
                f_fd, df_fd = value_and_jacobian(
                    (x) -> StiSche_EOM(x, _p2, ϕ, _model_list, config),
                    AutoFiniteDiff(),
                    state,
                )

                f_ad, df_ad = value_and_jacobian(
                    (x) -> Array(StiSche_EOM(x, _p2, ϕ, _model_list, config)),
                    backend[2],
                    state,
                )

                @test f_fd ≈ f_ad
                @test df_fd ≈ df_ad atol = 1e-3
            end
        end
    end
    SpaceIndices.destroy()
end

@testset "Stiefel-Scheifele Time Differentiability" begin
    SpaceIndices.init()

    time_types = (
        ("PhysicalTime", _state_stische_pt, _config_pt, _ϕ_pt),
        ("LinearTime", _state_stische_lt, _config_lt, _ϕ_lt),
    )

    for backend in _BACKENDS
        if backend[1] == "Zygote"
            continue
        end
        if backend[1] == "Enzyme"
            backend = (
                "Enzyme",
                AutoEnzyme(;
                    mode=Enzyme.set_runtime_activity(Enzyme.Forward),
                    function_annotation=Enzyme.Duplicated,
                ),
            )
        end
        for (time_type, state, config, ϕ) in time_types
            testname =
                "Stiefel-Scheifele Time Differentiability " * backend[1] * " " * time_type
            @testset "$testname" begin
                f_fd, df_fd = value_and_derivative(
                    (x) -> StiSche_EOM(Array(state), _p2, x, _model_list, config),
                    AutoFiniteDiff(),
                    ϕ,
                )

                f_ad, df_ad = value_and_derivative(
                    (x) -> Array(StiSche_EOM(Array(state), _p2, x, _model_list, config)),
                    backend[2],
                    ϕ,
                )

                @test f_fd ≈ f_ad
                @test df_fd ≈ df_ad atol = 1e-3
            end
        end
    end
    SpaceIndices.destroy()
end

@testset "Stiefel-Scheifele Parameter Differentiability" begin
    SpaceIndices.init()

    time_types = (
        ("PhysicalTime", _state_stische_pt, _config_pt, _ϕ_pt),
        ("LinearTime", _state_stische_lt, _config_lt, _ϕ_lt),
    )

    function dynamics_params(
        x::AbstractArray{T}, state::AbstractArray, config, ϕ
    ) where {T<:Number}
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

        return StiSche_EOM(Array(state), _p2, ϕ, model_list, config)
    end

    for backend in _BACKENDS
        if backend[1] == "Zygote"
            continue
        end
        if backend[1] == "Enzyme"
            backend = (
                "Enzyme",
                AutoEnzyme(;
                    mode=Enzyme.set_runtime_activity(Enzyme.Forward),
                    function_annotation=Enzyme.Duplicated,
                ),
            )
        end
        for (time_type, state, config, ϕ) in time_types
            testname =
                "Stiefel-Scheifele Parameter Differentiability " *
                backend[1] *
                " " *
                time_type
            @testset "$testname" begin
                f_fd, df_fd = value_and_jacobian(
                    (x) -> Array(dynamics_params(x, state, config, ϕ)),
                    AutoFiniteDiff(),
                    [0.2; 0.2],
                )

                f_ad, df_ad = value_and_jacobian(
                    (x) -> Array(dynamics_params(x, state, config, ϕ)),
                    backend[2],
                    [0.2; 0.2],
                )

                @test f_fd ≈ f_ad
                @test df_fd ≈ df_ad atol = 1e-3
            end
        end
    end
    SpaceIndices.destroy()
end
