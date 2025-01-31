@testset "Milankovich State Differentiability" begin
    SpaceIndices.init()
    for backend in _BACKENDS
        testname = "Milankovich Differentiability " * backend[1]
        @testset "$testname" begin
            f_fd, df_fd = value_and_jacobian(
                (x) -> Milankovich_EOM(x, _p2, _t, _model_list),
                AutoFiniteDiff(),
                Array(_state_mil),
            )

            f_ad, df_ad = value_and_jacobian(
                (x) -> Array(Milankovich_EOM(x, _p2, _t, _model_list)),
                backend[2],
                Array(_state_mil),
            )

            @test f_fd ≈ f_ad
            @test df_fd ≈ df_ad atol = 1e-5
        end
    end
    SpaceIndices.destroy()
end

@testset "Milankovich Time Differentiability" begin
    SpaceIndices.init()
    for backend in _BACKENDS
        testname = "Milankovich Differentiability " * backend[1]
        @testset "$testname" begin
            f_fd, df_fd = value_and_derivative(
                (x) -> Milankovich_EOM(Array(_state_mil), _p2, x, _model_list),
                AutoFiniteDiff(),
                _t,
            )

            f_ad, df_ad = value_and_derivative(
                (x) -> Array(Milankovich_EOM(Array(_state_mil), _p2, x, _model_list)),
                backend[2],
                _t,
            )

            @test f_fd ≈ f_ad
            @test df_fd ≈ df_ad atol = 1e-5
        end
    end
    SpaceIndices.destroy()
end
