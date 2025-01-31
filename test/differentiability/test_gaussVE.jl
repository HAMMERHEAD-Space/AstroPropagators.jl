@testset "Gauss VE State Differntiability" begin
    SpaceIndices.init()
    for backend in _BACKENDS
        testname = "Gauss VE Differntiability " * backend[1]
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

@testset "Gauss VE Time Differntiability" begin
    SpaceIndices.init()
    for backend in _BACKENDS
        testname = "Gauss VE Differntiability " * backend[1]
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
