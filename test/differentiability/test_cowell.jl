@testset "Cowell State Differntiability" begin
    SpaceIndices.init()
    for backend in _BACKENDS
        testname = "Cowell Differntiability " * backend[1]
        @testset "$testname" begin
            f_fd, df_fd = value_and_jacobian(
                (x) -> Cowell_EOM(x, _p, _t, _model_list), AutoFiniteDiff(), _state
            )

            f_ad, df_ad = value_and_jacobian(
                (x) -> Array(Cowell_EOM(x, _p, _t, _model_list)), backend[2], _state
            )

            @test f_fd ≈ f_ad
            @test df_fd ≈ df_ad atol = 1e-5
        end
    end
    SpaceIndices.destroy()
end

@testset "Cowell Time Differntiability" begin
    SpaceIndices.init()
    for backend in _BACKENDS
        testname = "Cowell Differntiability " * backend[1]
        @testset "$testname" begin
            f_fd, df_fd = value_and_derivative(
                (x) -> Cowell_EOM(_state, _p, x, _model_list), AutoFiniteDiff(), _t
            )

            f_ad, df_ad = value_and_derivative(
                (x) -> Array(Cowell_EOM(_state, _p, x, _model_list)), backend[2], _t
            )

            @test f_fd ≈ f_ad
            @test df_fd ≈ df_ad atol = 1e-5
        end
    end
    SpaceIndices.destroy()
end
