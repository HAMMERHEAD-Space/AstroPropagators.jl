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
