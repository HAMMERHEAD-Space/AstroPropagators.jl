@testset "Cross-Formulation Consistency (Keplerian, 1 day)" begin
    env = setup_keplerian()
    (; model_list, p, μ, grav_model) = env
    u0 = KEPLERIAN_U0
    tspan = (0.0, 86400.0)
    duration = 86400.0

    ref = run_cowell(u0, model_list, p, tspan)
    @test run_gaussve(u0, model_list, p, tspan) ≈ ref rtol=1e-6
    @test run_milankovich(u0, model_list, p, tspan) ≈ ref rtol=1e-6
    @test run_usm7(u0, model_list, p, tspan) ≈ ref rtol=1e-6
    @test run_usm6(u0, model_list, p, tspan) ≈ ref rtol=1e-6
    @test run_usmem(u0, model_list, p, tspan) ≈ ref rtol=1e-6
    @test run_modeq(u0, model_list, p, tspan) ≈ ref rtol=1e-6
    @test run_edromo(u0, model_list, μ, grav_model, duration) ≈ ref rtol=1e-6
    @test run_edromo(u0, model_list, μ, grav_model, duration; flag_time=ConstantTime()) ≈
        ref rtol=1e-6
    @test run_edromo(u0, model_list, μ, grav_model, duration; flag_time=LinearTime()) ≈ ref rtol=1e-6
    @test run_ks(u0, model_list, μ, grav_model, duration) ≈ ref rtol=1e-6
    @test run_ks(u0, model_list, μ, grav_model, duration; flag_time=LinearTime()) ≈ ref rtol=1e-6
    @test run_stische(u0, model_list, μ, grav_model, duration) ≈ ref rtol=1e-6
    @test run_stische(u0, model_list, μ, grav_model, duration; flag_time=LinearTime()) ≈ ref rtol=1e-6
end

@testset "Cross-Formulation Consistency (HF, 1 day)" begin
    env = setup_full_hf(; srp_coeff=0.2)
    (; model_list, p, μ, grav_model) = env
    u0 = KEPLERIAN_U0
    tspan = (0.0, 86400.0)
    duration = 86400.0

    ref = run_cowell(u0, model_list, p, tspan)
    @test run_gaussve(u0, model_list, p, tspan) ≈ ref rtol=1e-2
    @test run_milankovich(u0, model_list, p, tspan) ≈ ref rtol=1e-3
    @test run_usm7(u0, model_list, p, tspan) ≈ ref rtol=1e-3
    @test run_usm6(u0, model_list, p, tspan) ≈ ref rtol=1e-3
    @test run_usmem(u0, model_list, p, tspan) ≈ ref rtol=1e-3
    @test run_modeq(u0, model_list, p, tspan) ≈ ref rtol=1e-2
    @test run_edromo(u0, model_list, μ, grav_model, duration) ≈ ref rtol=1e-3
    @test run_ks(u0, model_list, μ, grav_model, duration) ≈ ref rtol=1e-3
    @test run_stische(u0, model_list, μ, grav_model, duration) ≈ ref rtol=1e-3
end
