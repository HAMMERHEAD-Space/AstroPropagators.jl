# Round-Trip Closure (RTC) Tests
#
# Forward-propagate for duration T, then reverse (via velocity negation) for T.
# The recovered initial state should match the original to solver tolerance.
# Only valid for conservative (time-independent) force models.
# References: Atallah et al. 2019, Amato et al. 2019 (THALASSA)

# ── Keplerian RTC (all propagators) ──────────────────────────────────
#
# Two-body dynamics are time-independent, so velocity negation gives an
# exact backward trajectory. Tolerance accounts for two propagation legs
# plus root-finding in regularized callbacks.

@testset "Keplerian Round-Trip Closure" begin
    env = setup_keplerian()
    (; model_list, p, μ, grav_model) = env
    u0 = KEPLERIAN_U0
    tspan = (0.0, 86400.0)
    duration = 86400.0

    @testset "Cowell" begin
        mid = run_cowell(u0, model_list, p, tspan)
        recovered = negate_velocity(run_cowell(negate_velocity(mid), model_list, p, tspan))
        @test recovered ≈ u0 rtol=1e-10
    end

    @testset "GaussVE" begin
        mid = run_gaussve(u0, model_list, p, tspan)
        recovered = negate_velocity(run_gaussve(negate_velocity(mid), model_list, p, tspan))
        @test recovered ≈ u0 rtol=1e-8
    end

    @testset "Milankovich" begin
        mid = run_milankovich(u0, model_list, p, tspan)
        recovered = negate_velocity(
            run_milankovich(negate_velocity(mid), model_list, p, tspan)
        )
        @test recovered ≈ u0 rtol=1e-10
    end

    @testset "USM7" begin
        mid = run_usm7(u0, model_list, p, tspan)
        recovered = negate_velocity(run_usm7(negate_velocity(mid), model_list, p, tspan))
        @test recovered ≈ u0 rtol=1e-10
    end

    @testset "USM6" begin
        mid = run_usm6(u0, model_list, p, tspan)
        recovered = negate_velocity(run_usm6(negate_velocity(mid), model_list, p, tspan))
        @test recovered ≈ u0 rtol=1e-10
    end

    @testset "USMEM" begin
        mid = run_usmem(u0, model_list, p, tspan)
        recovered = negate_velocity(run_usmem(negate_velocity(mid), model_list, p, tspan))
        @test recovered ≈ u0 rtol=1e-10
    end

    @testset "EDromo PhysicalTime" begin
        mid = run_edromo(u0, model_list, μ, grav_model, duration)
        recovered = negate_velocity(
            run_edromo(negate_velocity(mid), model_list, μ, grav_model, duration)
        )
        @test recovered ≈ u0 rtol=1e-8
    end

    @testset "EDromo ConstantTime" begin
        mid = run_edromo(u0, model_list, μ, grav_model, duration; flag_time=ConstantTime())
        recovered = negate_velocity(
            run_edromo(
                negate_velocity(mid),
                model_list,
                μ,
                grav_model,
                duration;
                flag_time=ConstantTime(),
            ),
        )
        @test recovered ≈ u0 rtol=1e-8
    end

    @testset "EDromo LinearTime" begin
        mid = run_edromo(u0, model_list, μ, grav_model, duration; flag_time=LinearTime())
        recovered = negate_velocity(
            run_edromo(
                negate_velocity(mid),
                model_list,
                μ,
                grav_model,
                duration;
                flag_time=LinearTime(),
            ),
        )
        @test recovered ≈ u0 rtol=1e-8
    end

    @testset "KS PhysicalTime" begin
        mid = run_ks(u0, model_list, μ, grav_model, duration)
        recovered = negate_velocity(
            run_ks(negate_velocity(mid), model_list, μ, grav_model, duration)
        )
        @test recovered ≈ u0 rtol=1e-8
    end

    @testset "KS LinearTime" begin
        mid = run_ks(u0, model_list, μ, grav_model, duration; flag_time=LinearTime())
        recovered = negate_velocity(
            run_ks(
                negate_velocity(mid),
                model_list,
                μ,
                grav_model,
                duration;
                flag_time=LinearTime(),
            ),
        )
        @test recovered ≈ u0 rtol=1e-8
    end

    @testset "StiSche PhysicalTime" begin
        mid = run_stische(u0, model_list, μ, grav_model, duration)
        recovered = negate_velocity(
            run_stische(negate_velocity(mid), model_list, μ, grav_model, duration)
        )
        @test recovered ≈ u0 rtol=1e-8
    end

    @testset "StiSche LinearTime" begin
        mid = run_stische(u0, model_list, μ, grav_model, duration; flag_time=LinearTime())
        recovered = negate_velocity(
            run_stische(
                negate_velocity(mid),
                model_list,
                μ,
                grav_model,
                duration;
                flag_time=LinearTime(),
            ),
        )
        @test recovered ≈ u0 rtol=1e-8
    end
end

# ── Gravity-Only HF RTC (standard propagators only) ──────────────────
#
# Spherical harmonics + third-body (Sun/Moon), no drag or SRP.
# Uses backward tspan for standard propagators that integrate in physical time.
# Regularized propagators are excluded here since their independent variable
# (fictitious time) does not support backward integration; they are validated
# through cross-formulation consistency with Cowell.

@testset "Gravity-Only HF Round-Trip Closure" begin
    env = setup_gravity_only()
    (; model_list, p, μ) = env
    u0 = KEPLERIAN_U0
    tspan_fwd = (0.0, 86400.0)
    tspan_bwd = (86400.0, 0.0)

    @testset "Cowell" begin
        mid = run_cowell(u0, model_list, p, tspan_fwd)
        recovered = run_cowell(mid, model_list, p, tspan_bwd)
        @test recovered ≈ u0 rtol=1e-8
    end

    @testset "Milankovich" begin
        mid = run_milankovich(u0, model_list, p, tspan_fwd)
        recovered = run_milankovich(mid, model_list, p, tspan_bwd)
        @test recovered ≈ u0 rtol=1e-8
    end

    @testset "USM7" begin
        mid = run_usm7(u0, model_list, p, tspan_fwd)
        recovered = run_usm7(mid, model_list, p, tspan_bwd)
        @test recovered ≈ u0 rtol=1e-8
    end

    @testset "USM6" begin
        mid = run_usm6(u0, model_list, p, tspan_fwd)
        recovered = run_usm6(mid, model_list, p, tspan_bwd)
        @test recovered ≈ u0 rtol=1e-8
    end

    @testset "USMEM" begin
        mid = run_usmem(u0, model_list, p, tspan_fwd)
        recovered = run_usmem(mid, model_list, p, tspan_bwd)
        @test recovered ≈ u0 rtol=1e-8
    end

    @testset "GaussVE" begin
        mid = run_gaussve(u0, model_list, p, tspan_fwd)
        recovered = run_gaussve(mid, model_list, p, tspan_bwd)
        @test recovered ≈ u0 rtol=1e-6
    end
end
