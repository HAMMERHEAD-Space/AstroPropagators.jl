@testset "eom / eom! dispatch" begin
    env = setup_keplerian()
    (; model_list, p, μ, grav_model) = env
    u0 = KEPLERIAN_U0
    t = 0.0

    @testset "Standard — eom vs eom!" begin
        for (prop, label) in (
            (CowellPropagator(), "Cowell"),
            (GaussVEPropagator(), "GaussVE"),
            (MilankovichPropagator(), "Milankovich"),
            (USM7Propagator(), "USM7"),
            (USM6Propagator(), "USM6"),
            (USMEMPropagator(), "USMEM"),
        )
            @testset "$label" begin
                u = if prop isa CowellPropagator
                    u0
                elseif prop isa GaussVEPropagator
                    Array(Keplerian(Cartesian(u0), μ))
                elseif prop isa MilankovichPropagator
                    Array(Milankovich(Cartesian(u0), μ))
                elseif prop isa USM7Propagator
                    Array(USM7(Cartesian(u0), μ))
                elseif prop isa USM6Propagator
                    Array(USM6(Cartesian(u0), μ))
                else
                    Array(USMEM(Cartesian(u0), μ))
                end
                p_vec = ComponentVector(; JD=COMMON_JD, μ=μ)
                du_oop = eom(prop, u, p_vec, t, model_list)
                du_ip = similar(u)
                eom!(prop, du_ip, u, p_vec, t, model_list)
                @test collect(du_oop) ≈ du_ip
            end
        end
    end

    @testset "Regularized — eom vs eom!" begin
        W = _regularized_W(u0, grav_model, μ)
        p_vec = ComponentVector(; JD=COMMON_JD, μ=μ)

        @testset "EDromo" begin
            config = RegularizedCoordinateConfig(
                u0, μ; W=W, t₀=0.0, flag_time=PhysicalTime()
            )
            ϕ₀ = compute_initial_phi(u0, μ, config)
            u = Array(EDromo(Cartesian(u0), μ, ϕ₀, config))
            du_oop = eom(EDromoPropagator(), u, p_vec, ϕ₀, model_list, config)
            du_ip = similar(u)
            eom!(EDromoPropagator(), du_ip, u, p_vec, ϕ₀, model_list, config)
            @test collect(du_oop) ≈ du_ip
        end

        @testset "KS" begin
            config = RegularizedCoordinateConfig(
                u0, μ; W=W, t₀=0.0, flag_time=PhysicalTime()
            )
            u = Array(KustaanheimoStiefel(Cartesian(u0), μ, config))
            du_oop = eom(KSPropagator(), u, p_vec, 0.0, model_list, config)
            du_ip = similar(u)
            eom!(KSPropagator(), du_ip, u, p_vec, 0.0, model_list, config)
            @test collect(du_oop) ≈ du_ip
        end

        @testset "StiSche" begin
            config = RegularizedCoordinateConfig(
                u0, μ; W=W, t₀=0.0, flag_time=PhysicalTime()
            )
            ϕ₀ = compute_initial_phi(u0, μ, config)
            u = Array(StiefelScheifele(Cartesian(u0), μ, ϕ₀, config))
            du_oop = eom(StiSchePropagator(), u, p_vec, ϕ₀, model_list, config)
            du_ip = similar(u)
            eom!(StiSchePropagator(), du_ip, u, p_vec, ϕ₀, model_list, config)
            @test collect(du_oop) ≈ du_ip
        end
    end
end

@testset "propagate API — Standard Propagators" begin
    env = setup_keplerian()
    (; model_list, p, μ) = env
    u0 = KEPLERIAN_U0
    tspan = (0.0, 86400.0)

    ref = run_cowell(u0, model_list, p, tspan)

    @testset "CowellPropagator" begin
        sol = propagate(CowellPropagator(), u0, p, model_list, tspan)
        @test sol.u[end] ≈ ref
    end

    @testset "GaussVEPropagator" begin
        u0_koe = Array(Keplerian(Cartesian(u0), μ))
        sol = propagate(GaussVEPropagator(), u0_koe, p, model_list, tspan)
        @test Array(Cartesian(Keplerian(sol.u[end]), μ)) ≈ ref
    end

    @testset "MilankovichPropagator" begin
        u0_mil = Array(Milankovich(Cartesian(u0), μ))
        sol = propagate(MilankovichPropagator(), u0_mil, p, model_list, tspan)
        @test Array(Cartesian(Milankovich(sol.u[end]), μ)) ≈ ref
    end

    @testset "USM7Propagator" begin
        u0_usm = Array(USM7(Cartesian(u0), μ))
        sol = propagate(USM7Propagator(), u0_usm, p, model_list, tspan)
        @test Array(Cartesian(USM7(sol.u[end]), μ)) ≈ ref
    end

    @testset "USM6Propagator" begin
        u0_usm = Array(USM6(Cartesian(u0), μ))
        sol = propagate(USM6Propagator(), u0_usm, p, model_list, tspan)
        @test Array(Cartesian(USM6(sol.u[end]), μ)) ≈ ref rtol = 1e-3
    end

    @testset "USMEMPropagator" begin
        u0_usm = Array(USMEM(Cartesian(u0), μ))
        sol = propagate(USMEMPropagator(), u0_usm, p, model_list, tspan)
        @test Array(Cartesian(USMEM(sol.u[end]), μ)) ≈ ref rtol = 1e-3
    end

    @testset "kwargs forwarding" begin
        sol = propagate(
            CowellPropagator(),
            u0,
            p,
            model_list,
            tspan;
            solver=Vern9(),
            abstol=1e-10,
            reltol=1e-10,
            saveat=3600.0,
        )
        @test length(sol.u) == 25
    end
end

@testset "propagate API — Regularized Propagators" begin
    env = setup_keplerian()
    (; model_list, p, μ, grav_model) = env
    u0 = KEPLERIAN_U0
    duration = 86400.0

    ref = run_cowell(u0, model_list, p, (0.0, duration))

    W = _regularized_W(u0, grav_model, μ)

    @testset "EDromoPropagator" begin
        config = RegularizedCoordinateConfig(u0, μ; W=W, t₀=0.0, flag_time=PhysicalTime())
        ϕ₀ = compute_initial_phi(u0, μ, config)
        u0_ed = Array(EDromo(Cartesian(u0), μ, ϕ₀, config))
        tspan = (ϕ₀, ϕ₀ + 6π)
        sol = propagate(
            EDromoPropagator(),
            u0_ed,
            p,
            model_list,
            tspan,
            config;
            callback=build_termination_callback(duration, EDromo, config),
        )
        final_cart = Array(Cartesian(EDromo(sol.u[end]), μ, sol.t[end], config))
        @test final_cart ≈ ref rtol = 1e-3
    end

    @testset "KSPropagator" begin
        config = RegularizedCoordinateConfig(u0, μ; W=W, t₀=0.0, flag_time=PhysicalTime())
        u0_ks = Array(KustaanheimoStiefel(Cartesian(u0), μ, config))
        tspan = (0.0, 9π)
        sol = propagate(
            KSPropagator(),
            u0_ks,
            p,
            model_list,
            tspan,
            config;
            callback=build_termination_callback(duration, KustaanheimoStiefel, config),
        )
        final_cart = Array(Cartesian(KustaanheimoStiefel(sol.u[end]), μ, config))
        @test final_cart ≈ ref rtol = 1e-3
    end

    @testset "StiSchePropagator" begin
        config = RegularizedCoordinateConfig(u0, μ; W=W, t₀=0.0, flag_time=PhysicalTime())
        ϕ₀ = compute_initial_phi(u0, μ, config)
        u0_ss = Array(StiefelScheifele(Cartesian(u0), μ, ϕ₀, config))
        tspan = (ϕ₀, ϕ₀ + 6π)
        sol = propagate(
            StiSchePropagator(),
            u0_ss,
            p,
            model_list,
            tspan,
            config;
            callback=build_termination_callback(duration, StiefelScheifele, config),
        )
        final_cart = Array(Cartesian(StiefelScheifele(sol.u[end]), μ, sol.t[end], config))
        @test final_cart ≈ ref rtol = 1e-3
    end
end

@testset "propagate! API — Standard Propagators" begin
    env = setup_keplerian()
    (; model_list, p, μ) = env
    u0 = KEPLERIAN_U0
    tspan = (0.0, 86400.0)

    ref = run_cowell(u0, model_list, p, tspan)

    @testset "CowellPropagator" begin
        sol = propagate!(CowellPropagator(), u0, p, model_list, tspan)
        @test sol.u[end] ≈ ref
    end

    @testset "GaussVEPropagator" begin
        u0_koe = Array(Keplerian(Cartesian(u0), μ))
        sol = propagate!(GaussVEPropagator(), u0_koe, p, model_list, tspan)
        @test Array(Cartesian(Keplerian(sol.u[end]), μ)) ≈ ref
    end

    @testset "MilankovichPropagator" begin
        u0_mil = Array(Milankovich(Cartesian(u0), μ))
        sol = propagate!(MilankovichPropagator(), u0_mil, p, model_list, tspan)
        @test Array(Cartesian(Milankovich(sol.u[end]), μ)) ≈ ref
    end

    @testset "USM7Propagator" begin
        u0_usm = Array(USM7(Cartesian(u0), μ))
        sol = propagate!(USM7Propagator(), u0_usm, p, model_list, tspan)
        @test Array(Cartesian(USM7(sol.u[end]), μ)) ≈ ref
    end

    @testset "USM6Propagator" begin
        u0_usm = Array(USM6(Cartesian(u0), μ))
        sol = propagate!(USM6Propagator(), u0_usm, p, model_list, tspan)
        @test Array(Cartesian(USM6(sol.u[end]), μ)) ≈ ref rtol = 1e-3
    end

    @testset "USMEMPropagator" begin
        u0_usm = Array(USMEM(Cartesian(u0), μ))
        sol = propagate!(USMEMPropagator(), u0_usm, p, model_list, tspan)
        @test Array(Cartesian(USMEM(sol.u[end]), μ)) ≈ ref rtol = 1e-3
    end
end

@testset "propagate! API — Regularized Propagators" begin
    env = setup_keplerian()
    (; model_list, p, μ, grav_model) = env
    u0 = KEPLERIAN_U0
    duration = 86400.0

    ref = run_cowell(u0, model_list, p, (0.0, duration))

    W = _regularized_W(u0, grav_model, μ)

    @testset "EDromoPropagator" begin
        config = RegularizedCoordinateConfig(u0, μ; W=W, t₀=0.0, flag_time=PhysicalTime())
        ϕ₀ = compute_initial_phi(u0, μ, config)
        u0_ed = Array(EDromo(Cartesian(u0), μ, ϕ₀, config))
        tspan = (ϕ₀, ϕ₀ + 6π)
        sol = propagate!(
            EDromoPropagator(),
            u0_ed,
            p,
            model_list,
            tspan,
            config;
            callback=build_termination_callback(duration, EDromo, config),
        )
        final_cart = Array(Cartesian(EDromo(sol.u[end]), μ, sol.t[end], config))
        @test final_cart ≈ ref rtol = 1e-3
    end

    @testset "KSPropagator" begin
        config = RegularizedCoordinateConfig(u0, μ; W=W, t₀=0.0, flag_time=PhysicalTime())
        u0_ks = Array(KustaanheimoStiefel(Cartesian(u0), μ, config))
        tspan = (0.0, 9π)
        sol = propagate!(
            KSPropagator(),
            u0_ks,
            p,
            model_list,
            tspan,
            config;
            callback=build_termination_callback(duration, KustaanheimoStiefel, config),
        )
        final_cart = Array(Cartesian(KustaanheimoStiefel(sol.u[end]), μ, config))
        @test final_cart ≈ ref rtol = 1e-3
    end

    @testset "StiSchePropagator" begin
        config = RegularizedCoordinateConfig(u0, μ; W=W, t₀=0.0, flag_time=PhysicalTime())
        ϕ₀ = compute_initial_phi(u0, μ, config)
        u0_ss = Array(StiefelScheifele(Cartesian(u0), μ, ϕ₀, config))
        tspan = (ϕ₀, ϕ₀ + 6π)
        sol = propagate!(
            StiSchePropagator(),
            u0_ss,
            p,
            model_list,
            tspan,
            config;
            callback=build_termination_callback(duration, StiefelScheifele, config),
        )
        final_cart = Array(Cartesian(StiefelScheifele(sol.u[end]), μ, sol.t[end], config))
        @test final_cart ≈ ref rtol = 1e-3
    end
end

@testset "propagate vs propagate! consistency" begin
    env = setup_keplerian()
    (; model_list, p) = env
    u0 = KEPLERIAN_U0
    tspan = (0.0, 86400.0)

    sol_oop = propagate(CowellPropagator(), u0, p, model_list, tspan)
    sol_ip = propagate!(CowellPropagator(), u0, p, model_list, tspan)
    @test sol_oop.u[end] ≈ sol_ip.u[end]
end
