function KS_ODE!(
    du::AbstractArray,
    u::AbstractArray,
    p::ComponentVector,
    ϕ::Number,
    models::AbstractDynamicsModel,
)
    KSp1, KSp2, KSp3, KSp4, KSv1, KSv2, KSv3, KSv4, h, τ = u

    ##################################################
    #* 1. Cartesian Coordinates
    ##################################################

    ##################################################
    #* 2. Potential Based Perturbations
    ##################################################

    ##################################################
    #* 3. Potential Based Perturbations
    ##################################################

    ##################################################
    #* 4. Non-Potential Based Perturbations
    ################################0##################    
    P = non_potential_accel(
        [rV; vV],
        p,
        tt,
        eop_data;
        drag_model=drag_model,
        atmosphere_type=atmosphere_type,
        srp_model=srp_model,
        shadow_model=shadow_model,
        lunar_3rd_body=lunar_3rd_body,
        solar_3rd_body=solar_3rd_body,
    )

    ##################################################
    #* 5. Equations of Motion
    ################################0##################    
    F = P * ∇Uᵣ

    L = KS_matrix()

    dh = -r_mag * ∇Uₜ - 2.0 * dot(@view(u[5:8]), L' * P)

    if time_flag == :Sundman
        dt = r_mag
    elseif time_flag == :Linear
        lte1 = (GE - 2.0 * r_mag * U) / (2.0 * h)
        lte20 = 2.0 * (L' * -F)
        lte2 = (r_mag / (4.0 * h)) * dot(@view(u[1:4]), lte20)
        lte3 = dh / (h^2) * dot(@view(u[1:4]), @view(u[5:8]))
        dt = lte1 - lte2 - lte3
    end

    du .= [
        @view(u[5:8])
        -0.5 * (@view(u[1:4]) * (h + U) - r_mag * (L' * F))
        dh
        dt
    ]

    return nothing
end

function KS_matrix(u::AbstractArray)
    return @SMatrix{4, 4}(
        [
            u[1] u[2] u[3] u[4]
            -u[2] u[1] u[4] -u[3]
            -u[3] -u[4] u[1] u[2]
            u[4] -u[3] u[2] -u[1]
        ]
    )
end

function KS2cart(u::AbstractArray)
    x = u[1]^2 - u[2]^2 - u[3]^2 + u[4]^2
    y = 2.0 * (u[1] * u[2] - u[3] * u[4])
    z = 2.0 * (u[1] * u[3] + u[2] * u[4])
    r = norm(@view(u[1:4]))

    ẋ = 2.0 * (u[1] * u[5] - u[2] * u[6] - u[3] * u[7] + u[4] * u[8]) / r
    ẏ = 2.0 * (u[2] * u[5] + u[1] * u[6] - u[4] * u[7] - u[3] * u[8]) / r
    ż = 2.0 * (u[3] * u[5] + u[4] * u[6] + u[1] * u[7] + u[2] * u[8]) / r

    return [x; y; z; ẋ; ẏ; ż]
end

function cart2KS(
    u::AbstractArray,
    t::Number,
    μ::Number,
    DU::Number,
    TU::Number,
    u0::AbstractArray,
    U::Number,
    time_flag::Symbol,
)

    ##################################################
    #* 1. Initialize Position & Velocity in ℛ⁴, Non-Dimensional
    ##################################################   
    x0 = [@view(u[1:3]) / DU; 0.0]
    r = norm(x0)
    #TODO: IS THIS RIGHT?
    t0 = t * TU

    v0 = [@view(u[4:6]) / (DU * TU); 0.0]

    U0 = U / ((DU * TU)^2)
    KSq = μ / (DU^3 * TU^2)

    ##################################################
    #* 2. Initialize KS-Position & KS-Velocity
    ##################################################   
    if (x0[1] ≥ 0.0)
        KSp1 = 0.0
        KSp4 = √(0.5 * (r + x0[1]) - KSp1^2)
        KSp2 = (x0[2] * KSp1 + x0[3] * KSp4) / (r + x0[1])
        KSp3 = (x0[3] * KSp1 - x0[2] * KSp4) / (r + x0[1])
    else
        KSp2 = 0.0
        KSp3 = √(0.5 * (r - x0[1]) - KSp2^2)
        KSp1 = (x0[2] * KSp2 + x0[3] * KSp4) / (r + x0[1])
        KSp4 = (x0[3] * KSp2 - x0[2] * KSp3) / (r + x0[1])
    end

    KSv1 = 0.5 * dot([KSp1; KSp2; KSp3], @view(v0[1:3]))
    KSv2 = 0.5 * dot([-KSp2; KSp1; KSp4], @view(v0[1:3]))
    KSv3 = 0.5 * dot([-KSp3; -KSp4; KSp1], @view(v0[1:3]))
    KSv4 = 0.5 * dot([KSp4; -KSp3; KSp2], @view(v0[1:3]))

    ##################################################
    #* 3. Total Energy
    ##################################################   
    K = 0.5 * norm(v0)^2
    h = -KSq / r - K - U0

    ##################################################
    #* 4. Initial Time
    ##################################################   
    if time_flag == :Sundman
        τ = t0
    elseif time_flag == :Linear
        τ = t0 + dot([KSp1; KSp2; KSp3; KSp4], [KSv1; KSv2; KSv3; KSv4]) / h
    end

    return [KSp1; KSp2; KSp3; KSp4; KSv1; KSv2; KSv3; KSv4; h; τ]
end
