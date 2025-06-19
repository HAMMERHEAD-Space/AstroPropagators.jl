export STI_SCHE_ODE, STI_SCHE_ODE!

function STI_SCHE_ODE(
    du::AbstractArray,
    u::AbstractArray,
    p::ComponentVector,
    ϕ::Number,
    models::AbstractDynamicsModel,
)

    ##################################################
    #* 1. Auxiliary Quantities (1)
    ##################################################
    sϕ2, cϕ2 = sincos(0.5 * ϕ)

    ##################################################
    #* 2. Position and Time in Inertial Frame
    ##################################################    
    KSp = @view(u[1:4]) * cϕ2 + @view(u[5:8]) * sϕ2
    KSv = 0.5 * (-@view(u[1:4]) * sϕ2 + @view(u[5:8]) * cϕ2)

    rV = [
        KSp[1]^2 - KSp[2]^2 - KSp[3]^2 + KS[4]^2
        2.0 * (KSp[1] * KSp[2] - KSp[3] * KSp[4])
        2.0 * (KSp[1] * KSp[3] + KSp[2] * KSp[4])
    ]

    #TODO: IS THIS RIGHT?
    r_mag = dot(rV, rV)

    vV =
        (4.0 * u[9] / r_mag) .* [
            KSp[1] * KSv[1] - KSp[2] * KSv[2] - KSp[3] * KSv[3] + KSp[4] * KSv[4]
            KSp[2] * KSv[1] + KSp[2] * KSv[2] - KSp[4] * KSv[3] - KSp[3] * KSv[4]
            KSp[3] * KSv[1] + KSp[4] * KSv[2] + KSp[1] * KSv[3] + KSp[2] * KSv[4]
        ]

    if p.time_flag == :Physical
        t = u[10]
    elseif p.time_flag == :Linear
        t = u[10] - dot(KSp, KSv) / u[9]
    end

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
    #* 5. Auxiliary Quantities (2)
    ##################################################
    Lp = [
        KSp[1] * P[1] + KSp[2] * P[2] + KSp[3] * P[3]
        -KSp[2] * P[1] + KSp[1] * P[2] + KSp[4] * P[3]
        -KSp[3] * P[1] - KSp[4] * P[2] + KSp[1] * P[3]
        KSp[4] * P[1] - KSp[3] * P[2] + KSp[2] * P[3]
    ]

    ∇Uᵤ =
        -2.0 * [
            KSp[1] * ∇Uᵣ[1] + KSp[2] * ∇Uᵣ[2] + KSp[3] * ∇Uᵣ[3]
            -KSp[2] * ∇Uᵣ[1] + KSp[1] * ∇Uᵣ[2] + KSp[4] * ∇Uᵣ[3]
            -KSp[3] * ∇Uᵣ[1] - KSp[4] * ∇Uᵣ[2] + KSp[1] * ∇Uᵣ[3]
            KSp[4] * ∇Uᵣ[1] - KSp[3] * ∇Uᵣ[2] + KSp[2] * ∇Uᵣ[3]
        ]

    ##################################################
    #* 5. Equations of Motion
    ##################################################
    dh = -r_mag / (8.0 * u[9]^2) * ∇Uₜ - (0.5 / u[9]) * dot(du, Lp)

    aux =
        (0.5 / u[9]^2) * (0.5 * U * KSp + r_mag / 4.0 * (∇Uᵤ - 2.0 * Lp)) +
        2.0 / u[9] * dh * KSv

    if p.time_flag == :Physical
        dt = 0.5 * r_mag / u[9]
    elseif p.time_flag == :Linear
        lte1 = (GE - 2.0 * r_mag * U) / (8.0 * z[9]^3)
        lte2 = (r_mag / (16.0 * u[9]^3)) * dot(KSp, KSv)
        lte3 = (2.0 / u[9]^2) * dh * dot(KSp, KSv)
        dt = lte1 - lte2 - lte3
    end

    du .= [
        aux * sϕ2
        -aux * cϕ2
        dh
        dt
    ]

    return nothing
end

function STI_SCHE_ODE!(
    du::AbstractArray,
    u::AbstractArray,
    p::ComponentVector,
    ϕ::Number,
    models::AbstractDynamicsModel,
)

    ##################################################
    #* 1. Auxiliary Quantities (1)
    ##################################################
    sϕ2, cϕ2 = sincos(0.5 * ϕ)

    ##################################################
    #* 2. Position and Time in Inertial Frame
    ##################################################    
    KSp = @view(u[1:4]) * cϕ2 + @view(u[5:8]) * sϕ2
    KSv = 0.5 * (-@view(u[1:4]) * sϕ2 + @view(u[5:8]) * cϕ2)

    rV = [
        KSp[1]^2 - KSp[2]^2 - KSp[3]^2 + KS[4]^2
        2.0 * (KSp[1] * KSp[2] - KSp[3] * KSp[4])
        2.0 * (KSp[1] * KSp[3] + KSp[2] * KSp[4])
    ]

    #TODO: IS THIS RIGHT?
    r_mag = dot(rV, rV)

    vV =
        (4.0 * u[9] / r_mag) .* [
            KSp[1] * KSv[1] - KSp[2] * KSv[2] - KSp[3] * KSv[3] + KSp[4] * KSv[4]
            KSp[2] * KSv[1] + KSp[2] * KSv[2] - KSp[4] * KSv[3] - KSp[3] * KSv[4]
            KSp[3] * KSv[1] + KSp[4] * KSv[2] + KSp[1] * KSv[3] + KSp[2] * KSv[4]
        ]

    if p.time_flag == :Physical
        t = u[10]
    elseif p.time_flag == :Linear
        t = u[10] - dot(KSp, KSv) / u[9]
    end

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
    #* 5. Auxiliary Quantities (2)
    ##################################################
    Lp = [
        KSp[1] * P[1] + KSp[2] * P[2] + KSp[3] * P[3]
        -KSp[2] * P[1] + KSp[1] * P[2] + KSp[4] * P[3]
        -KSp[3] * P[1] - KSp[4] * P[2] + KSp[1] * P[3]
        KSp[4] * P[1] - KSp[3] * P[2] + KSp[2] * P[3]
    ]

    ∇Uᵤ =
        -2.0 * [
            KSp[1] * ∇Uᵣ[1] + KSp[2] * ∇Uᵣ[2] + KSp[3] * ∇Uᵣ[3]
            -KSp[2] * ∇Uᵣ[1] + KSp[1] * ∇Uᵣ[2] + KSp[4] * ∇Uᵣ[3]
            -KSp[3] * ∇Uᵣ[1] - KSp[4] * ∇Uᵣ[2] + KSp[1] * ∇Uᵣ[3]
            KSp[4] * ∇Uᵣ[1] - KSp[3] * ∇Uᵣ[2] + KSp[2] * ∇Uᵣ[3]
        ]

    ##################################################
    #* 5. Equations of Motion
    ##################################################
    dh = -r_mag / (8.0 * u[9]^2) * ∇Uₜ - (0.5 / u[9]) * dot(du, Lp)

    aux =
        (0.5 / u[9]^2) * (0.5 * U * KSp + r_mag / 4.0 * (∇Uᵤ - 2.0 * Lp)) +
        2.0 / u[9] * dh * KSv

    if p.time_flag == :Physical
        dt = 0.5 * r_mag / u[9]
    elseif p.time_flag == :Linear
        lte1 = (GE - 2.0 * r_mag * U) / (8.0 * z[9]^3)
        lte2 = (r_mag / (16.0 * u[9]^3)) * dot(KSp, KSv)
        lte3 = (2.0 / u[9]^2) * dh * dot(KSp, KSv)
        dt = lte1 - lte2 - lte3
    end

    du .= [
        aux * sϕ2
        -aux * cϕ2
        dh
        dt
    ]

    return nothing
end
