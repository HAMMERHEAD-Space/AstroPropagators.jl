function EDROMO_ODE!(
    du::AbstractArray,
    u::AbstractArray,
    p::ComponentVector,
    ϕ::Number,
    models::NTuple{N,AstroForceModels.AbstractAstroForceModel},
) where {N}
    ζ1, ζ2, ζ3, ζ4, ζ5, ζ6, ζ7, ζ8 = u

    ##################################################
    #* 1. Auxiliary Quantities (1)
    ##################################################
    sϕ, cϕ = sincos(ϕ)

    ρ = 1.0 - ζ1 * cϕ - ζ2 * sϕ
    r_mag = ζ3 * ρ
    χ = ζ1 * sϕ - ζ2 * cϕ
    ϵ = √(1.0 - ζ1^2 - ζ2^2)

    cν = (cϕ - ζ1 + (χ * ζ2) / (ϵ + 1.0)) / ρ
    sν = (sϕ - ζ2 + (χ * ζ1) / (ϵ + 1.0)) / ρ

    ##################################################
    #* 2. Position and Time in Inertial Frame
    ##################################################    
    x = 2.0 .* [
        0.5 - ζ5^2 - ζ6^2
        ζ4 * ζ5 + ζ6 * ζ7
        ζ4 * ζ6 - ζ5 * ζ7
    ]

    y = 2.0 .* [
        ζ4 * ζ5 - ζ6 * ζ7
        0.5 - ζ4^2 - ζ6^2
        ζ5 * ζ6 + ζ4 * ζ7
    ]

    rV = r_mag * (x * cν + y * sν)

    if time_flag == :Physical
        tt = ζ8
    elseif time_flag == :Constant
        tt = ζ8 - ζ3^1.5 * (χ - ϕ)
    elseif time_flag == :Linear
        tt = ζ8 - ζ3^1.5 * χ
    end

    ##################################################
    #* 3. Potential Based Perturbations
    ##################################################    

    ##################################################
    #* 4. Velocity in the Interital Frame
    ##################################################    
    i_vec = x * cν + y * sν
    j_vec = y * cν - x * sν

    v_rad = χ / (√(ζ3) * ρ)
    v_tan = √((1.0 - ζ1^2 - ζ2^2) / (ζ3 * ρ^2) - 2.0 * U)

    vV = v_rad * i_vec + v_tan * j_vec
    v = norm([v_rad; v_tan])

    #TODO: NEEDED?
    cg = v_rad / v
    sg = v_tan / v

    ##################################################
    #* 5. Non-Potential Based Perturbations
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
    #* 6. Auxiliary Quantities (2)
    ##################################################    
    F = P * ∇Uᵣ

    η = √(ϵ^2 - 2.0 * ζ3 * ρ^2 * U)
    aux0 = P[1] * χ + P[2] * η + ∇Uₜ * √(ζ3 * ρ)

    dζ3 = 2.0 * ζ3^3 * aux0
    L3 = dζ3 / (2.0 * ζ3)

    aux1 = ((2.0 * U - F[1] * r_mag) * (2.0 - ρ + ϵ) * r_mag) / (ϵ * (ϵ + 1.0))
    aux2 = (L3 * χ * (ρ - ϵ)) / (ϵ * (ϵ + 1.0))
    wz = (η - ϵ) / ρ + aux1 + aux2

    ##################################################
    #* 7. Equations of Motion
    ##################################################   

    aux3 = (F[1] * r_mag - 2.0 * U) * r_mag
    aux4 = (F[3] * r_mag^2) / (2.0 * η)

    if time_flag == :Physical
        dt = √(ζ3) * r_mag
    elseif time_flag == :Constant
        dt = ζ3^1.5 * (aux3 + (χ - 1.5 * ϕ) * dζ3 / ζ3)
    elseif time_flag == :Linear
        dt = ζ3^1.5 * (1.0 + aux3 + 2.0 * L3 * χ)
    end

    du .= [
        aux3 * sϕ + L3 * (1.0 + ρ) * cϕ - ζ1
        -aux3 * cϕ + L3 * (1.0 + ρ) * sϕ - ζ2
        dζ3
        aux4 * (ζ7 * cν - ζ6 * sν) + 0.5 * wz * ζ5
        aux4 * (ζ6 * cν + ζ7 * sν) - 0.5 * wz * ζ4
        -aux4 * (ζ5 * cν - ζ4 * sν) + 0.5 * wz * ζ7
        -aux4 * (ζ4 * cν + ζ5 * sν) - 0.5 * wz * ζ6
        dt
    ]

    return nothing
end
