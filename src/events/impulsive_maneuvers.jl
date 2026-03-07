export impulsive_burn!
"""
    impulsive_burn!(
        integrator::SciMLBase.DEIntegrator,
        ΔV::AbstractVector;
        coordinate_set::AstroCoords.AstroCoord=Cartesian,
        frame::AbstractThrustFrame=InertialFrame(),
    )

Apply an instantaneous velocity change (impulsive burn) to the integrator state.

The `ΔV` vector is interpreted in the reference frame specified by `frame`:

- `InertialFrame()`: `ΔV` components are inertial (ECI) — no rotation needed (default)
- `RTNFrame()`: `ΔV = [ΔV_R, ΔV_T, ΔV_N]` in the Radial–Transverse–Normal frame
- `VNBFrame()`: `ΔV = [ΔV_V, ΔV_N, ΔV_B]` in the Velocity–Normal–Binormal frame

# Arguments
- `integrator::SciMLBase.DEIntegrator`: The differential equation integrator object.
- `ΔV::AbstractVector`: The 3-component delta-V of the impulsive burn [km/s].

# Keyword Arguments
- `coordinate_set`: The coordinate set the propagation is using (default: `Cartesian`).
- `frame::AbstractThrustFrame`: Reference frame in which `ΔV` is expressed (default: `InertialFrame()`).
"""
function impulsive_burn!(
    integrator::T,
    ΔV::AbstractVector;
    coordinate_set::V=Cartesian,
    frame::AbstractThrustFrame=InertialFrame(),
) where {T<:SciMLBase.DEIntegrator,V<:typeof(AstroCoords.AstroCoord)}
    cart_state = Cartesian(coordinate_set(integrator.u), integrator.p.μ)
    ΔV_inertial = transform_to_inertial(SVector{3}(ΔV[1], ΔV[2], ΔV[3]), cart_state, frame)
    new_state =
        cart_state + SVector{6}(0, 0, 0, ΔV_inertial[1], ΔV_inertial[2], ΔV_inertial[3])
    new_cart_state = Cartesian(new_state...)

    integrator.u = params(coordinate_set(new_cart_state, integrator.p.μ))

    return nothing
end
