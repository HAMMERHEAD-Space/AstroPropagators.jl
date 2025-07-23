export impulsive_burn!
"""
    function impulsive_burn!(
        integrator::SciMLBase.DEIntegrator, 
        ΔV::AbstractVector; 
        coordinate_set::AstroCoords.AstroCoord=Cartesian
    )

Computes new state from an input impulsive burn. The supplied burn should be in the inertial frame
and has to be computed before-hand. Both will change in future iterations.

Arguments:
-`integrator::SciMLBase.DEIntegrator`: The differential equation integrator object.
-`ΔV::AbstractVector`: The deltaV of the impulsive burn.
-`coordinate_set::AstroCoords.AstroCoord=Cartesian`: The coordinate set the propagation is occurring in.
Returns:
-`nothing`
"""
function impulsive_burn!(
    integrator::T, ΔV::AbstractVector; coordinate_set::V=Cartesian
) where {T<:SciMLBase.DEIntegrator,V<:typeof(AstroCoords.AstroCoord)}
    cart_state = Cartesian(coordinate_set(integrator.u), integrator.p.μ)
    new_state = cart_state + SVector{6}(0, 0, 0, ΔV[1], ΔV[2], ΔV[3])
    new_cart_state = Cartesian(new_state...)

    integrator.u = params(coordinate_set(new_cart_state, integrator.p.μ))

    return nothing
end
