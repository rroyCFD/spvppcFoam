// explicit convection function
volVectorField explConvection
(
    surfaceScalarField ssf_,
    volVectorField vvf_
)
{
    // cross product of the convection operation
    return (fvc::div(ssf_, vvf_, "div(phi,U)") - fvc::div(ssf_)*vvf_ );
}


volVectorField explConvection
(
    volVectorField vvfFlux_,
    volVectorField vvf_
)
{
    const fvMesh& mesh = vvfFlux_.mesh();

    return explConvection( (fvc::interpolate(vvfFlux_) & mesh.Sf()), vvf_);
}