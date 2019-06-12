// Define polynomial laplace filter for regularization
volVectorField polyLaplaceFilter
(
    volVectorField vecField,
    surfaceScalarField deltaNf2,
    dimensionedScalar d1,
    dimensionedScalar d2
)
{
    volVectorField Lap = fvc::laplacian(deltaNf2, vecField, "laplacian(filterField)");
    return (vecField - d1*Lap  + d2*fvc::laplacian(deltaNf2, Lap, "laplacian(filterField)"));
}


surfaceScalarField polyLaplaceFilter
(
    surfaceScalarField fluxField,
    surfaceScalarField deltaNf2,
    dimensionedScalar d1,
    dimensionedScalar d2
)
{
    const fvMesh& mesh = fluxField.mesh();

    volVectorField vecField = fvc::reconstruct(fluxField);
    volVectorField Lap = fvc::laplacian(deltaNf2, vecField, "laplacian(filterField)");

    vecField += (- d1*Lap  + d2*fvc::laplacian(deltaNf2, Lap, "laplacian(filterField)"));
    return (fvc::interpolate(vecField) & mesh.Sf());
}


volScalarField polyLaplaceFilter
(
    volScalarField scalField,
    surfaceScalarField deltaNf2,
    dimensionedScalar d1,
    dimensionedScalar d2
)
{
    volScalarField Lap = fvc::laplacian(deltaNf2, scalField, "laplacian(filterField)");
    return (scalField - d1*Lap  + d2*fvc::laplacian(deltaNf2, Lap, "laplacian(filterField)"));
}
//----------------------------------------------------------------------------//


surfaceScalarField fluxInterpolationCorrection
(
   volVectorField vecField,
   surfaceVectorField faceDistOwner,
   surfaceVectorField faceDistNeighbour
)
{
    const fvMesh& mesh = vecField.mesh();
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const vectorField& Sf = mesh.Sf(); // mesh.faceAreas();

    volTensorField gradField = fvc::grad(vecField);

    surfaceScalarField ssf = (fvc::interpolate(vecField) & mesh.Sf());

    // correction for internal faces only
    forAll(owner, facei)
    {
        ssf[facei] +=
        (0.5 *
            (
                (gradField[owner[facei]] & faceDistOwner[facei])
              - (gradField[neighbour[facei]] & faceDistNeighbour[facei])
            ) & Sf[facei]
        );
    }

    return ssf;
}



surfaceScalarField polyLaplaceFilterFluxCorrection
(
    surfaceScalarField fluxField,
    surfaceScalarField deltaNf2,
    dimensionedScalar d1,
    dimensionedScalar d2,
    surfaceVectorField faceDistOwner,
    surfaceVectorField faceDistNeighbour
)
{
    volVectorField vecField(fvc::reconstruct(fluxField));
    volVectorField Lap = fvc::laplacian(deltaNf2, vecField, "laplacian(filterField)");
    vecField += (- d1*Lap + d2*fvc::laplacian(deltaNf2, Lap, "laplacian(filterField)"));

    return (
        fluxInterpolationCorrection(
            vecField,
            faceDistOwner,
            faceDistNeighbour
        )
    );
}
