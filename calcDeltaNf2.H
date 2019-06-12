// calculate staggared deltaNf2 field
autoPtr<surfaceScalarField> dNf2Ptr_;
bool calcDNf2(false);

IOobject deltaNf2Header
(
    "deltaNf2",
    runTime.constant(),
    mesh,
    IOobject::NO_READ
);

if(deltaNf2Header.typeHeaderOk<surfaceScalarField>(true))
{
    Info << "Reading grid delta squared field..." << endl;

    dNf2Ptr_.reset
    (
        new surfaceScalarField
        (
            IOobject
            (
                "deltaNf2",
                runTime.constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        )
    );

}
else
{
    Info << "Calculating grid delta squared field..." << endl;
    calcDNf2 = true;

    dNf2Ptr_.reset
    (
        new surfaceScalarField
        (
            IOobject
            (
                "deltaNf2",
                runTime.constant(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimArea, 0)
        )
    );
}

// dereference pointer and assign to a field
surfaceScalarField deltaNf2 = *dNf2Ptr_;


// calculate when necessary
if(calcDNf2)
{
    deltaNf2.primitiveFieldRef() = (mesh.delta() & mesh.delta());

    forAll(mesh.boundaryMesh(), patchI)
    {
        fvsPatchScalarField& Asf = deltaNf2.boundaryFieldRef()[patchI];
        const fvPatch& cPatch = Asf.patch();

        if(cPatch.type() != "empty")
        {
            Asf = cPatch.delta() & cPatch.delta();
        }
    }

    deltaNf2.write();
}
//----------------------------------------------------------------------------//
