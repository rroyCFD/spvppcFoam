// filter the extrapolated pressure
{
    Info << "Filtering the extrapolated pressure" << endl;
    volScalarField pef
    (
        IOobject
        (
            "pef",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        C4Regularization.polyLaplaceFilter(pe),
        p.boundaryField().types()
    );
    pef.correctBoundaryConditions();

    volScalarField pePrime
    (
        IOobject
        (
            "pePrime",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (pe - pef),
        p.boundaryField().types()
    );
    pePrime.correctBoundaryConditions();


    pe = pef + C4Regularization.polyLaplaceFilter(pePrime);
    pe.correctBoundaryConditions();

    if(runTime.outputTime())
    {
        pe.write();
    }
    Info << endl;
}