// filter the extrapolated pressure
// p^star = filter(pe) + filter(pe - filter(pe))
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


    tmp<volScalarField> tpePrime
    (
        new volScalarField
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
        )
     );
    volScalarField pePrime = tpePrime.ref();
    pePrime.correctBoundaryConditions();

    // 4th order filter; add filter of the (pe-pef) field
    pef += C4Regularization.polyLaplaceFilter(pePrime);
    pef.correctBoundaryConditions();

    tpePrime.clear();

    // bound the filtered pressure field
    dimensionedScalar peMax ("peMax", pe.dimensions(), gMax(pe));
    dimensionedScalar peMin ("peMin", pe.dimensions(), gMin(pe));
    Info << peMax << nl << peMin << endl;
    pef = max(min(pef, peMax), peMin);

    // write filtered and residual field
    if(runTime.outputTime())
    {
        pef.write();

        // residual field
        volScalarField pefPrime
        (
            IOobject
            (
                "pefPrime",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            (pe - pef),
            p.boundaryField().types()
        );
        pefPrime.correctBoundaryConditions();
        pefPrime.write();

        dimensionedScalar pefPrimeMax
            ("pefPrimeMax", pefPrime.dimensions(), gMax(pefPrime));
        dimensionedScalar pefPrimeMin
            ("pefPrimeMin", pefPrime.dimensions(), gMin(pefPrime));
        Info << pefPrimeMax << nl << pefPrimeMin << endl;
    }

    // assign the extrapolate-filtered field to the extrapolated field
    pe = pef;
    pe.correctBoundaryConditions();

    Info << endl;
}