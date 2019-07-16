// create error fields
#include "vector.H"
volVectorField Ustart
(
    IOobject
    (
        "Ustart",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE // AUTO_WRITE
    ),
    // mesh,
    // dimensionedVector("Ustart",dimVelocity,vector::zero)
    U
);

volScalarField pstart
(
    IOobject
    (
        "pstart",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE // AUTO_WRITE
    ),
    // mesh,
    // dimensionedScalar("pstart",dimensionSet(0, 2, -2, 0, 0, 0, 0), 0.0)
    p
);

volVectorField Uerror
 (
     IOobject
     (
         "Uerror",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
     ),
     mesh,
     dimensionedVector("Uerror",dimVelocity,vector::zero)
 );

volScalarField perror
 (
     IOobject
     (
         "perror",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
     ),
     mesh,
     dimensionedScalar("perror",dimensionSet(0, 2, -2, 0, 0, 0, 0), 0.0)
 );


 dimensionedScalar ULinfErr = max(mag(Uerror));          // infinity norm
 dimensionedScalar UL2err   = sqrt(sum(magSqr(Uerror))); // 2nd norm

 dimensionedScalar pLinfErr = max(mag(perror));          // infinity norm
 dimensionedScalar pL2err   = sqrt(sum(magSqr(perror))); // 2nd norm

 dimensionedScalar pRefValueAna("pRefValueAna", p.dimensions(), 0);
 dimensionedScalar pStartRefValueAna("pStartRefValueAna", p.dimensions(), 0);
 dimensionedScalar pCurrentRefValue
 ("pCurrentRefValue", p.dimensions(), getRefCellValue(p, pRefCell));


