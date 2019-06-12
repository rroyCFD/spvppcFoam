// Coefficients for filtering
dimensionedScalar d1
(
    "d1",
    dimless,
    -(3.0/8.0)    // epsilon = 3
    // -(1.0/6.0) // epsilon = 2
);

dimensionedScalar d2
(
    "d2",
    dimless,
    (3.0/80.0)     // epsilon = 3
    // (1.0/240.0) // epsilon = 2
);


volVectorField convection
(
    IOobject
    (
        "convection",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector ("zero",dimAcceleration,vector::zero)
);
//----------------------------------------------------------------------------//