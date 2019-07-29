Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field pp\n" << endl;
volScalarField pp
(
    IOobject
    (
        "pp",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, piso.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(pp.name());

const bool pRefOn (p.needReference());
if(pRefOn)
{
    Info << "pRefCell: " << pRefCell << "\tpRefValue: " << pRefValue << tab
         << "initial field ref value: " << getRefCellValue(p, pRefCell) << endl;
}

// add transport and turbulence model
singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

#include "createMRF.H"
#include "createFvOptions.H"

// Check the convection-term-regularization on or off (projection)
bool regOn (piso.dict().lookupOrDefault<bool>("regularizationOn", true));
Info << "regularizationOn: " << regOn << endl;


// least square estimation of pressure gradient
volSymmTensorField LHS = leastSquareGradientCoeffs(mesh, runTime);
LHS.write();
