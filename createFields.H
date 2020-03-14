// Read fields: p, pp, U and phi(if present) and their old time stepss
Info<< "Reading field p" << endl;
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

label pRefCell = 0;
scalar pRefValue = 0.0;

const bool pRefOn (setRefCell(p, piso.dict(), pRefCell, pRefValue));
if(pRefOn)
{

    Info << "pRefCell: " << pRefCell << "\tpRefValue: " << pRefValue << tab
         << "initial field ref value: " << getRefCellValue(p, pRefCell) << endl;

    if( mag(getRefCellValue(p, pRefCell) - pRefValue) > SMALL)
    {
        dimensionedScalar pCorr
        (
            "pCorr",
            p.dimensions(),
            (getRefCellValue(p, pRefCell) - pRefValue)
        );
        p -= pCorr;
        p.correctBoundaryConditions();

        Info << "pRefCell: " << pRefCell << "\tpRefValue: " << pRefValue << tab
        << "corrected field ref value: " <<  getRefCellValue(p, pRefCell) << endl;
    }
}

Info<< "Reading field pp" << endl;
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
mesh.setFluxRequired(pp.name());

Info<< "\nReading field U" << endl;
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

Info    << "\tU.nOldTimes(): "   << U.nOldTimes()   << nl
        << "\tphi.nOldTimes(): " << phi.nOldTimes() << nl
        << "\tp.nOldTimes(): "   << p.nOldTimes()   <<nl << endl;

// add transport and turbulence model
singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

#include "createMRF.H"
#include "createFvOptions.H"

// Check for the RhieChow correction
bool RhieChow (piso.dict().lookupOrDefault<bool>("RhieChow", true));
Info << "RhieChow: " << RhieChow << endl;

