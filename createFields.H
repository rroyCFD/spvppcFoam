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
setRefCell(p, piso.dict(), pRefCell, pRefValue);

const bool pRefOn (p.needReference());
if(pRefOn)
{
    Info << "pRefCell: " << pRefCell << "\tpRefValue: " << pRefValue << tab
         << "initial field ref value: " << getRefCellValue(p, pRefCell) << endl;
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

// Check the convection-term-regularization on or off (projection)
bool regOn (piso.dict().lookupOrDefault<bool>("regularizationOn", true));
Info << "regularizationOn: " << regOn << endl;

// Check for the RhieChow correction
bool RhieChow (piso.dict().lookupOrDefault<bool>("RhieChow", true));
Info << "RhieChow: " << RhieChow << endl;

// Check the whether to set extrapolated-filtered fields to make divergence free
bool filterFieldDivFree(true);
if(regOn)
{
    //filterFieldDivFree = readBool(piso.dict().lookup("filterFieldDivFree"));
    filterFieldDivFree =
            piso.dict().lookupOrDefault<bool>("filterFieldDivFree", true);
    Info << "Extrapolated-filtered fields made divergence free: "
         << filterFieldDivFree << endl;
}

// Check for kinetic energy analysis switch
bool KEAnalysis (piso.dict().lookupOrDefault<bool>("KEAnalysis", false));
Info << "kinetic energy analysis: " << KEAnalysis << endl;
