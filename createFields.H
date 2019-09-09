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


//volScalarField LESdelta
//    (
//        IOobject
//        (
//            "LESdelta",
//            runTime.timeName(),
//            mesh,
//            IOobject::NO_READ,
//            IOobject::AUTO_WRITE
//        ),
//        mesh,
//        dimensionedScalar("LESdelta", dimLength, SMALL)
//    );

//if (mesh.objectRegistry::foundObject<volScalarField>("delta"))
//{
//    const volScalarField& delta =
//            mesh.objectRegistry::lookupObject<volScalarField>("delta");

//    LESdelta.primitiveFieldRef() = delta.primitiveField();
//}
