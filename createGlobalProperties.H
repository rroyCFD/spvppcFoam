// Declare a Kinetic Energy term-by-term analysis object (if present)
//------------------------------------------------------------------------//

kineticEnergyAnalysis* KE;

// Check for kinetic energy analysis switch
bool KEAnalysis (piso.dict().lookupOrDefault<bool>("KEAnalysis", false));
Info << "kinetic energy analysis: " << KEAnalysis << endl;

if(KEAnalysis)
{
    KE = new kineticEnergyAnalysis(U, phi, p);
    KE->setPropertiesOutput();
    KE->analyzeKEBalance();
    KE->getAddtionalKETerms();
    KE->writeAvgValues();
}
else
{
    delete KE;
}

//------------------------------------------------------------------------//


// Declare a TGV object pointer and assign if TGV present
TaylorGreenVortex* TGVPtr(nullptr);
bool TGVOn (false);

fileName dictPath(runTime.constantPath());
if(Pstream::parRun())
{
    dictPath = runTime.path()+"/../constant/";
    Info << "Constant dictionary path: " << dictPath << endl;
}

// Check for Taylor Green Vortex
IOobject TGVPropertiesHeader
(
    "TaylorGreenVortexProperties",
    dictPath, // runTime.constant(),
    mesh,
    IOobject::NO_READ
);

// debug: decomposed case execution: points to processor0/constant folder
// Info << TGVPropertiesHeader.path() << endl;
// Info << TGVPropertiesHeader.rootPath() << endl;
// Info << TGVPropertiesHeader.objectPath() << endl;

if(TGVPropertiesHeader.typeHeaderOk<dictionary>(true))
{
    TGVOn = true;

    IOdictionary TGVProperties_
    (
        IOobject
        (
            "TaylorGreenVortexProperties",
            dictPath,
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const word dimension_(TGVProperties_.lookup("dimension"));
    // dimension_ = (TGVProperties_.lookupOrDefault<word>("dimension", "1D"))

    if(dimension_ == "2D")
    {
         TGVPtr = new TaylorGreenVortex2D(U, phi, p, pRefCell);
    }
    else if(dimension_ == "3D")
    {
        TGVPtr = new TaylorGreenVortex3D(U, phi, p);
    }
    else
    {
        FatalErrorInFunction
            << "Taylor-Green Vortex dimension " << dimension_ << " is invalid!"
            << "\n Options are 2D or 3D." << abort(FatalError);
    }


    // Initialize the fields if required
    Switch initialize_(false);
    initialize_.readIfPresent("initialize", TGVProperties_);
    Info << "TGV initialize: " << initialize_ << endl;

    if(initialize_ && runTime.startTime().value() == 0)
    {
        TGVPtr->setInitialFieldsAsAnalytical(U, phi, p);
        U.write();
        phi.write();
        p.write();
    }
}
else
{
    delete TGVPtr;
}
//------------------------------------------------------------------------//