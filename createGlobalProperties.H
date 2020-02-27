// Declare a Kinetic Energy term-by-term analysis object (if present)
//------------------------------------------------------------------------//

kineticEnergyAnalysis* KE;
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
TaylorGreenVortex* TGVPtr;
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
    TGVPtr = new TaylorGreenVortex(U, phi, p, pRefCell);

    if(runTime.startTime().value() == 0)
    {
        TGVPtr->setInitialFieldsAsAnalytical();
    }

    TGVPtr->setupProperties();

    TGVOn = true;
    Info << "TGV Object Constructor" << endl;
}
else
{
    delete TGVPtr;
}
//------------------------------------------------------------------------//