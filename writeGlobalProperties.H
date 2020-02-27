// kinetic energy analysis
if(KEAnalysis)
{
     Info << "Time = " << runTime.elapsedCpuTime() << " s" << endl;
     KE->analyzeKEBalance();
     KE->getAddtionalKETerms();
     KE->writeAvgValues();
}


// Taylor-Green vortex
if(TGVOn)
{
   TGVPtr->calcProperties();
   TGVPtr->writeProperties();
}