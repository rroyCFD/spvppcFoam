// tmp<surfaceScalarField> pfSnGradAf = fvc::snGrad(p)* mesh.magSf();
volVectorField pfGrad("pfGrad",fvc::reconstruct(fvc::snGrad(p)* mesh.magSf()));
pfGrad.write(runTime.outputTime());

volVectorField pGradDiff("pGradDiff", (fvc::grad(p)-pfGrad));
pGradDiff.write(runTime.outputTime());

volScalarField pGradDiffKE ("pGradDiffKE", U & pGradDiff);
pGradDiffKE.write(runTime.outputTime());

scalar sumPGradDiffKE = gSum(pGradDiffKE.primitiveFieldRef() * mesh.V());
Info << "sum of pGradDiff KE: " << sumPGradDiffKE << endl;




//volVectorField ppGradDiff(
//    "ppGradDiff",
//    (fvc::grad(pp)-fvc::reconstruct(fvc::snGrad(pp)* mesh.magSf()))
//);
//ppGradDiff.write();

//volScalarField ppGradDiffKE ("ppGradDiffKE", U & pGradDiff);
//ppGradDiffKE.write();

//scalar sumPPGradDiffKE = gSum(pGradDiffKE.primitiveFieldRef() * mesh.V());
//Info << "sum of pGradDiff KE: " << sumPGradDiffKE << endl;

