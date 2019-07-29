// // tmp<surfaceScalarField> pfSnGradAf = fvc::snGrad(p)* mesh.magSf();
// volVectorField pfGrad("pfGrad",fvc::reconstruct(fvc::snGrad(p)* mesh.magSf()));
// pfGrad.write(runTime.outputTime());

// volVectorField pGradDiff("pGradDiff", (fvc::grad(p)-pfGrad));
// pGradDiff.write(runTime.outputTime());

// volScalarField pGradDiffKE ("pGradDiffKE", U & pGradDiff);
// pGradDiffKE.write(runTime.outputTime());

// scalar sumPGradDiffKE = gSum(pGradDiffKE.primitiveFieldRef() * mesh.V());
// Info << "sum of pGradDiff KE: " << sumPGradDiffKE << endl;



volVectorField ppGradDiff(
   "ppGradDiff",
   (fvc::grad(pp)-fvc::reconstruct(fvc::snGrad(pp)* mesh.magSf()))
);
// ppGradDiff.write(runTime.outputTime());

volScalarField ppGradDiffKE ("ppGradDiffKE", U & ppGradDiff);
// volScalarField ppGradDiffKE ("ppGradDiffKE", pp * mag(ppGradDiff));
ppGradDiffKE.write(runTime.outputTime());

scalar sumPPGradDiffKE = gSum(ppGradDiffKE.primitiveFieldRef() * mesh.V());
Info << "sum of ppGradDiff KE: " << sumPPGradDiffKE << endl;

