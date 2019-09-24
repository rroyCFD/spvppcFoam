// Least square gradient and error fields
if(runTime.write())
{
    Info << "Writing pressure gradients and coooresponding errors" << endl;

    // gradient estimated by least-square minimization
    tmp<volVectorField> pGradLS = getLSGrad(LHS, p);
    pGradLS.ref().write();

    // error incurred in least-square minimization method
    tmp<volScalarField> errorGrad = getLSError(pGradLS, p);
    errorGrad.ref().write();

    // gradient estimated by Gauss midPoint scheme
    volVectorField pGrad("pGrad", fvc::grad(p));
    pGrad.write();

    // difference between Gauss-midpoint and least-Square minimization
    volVectorField pGradDiff ("pGradDiff", (pGrad - pGradLS));
    pGradDiff.write();

    // gradient by snGrad reconstruction
    volVectorField pSnGrad
        ("pSnGrad", fvc::reconstruct(fvc::snGrad(p)* mesh.magSf()));
    pSnGrad.write();
}