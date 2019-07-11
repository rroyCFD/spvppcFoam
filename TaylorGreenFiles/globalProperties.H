if(TaylorGreenVortex)
{
    // Dissipation rate weighted average
//    epsilon = 2*(magSqr(symm(fvc::grad(U)))*laminarTransport.nu())().weightedAverage(mesh.V()) /(Uinit*Uinit*Uinit/L);

    volScalarField  dissipation = (U & fvc::laplacian(turbulence->nuEff(), U));
    epsilon = dissipation.weightedAverage(mesh.V())/(Uinit*Uinit*Uinit/L);

    // Kinetic enery weighted average
    Ek = 0.5*magSqr(U)().weightedAverage(mesh.V())/ (Uinit*Uinit);

    if(Pstream::master())
    {
        if(dimension =="2D")
        {
             globalPropertiesFile
             << runTime.value()/tc.value() << "\t" << Ek.value() << "\t" << epsilon.value()
             << "\t" << UL2err.value() << "\t" << pL2err.value()
             << "\t" << ULinfErr.value() << "\t" << pLinfErr.value()<< endl;
        }
        else if (dimension =="3D")
        {
            globalPropertiesFile << runTime.value()/tc.value() << "\t" << Ek.value() << "\t" << epsilon.value()<< endl;
        }
    }
}

