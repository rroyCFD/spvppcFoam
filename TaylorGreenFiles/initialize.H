if(TaylorGreenVortex)
{
    // Colloacated coordinate points
    volScalarField x_c = mesh.C().component(vector::X);
    volScalarField y_c = mesh.C().component(vector::Y);
    volScalarField z_c = mesh.C().component(vector::Z);

    if(dimension =="2D")
    {
         Ustart =  Uinit*( vector(1,0,0) * sin(x_c/L) * cos(y_c/L)
                         - vector(0,1,0) * cos(x_c/L) * sin(y_c/L)
                         + vector(0,0,1) * scalar(0.));

         pstart = sqr(Uinit)/4 * (cos(2.*x_c/L) + cos(2.*y_c/L));

    } else if (dimension == "3D")
    {
        Ustart = Uinit*( vector(1,0,0) * sin(x_c/L) * cos(y_c/L) * cos(z_c/L)
                       - vector(0,1,0) * cos(x_c/L) * sin(y_c/L) * cos(z_c/L)
                       + vector(0,0,1) * scalar(0.));

        pstart = sqr(Uinit)/16 * (cos(2.*x_c/L) + cos(2.*y_c/L)) * (cos(2.*y_c/L) + 2);
    }

     if(pRefOn)
     {
         Info << "Ref. cell coordinate: " << mesh.C()[pRefCell] << endl;
         pStartRefValueAna.value() = getRefCellValue(pstart, pRefCell);
         Info << "analytical pressure value at ref. cell "
              << pStartRefValueAna.value() << endl;
     }

    // Initialize velocity and pressure field at time 0
    if (runTime.timeName() == "0")
    {
        // dimensionedScalar minNu = min(laminarTransport.nu());
        if( nu.value()> 0.0 )
        {
            // Reynolds number
            dimensionedScalar Re = Uinit*L/nu;
            Info<< "\nInitializing Taylor-Greem Vortex flow with Re = "<< Re.value() << "\n" << endl;
        }

        U = Ustart;
        U.correctBoundaryConditions();
        Info<< "Writing field U" << endl;
        U.write();

        p = pstart;
        p.correctBoundaryConditions();
        Info << "pressure value at ref. cell " << getRefCellValue(p, pRefCell) << endl;
        Info<< "Writing field p" << endl;
        p.write();
    }
}
