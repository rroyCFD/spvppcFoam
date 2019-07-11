if(TaylorGreenVortex && runTime.timeName() == "0")
{
    // Staggered coordinate points
    surfaceScalarField x_s = mesh.Cf().component(vector::X);
    surfaceScalarField y_s = mesh.Cf().component(vector::Y);
    surfaceScalarField z_s = mesh.Cf().component(vector::Z);

    surfaceScalarField phiStart = phi;

    if(dimension =="2D")
    {
        // phiStart =  Uinit *
        //             (     vector(1,0,0) * sin(x_s/L) * cos(y_s/L)
        //                 - vector(0,1,0) * cos(x_s/L) * sin(y_s/L)
        //                 + vector(0,0,1) * scalar(0.)
        //             ) & mesh.Sf();

    }
    else if (dimension == "3D")
    {
        phiStart = Uinit *
                    (    vector(1,0,0) * sin(x_s/L) * cos(y_s/L) * cos(z_s/L)
                       - vector(0,1,0) * cos(x_s/L) * sin(y_s/L) * cos(z_s/L)
                       + vector(0,0,1) * scalar(0.)
                    ) & mesh.Sf();
    }

    if (runTime.timeName() == "0")
    {
        Info << "Writing analytical field phi" << endl;
        phi = phiStart;
        phi.write();
    }
}

