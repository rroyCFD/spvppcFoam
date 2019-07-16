if(TaylorGreenVortex && dimension =="2D")
{
    // Velocity error norms--------------------------------------------------------------//
    Info << "Calculating Uerror" << endl;
    Uerror = Ustart*Foam::exp(-2.0*nu.value()*runTime.value()) - U;

    ULinfErr.value() = gMax(mag(Uerror)());   // infinity norm
    UL2err = sqrt(magSqr(Uerror)().weightedAverage(mesh.V())); // 2nd norm

    // Pressure error norms ----------------------------------------------------------//
    Info << "Calculating perror" << endl;

    if(pRefOn) {
        // Reference analytical pressure at pRefCell
        pRefValueAna  = pStartRefValueAna * Foam::exp(-4.*nu.value()*runTime.value());
        pCurrentRefValue.value() = getRefCellValue(p, pRefCell);

        Info << "Value at pRefCell, Analytical :"<< pRefValueAna.value() << tab
             << "Current: " << pCurrentRefValue.value() << endl;

        perror = pstart*Foam::exp(-4.0*nu.value()*runTime.value()) - p
                - pRefValueAna + pCurrentRefValue;
    }
    else
    {
        perror = pstart*Foam::exp(-4.0*nu.value()*runTime.value()) - p;
    }

    pLinfErr.value() = gMax(mag(perror)()); // infinity norm
    pL2err = sqrt(magSqr(perror)().weightedAverage(mesh.V()));  // 2nd norm

    Info << "Linf error norm | Velocity: " << ULinfErr.value()
         << "\tPressure: " << pLinfErr.value() << endl;
    Info << "L2   error norm | Velocity: " << UL2err.value()
         << "\tPressure: " << pL2err.value()   << endl;


    Info << "Pressure decay coeffs: " << setprecision(10)
         << Foam::exp(-4.0*nu.value()*runTime.value()) << endl;
    Info << "Velocity decay coeffs: " << setprecision(10)
         << Foam::exp(-2.0*nu.value()*runTime.value()) << endl;

    // write the error induced during startTime and output time
    if( (runTime.value() == runTime.startTime().value()) || runTime.outputTime())
    {
        Uerror.write();
        perror.write();
    }
}
