IOdictionary tgvProperties
(
    IOobject
    (
        "TaylorGreenVortexProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

const bool TaylorGreenVortex(true);

// Read the dimension of the case
const word dimension
(
    tgvProperties.lookupOrDefault<word>("dimension", "3D")
);

//  Read initial Velocity V_0
dimensionedScalar Uinit
(
    "Uinit",
    dimVelocity,
    tgvProperties
);

//  Read length scale
dimensionedScalar L
(
    "L",
    dimLength,
    tgvProperties
);

// characteristic time
dimensionedScalar tc
(
    "tc",
    dimensionSet(0, 0, 1, 0, 0),
    0.0
);
tc = L/Uinit;

if(TaylorGreenVortex)
{
    Info << "TaylorGreenVortex: "   << tab
        << "dimension: "            << dimension << tab
        << "velocity Scale: "       << Uinit.value() << tab
        << "length Scale: "         << L.value() << nl << endl;
}

// Declare variables -------------------------------------------- //
dimensionedScalar epsilon
(
    "epsilon",
    dimensionSet(0, 0, 0, 0, 0),
    0.0
);

dimensionedScalar Ek
(
    "Ek",
    dimensionSet(0, 0, 0, 0, 0),
    0.0
);


// create output file
fileName outFilePath = runTime.path()/("TaylorGreenVortexProperties.dat_"+runTime.timeName());

if(Pstream::parRun() && Pstream::master())
{
    outFilePath = runTime.path()/"../"/("TaylorGreenVortexProperties.dat_"+runTime.timeName());
}

OFstream globalPropertiesFile(outFilePath);

// write to file
if(TaylorGreenVortex && Pstream::master())
{
   if(dimension =="2D" )
    {
        globalPropertiesFile << "time" << "\t" << "Ek" << "\t" << "epsilon\t\t"
        << "L2(U)\t" << "L2(p)\t" << "Linf(U)\t" << "Linf(p)" << endl;
    }
    else if (dimension =="3D" && Pstream::master())
    {
        globalPropertiesFile << "time" << "\t" << "Ek" << "\t" << "epsilon" << endl;
    }
    else {
        globalPropertiesFile << "TG Dimension mismatch: check for dimension in the transportproperties file" << endl;
    }
} else
{
    rm(outFilePath);
}

volScalarField nuField = laminarTransport.nu();

dimensionedScalar nu
(
    "nu",
    dimensionSet(0, 2, -1, 0, 0, 0, 0),
    nuField[0]
);
