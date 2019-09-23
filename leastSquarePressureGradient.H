/*
Pressure gradient error (difference between collocated and reconstructed surface
normal) is minimized using least square approach. Hence, the collocated gradient
is derived from snGrad values.

Reference: Equation 54 at Mahesh, K., Constantinescu, G., & Moin, P. (2003).
A numerical method for large-eddy simulation in complex geometries.
Computational Fluid and Solid Mechanics 2003, 197, 23–30.
https://doi.org/10.1016/B978-008044046-0.50009-9

*/


// void leastSquareGradientField(volVectorField& pGradLS, volScalarField& p)
// {
//     surfaceScalarField pGradf = fvc::snGrad(p);

//     // loop over all faces to build the symmTensorField(LHS) and vectorField (RHS)
//     // then, loop over all cell to perform inverse and inner product to get the gradient
//     const labelUList& owner = mesh.owner();
//     const labelUList& neighbour = mesh.neighbour();

//     const surfaceVectorField& Sf = mesh.Sf();

//     vectorField RHS(mesh.nCells(), vector::zero);

//     // Create the base form of the dd-tensor
//     // including components for the "empty" directions
//     symmTensor dd0(sqr((Vector<label>::one - mesh.geometricD())/2));
//     // Info << dd0 << endl;

//     symmTensorField LHS(mesh.nCells(), dd0); //symmTensor::zero

//     forAll(owner, f)
//     {
//         vector rhs = pGradf[f]*Sf[f];
//         symmTensor lhs = symm(Sf[f] * Sf[f])/mag(Sf[f]);

//         label oci = owner[f];     // owner cell index
//         label nci = neighbour[f]; // neighbour cell index

//         RHS[oci] += rhs;
//         RHS[nci] += rhs;

//         LHS[oci] += lhs;
//         LHS[nci] += lhs;
//         // watch(f);
//     }

//     forAll(mesh.boundaryMesh(), patchI)
//     {
//         fvsPatchScalarField& pGradPatch = pGradf.boundaryFieldRef()[patchI];
//         const fvPatch& cPatch = pGradPatch.patch();

//         // if(cPatch.type() == "empty")
//         // {
//         //    continue;
//         // }
//         // watch(patchI);

//         const labelUList& faceCells = cPatch.faceCells();
//         const vectorField& APatch = cPatch.Sf();

//         forAll(pGradPatch, patchFaceI)
//         {
//             label oci = faceCells[patchFaceI]; // owner cell index
//             vector Af = APatch[patchFaceI];

//             RHS[oci] += pGradPatch[patchFaceI]*Af;

//             LHS[oci] += symm(Af*Af)/mag(Af);
//         }
//     }

//     forAll(pGradLS, c)
//     {
//         // inverse and substract the empty direction tensor
//         // LHS[c] = inv(LHS[c]) - dd0;
//         pGradLS[c] = (inv(LHS[c]) - dd0) & RHS[c];
//     }

//     return;
// }



volSymmTensorField leastSquareGradientCoeffs(const fvMesh& mesh, const Time& runTime)
{
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    const surfaceVectorField& Sf = mesh.Sf();

    // Create the base tensor including components for the "empty" directions
    symmTensor dd0(sqr((Vector<label>::one - mesh.geometricD())/2));
    // symmTensorField LHS(mesh.nCells(), dd0); //symmTensor::zero

    volSymmTensorField LHS
    (
        IOobject
        (
            "LHS",
            runTime.constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("",dimless/dimArea,dd0)
    );

    forAll(owner, faceI)
    {
        symmTensor lhs = symm(Sf[faceI] * Sf[faceI])/mag(Sf[faceI]);

        label oci = owner[faceI];     // owner cell index
        label nci = neighbour[faceI]; // neighbour cell index

        LHS[oci] += lhs;
        LHS[nci] += lhs;
    }

    forAll(mesh.boundary(), patchI)
    {
        const polyPatch& cPatch = mesh.boundary()[patchI].patch();

        const labelUList& faceCells = cPatch.faceCells();
        const vectorField& APatch = cPatch.faceAreas();

        forAll(APatch, patchFaceI)
        {
            vector Af = APatch[patchFaceI];
            LHS[faceCells[patchFaceI]] += symm(Af*Af)/mag(Af);
        }
    }

    forAll(LHS, c)
    {
        // inverse and substract the empty direction tensor
        LHS[c] = inv(LHS[c]) - dd0;
    }

    return LHS;
}


tmp<volVectorField> getLSGrad(volSymmTensorField& LHS, volScalarField& p)
{
    const fvMesh& mesh = p.mesh();

    //tmp<volVectorField> gradLS = LHS & fvc::surfaceSum(fvc::snGrad(p) * mesh.Sf());
    tmp<volVectorField> gradLS
    (
        new volVectorField
        (
            IOobject
            (
                p.name()+"GradLS",
                p.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            //mesh,
            //dimensionedVector("",p.dimesions()/dimLength,vector::zero)
            LHS & fvc::surfaceSum(fvc::snGrad(p) * mesh.Sf())
        )
    );

    return gradLS;
}


// Determine the least square error in pressure gradient calculation
// Find error between least square gradient and reconstructed surfaceNormal grad
// output of equation 54 on above refered publication
tmp<volScalarField> getLSError(tmp<volVectorField> gradLS, volScalarField& p)
{
    const fvMesh& mesh = p.mesh();

    tmp<volScalarField> terrorGrad
    (
        new volScalarField
            (
                IOobject
                (
                    "errorGrad",
                    p.instance(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("",p.dimensions()/dimLength,Zero)
            )
    );
    volScalarField& errorGrad = terrorGrad.ref();


    // loop over each cell to find the error
    volVectorField& vvf = gradLS.ref();
    surfaceScalarField snGradf = fvc::snGrad(p)*mesh.magSf();

    forAll(vvf, cellI)
    {
        // list of all faces in a cell
        const labelList& faceList =  mesh.cells()[cellI];

        // loop over all faces in each cell
        forAll(faceList, faceI)
        {
            errorGrad[cellI] += ((vvf[cellI] & mesh.Sf()[faceI]) - snGradf[faceI]);
        }
    }

    return terrorGrad;
}

