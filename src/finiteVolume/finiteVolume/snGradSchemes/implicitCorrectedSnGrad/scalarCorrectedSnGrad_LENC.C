/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "snGradScheme.H"
#include "scalarCorrectedSnGrad.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "linear.H"
#include "fvcGrad.H"
#include "gaussGrad.H"
#include "fvc.H"
#include "skewCorrectionVectors.H"


// * * * * * * * * * * * * * * *  * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
scalarCorrectedSnGrad<Type>::correction
(
    const GeometricField<Type, fvPatchField, volMesh>&
) const
{
    NotImplemented;
    return nullptr;
    //return tnonOrthCorr_();
}




template<class Type>
tmp<surfaceScalarField> scalarCorrectedSnGrad<Type>::deltaCoeffs
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    //Info << "LENC-scheme used on: " << vf.name() << endl;
    // Implementation for implicitly limited 'Over-relaxed approach'
    const fvMesh& mesh = this->mesh();
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    tmp<surfaceScalarField> tdeltaCoeffs 
    (
        new surfaceScalarField
        (
            IOobject
            (
                "deltaCoeffs",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            mesh,
            dimensionedScalar("deltaCoeffs", dimless/dimLength, 0.0)
        )
    );
    surfaceScalarField& deltaCoeffs = tdeltaCoeffs.ref();
    
    // Initialize with mesh deltaCoeffs
    deltaCoeffs = mesh.nonOrthDeltaCoeffs();
    
    GeometricField<Type, fvsPatchField, surfaceMesh>& nonOrthCorr 
        = tnonOrthCorr_.ref();
    nonOrthCorr *= dimensionedScalar("unitVF", vf.dimensions(), 1.);

    // Hard-coded gradient calculation of 'vf' with 'leastSquares'-scheme
    IStringStream gradSchemeLS("leastSquares");
    tmp<GeometricField<typename outerProduct<vector, Type>::type, 
        fvPatchField, volMesh> > tgradc = fv::gradScheme<Type>::New
        (
            mesh, 
            gradSchemeLS
        )().grad(vf, "grad(" + vf.name() + ')');
    GeometricField<typename outerProduct<vector, Type>::type, 
        fvPatchField, volMesh>& gradc 
        = tgradc.ref();
    
    nonOrthCorr =
        mesh.nonOrthCorrectionVectors()
      & linear<typename outerProduct<vector, Type>::type>(mesh).interpolate
        (gradc);

    // Skewness Correction
    if (true)
    {
        // gamma is needed because skewCorrectionVectors should have the unit [m]
        // but it is dimless.
        // The unit of mesh.Sf() as well as mesh.magSf() is [m^2]
        dimensioned<scalar> unit("unit", dimensionSet(0, 1, 0, 0, 0, 0, 0), 1);
        IStringStream gradSchemeLS2("leastSquares");
        nonOrthCorr +=
            (mesh.Sf()/mesh.magSf()) & 
            (
                 unit*skewCorrectionVectors::New(mesh)()
               & linear<typename outerProduct<vector, vector>::type>(mesh).interpolate
                 (
                     fv::gradScheme<vector>::New
                     (
                         mesh, 
                         gradSchemeLS2
                     )().grad(gradc, "grad(grad(" + vf.name() + "))")
                 )
            );
    }

    // Convert the templated correction to scalar values and add to deltaCoeffs
    forAll(deltaCoeffs, faceI)
    {
        deltaCoeffs[faceI] += mag(nonOrthCorr[faceI]);
    }

    return tdeltaCoeffs;
}

}

}
