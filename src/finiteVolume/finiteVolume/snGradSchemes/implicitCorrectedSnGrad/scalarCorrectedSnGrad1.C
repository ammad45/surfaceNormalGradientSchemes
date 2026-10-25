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
}




template<class Type>
tmp<surfaceScalarField> scalarCorrectedSnGrad<Type>::deltaCoeffs
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
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
                vf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false // Do not register
            ),
            mesh,
            dimensionedScalar("deltaCoeffs", dimless/dimLength, 0.0)
        )
    );
    surfaceScalarField& deltaCoeffs = tdeltaCoeffs.ref();
    
    // Initialize with nonOrthDeltaCoeffs
    deltaCoeffs = mesh.nonOrthDeltaCoeffs();
    
    // non-orthogonal correction part
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tnonOrthCorr =
        mesh.nonOrthCorrectionVectors()
      & linear<typename outerProduct<vector, Type>::type>(mesh).interpolate
        (
            fv::gradScheme<Type>::New
            (
                mesh,
                mesh.gradScheme("grad(" + vf.name() + ')')
            )().grad(vf, "grad(" + vf.name() + ')')
        );
    
    const GeometricField<Type, fvsPatchField, surfaceMesh>& nonOrthCorr = tnonOrthCorr();

    forAll(owner, facei)
    {
        Type delta = vf[neighbour[facei]] - vf[owner[facei]];
        scalar magDelta = mag(delta);
        if (magDelta > SMALL)
        {
            scalar correction = mag(nonOrthCorr[facei]) / magDelta;
            deltaCoeffs[facei] += correction;
        }
        else
        {
            // If delta is < SMALL, add a large correction
            deltaCoeffs[facei] += GREAT;
        }
    }

    // Limiting - ensure coefficients are non-negative
    forAll(owner, facei)
    {
        deltaCoeffs[facei] = max(deltaCoeffs[facei], 0.0);
    }
    
    return tdeltaCoeffs;
}

}

}

