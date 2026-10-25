/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
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

Description
    Simple central-difference snGrad scheme with implicit bounded
    non-orthogonal correction

\*---------------------------------------------------------------------------*/

#include "SINCSnGrad.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
SINCSnGrad<Type>::~SINCSnGrad()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
SINCSnGrad<Type>::correction
(
    const GeometricField<Type, fvPatchField, volMesh>&
) const
{
    notImplemented
    (
        "SINCSnGrad<Type>::correction"
        "(const GeometricField<Type, fvPatchField, volMesh>&)"
    );
    return tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >();
}


template<class Type>
tmp<surfaceScalarField> SINCSnGrad<Type>::deltaCoeffs
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    // Implementation for implicitly limited 'Over-relaxed approach'

    const fvMesh& mesh = this->mesh();
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
            
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tdeltaCoeffs
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "deltaCoeffs",
                mesh.pointsInstance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false // Do not register
            ),
            mesh.nonOrthDeltaCoeffs()
        ) 
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& deltaCoeffs 
        = tdeltaCoeffs.ref();

    // non-orthogonal correction part
    GeometricField<Type, fvsPatchField, surfaceMesh> nonOrthCorr =
        mesh.nonOrthCorrectionVectors()
      & linear<typename outerProduct<vector, Type>::type>(mesh).interpolate
        (
            fv::gradScheme<Type>::New
            (
                mesh,
                mesh.gradScheme("grad(" + vf.name() + ')')
            )().grad(vf, "grad(" + vf.name() + ')')
        );
    
    forAll(owner, facei)
    {
        scalar delta = vf[neighbour[facei]] - vf[owner[facei]];
        if (mag(delta) > SMALL)
        {
            nonOrthCorr[facei] /= delta;
        }
        else
        {
            // If delta is < SMALL nonOrthCorr has to be GREAT.
            // The sign has to be calculated.
            nonOrthCorr[facei] = 
                GREAT*(pos(delta*nonOrthCorr[facei])*2 - 1);
        }
    }

    // Add orthogonal correction part
    deltaCoeffs += nonOrthCorr;
    
    // Limiting
    forAll(owner, facei)
    {
        deltaCoeffs[facei] = max(deltaCoeffs[facei], 0);
    }
    return tdeltaCoeffs;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
