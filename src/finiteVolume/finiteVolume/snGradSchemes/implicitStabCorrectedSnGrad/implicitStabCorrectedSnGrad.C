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
#include "implicitStabCorrectedSnGrad.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "linear.H"
#include "fvcGrad.H"
#include "gaussGrad.H"
#include "fvc.H"
#include "skewCorrectionVectors.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
implicitStabCorrectedSnGrad<Type>::correction
(
    const GeometricField<Type, fvPatchField, volMesh>&
) const
{
    return nullptr;
}


template<class Type>
tmp<surfaceScalarField> implicitStabCorrectedSnGrad<Type>::deltaCoeffs
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    const fvMesh& mesh = this->mesh();

    // Start from the standard over-relaxed non-orthogonal delta.
    tmp<surfaceScalarField> tdeltaCoeffs
    (
        new surfaceScalarField
        (
            IOobject
            (
                "implicitStabCorrectedDeltaCoeffs(" + vf.name() + ')',
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh.nonOrthDeltaCoeffs()
        )
    );
    surfaceScalarField& deltaCoeffs = tdeltaCoeffs.ref();

    typedef typename pTraits<Type>::cmptType cmptType;
    typedef typename outerProduct<vector, cmptType>::type CmptGradType;
    typedef typename outerProduct<vector, CmptGradType>::type HessianType;

    // Skewness (non-conjunctionality) vectors.
    const skewCorrectionVectors& scv =
        skewCorrectionVectors::New(mesh);

    const bool useSkew = scv.skew();

    tmp<surfaceVectorField> tSkewVec;
    if (useSkew)
    {
        const dimensionedScalar unitLength("unitLength", dimLength, 1.0);
        tSkewVec = tmp<surfaceVectorField>(unitLength*scv());
    }

    // Unit face normal for the Hessian projection.
    tmp<surfaceVectorField> tNhat(mesh.Sf()/mesh.magSf());

    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; ++cmpt)
    {
        tmp<GeometricField<CmptGradType, fvPatchField, volMesh>> tGradVol
        (
            gradScheme<cmptType>::New
            (
                mesh,
                mesh.gradScheme("grad(" + vf.name() + ')')
            )().grad(vf.component(cmpt))
        );

        tmp<GeometricField<CmptGradType, fvsPatchField, surfaceMesh>> tGradF
        (
            linear<CmptGradType>(mesh).interpolate(tGradVol())
        );

        // Non-orthogonal contribution
        tmp<GeometricField<cmptType, fvsPatchField, surfaceMesh>> tCorr
        (
            mesh.nonOrthCorrectionVectors() & tGradF()
        );
        tCorr.ref().setOriented();
        GeometricField<cmptType, fvsPatchField, surfaceMesh>& corr = tCorr.ref();

        // Non-conjunctional (skew) contribution
        if (useSkew)
        {
            tmp<GeometricField<HessianType, fvPatchField, volMesh>> tHvol
            (
                fvc::grad(tGradVol())
            );

            tmp<GeometricField<HessianType, fvsPatchField, surfaceMesh>> tHf
            (
                linear<HessianType>(mesh).interpolate(tHvol())
            );

            tmp<GeometricField<cmptType, fvsPatchField, surfaceMesh>> tSkew
            (
                tNhat() & (tHf() & tSkewVec())
            );
            tSkew.ref().setOriented();

            corr += tSkew();
        }

        tGradVol.clear();
        tGradF.clear();

        // ---------------- Internal faces ----------------
        // Purely stabilising: add |corr| so the effective delta is
        // always >= the orthogonal one.
        forAll(deltaCoeffs, facei)
        {
            deltaCoeffs[facei] += mag(corr[facei]);
        }

        // ---------------- Coupled patches ----------------
        surfaceScalarField::Boundary& deltaBf =
            deltaCoeffs.boundaryFieldRef();

        forAll(deltaBf, patchi)
        {
            fvsPatchScalarField& pDelta = deltaBf[patchi];

            if (!pDelta.coupled())
            {
                continue;
            }

            const fvsPatchField<cmptType>& pCorr =
                corr.boundaryField()[patchi];

            forAll(pDelta, patchFacei)
            {
                pDelta[patchFacei] += mag(pCorr[patchFacei]);
            }
        }
    }

    return tdeltaCoeffs;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
