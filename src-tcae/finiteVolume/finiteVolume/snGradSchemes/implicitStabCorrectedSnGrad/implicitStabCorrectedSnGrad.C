/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Ported to OpenFOAM-dev (TCAE 26.03) from the
     \\/     M anipulation  | ESI OpenFOAM-v2412 implementation.
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

#include "implicitStabCorrectedSnGrad.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "linear.H"
#include "fvcGrad.H"
#include "gaussGrad.H"
#include "skewCorrectionVectors.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<surfaceScalarField>
implicitStabCorrectedSnGrad<Type>::deltaCoeffs
(
    const VolField<Type>& vf
) const
{
    const fvMesh& mesh = this->mesh();

    // Start from the standard over-relaxed non-orthogonal delta.
    tmp<surfaceScalarField> tdeltaCoeffs
    (
        surfaceScalarField::New
        (
            "implicitStabCorrectedDeltaCoeffs(" + vf.name() + ')',
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
        tmp<VolField<CmptGradType>> tGradVol
        (
            gradScheme<cmptType>::New
            (
                mesh,
                mesh.schemes().grad("grad(" + vf.name() + ')')
            )().grad(vf.component(cmpt))
        );

        tmp<SurfaceField<CmptGradType>> tGradF
        (
            linear<CmptGradType>(mesh).interpolate(tGradVol())
        );

        // Non-orthogonal contribution
        tmp<SurfaceField<cmptType>> tCorr
        (
            mesh.nonOrthCorrectionVectors() & tGradF()
        );
        SurfaceField<cmptType>& corr = tCorr.ref();

        // Non-conjunctional (skew) contribution
        if (useSkew)
        {
            tmp<VolField<HessianType>> tHvol
            (
                fvc::grad(tGradVol())
            );

            tmp<SurfaceField<HessianType>> tHf
            (
                linear<HessianType>(mesh).interpolate(tHvol())
            );

            corr += tNhat() & (tHf() & tSkewVec());
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
