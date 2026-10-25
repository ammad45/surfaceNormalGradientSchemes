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
#include "implicitCorrectedSnGrad.H"
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

namespace
{
    //- Add a signed deltaCoeff correction without crossing zero.
    inline scalar limitDeltaIncrement(const scalar delta, const scalar dC)
    {
        return max(dC, -delta);
    }


    //- Keep laplacian face coefficients strictly positive.
    //  If the implicit correction would (nearly) cancel the orthogonal
    //  coefficient, revert to the uncorrected value on that face only.
    inline scalar stabiliseDeltaCoeff(const scalar delta, const scalar deltaOrth)
    {
        if (deltaOrth > ROOTVSMALL && delta < 0.01*deltaOrth)
        {
            return deltaOrth;
        }

        return max(delta, scalar(0));
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
implicitCorrectedSnGrad<Type>::correction
(
    const GeometricField<Type, fvPatchField, volMesh>&
) const
{
    return nullptr;
}


template<class Type>
tmp<surfaceScalarField> implicitCorrectedSnGrad<Type>::deltaCoeffs
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
                "implicitCorrectedDeltaCoeffs(" + vf.name() + ')',
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

    // Skew (non-conjunctionality) correction, evaluated once for the
    // mesh.  Note: skewCorrectionVectors_ carries dimless dimensions
    // even though the numerical values are in metres; multiply by a
    // unit length dimensioned scalar to restore length dimensions.
    const skewCorrectionVectors& scv =
        skewCorrectionVectors::New(mesh);

    const bool useSkew = scv.skew();

    tmp<surfaceVectorField> tSkewVec;
    if (useSkew)
    {
        const dimensionedScalar unitLength("unitLength", dimLength, 1.0);
        tSkewVec = tmp<surfaceVectorField>(unitLength*scv());
    }

    // Unit face normal, used to project the Hessian correction back to
    // a surface-normal gradient.
    tmp<surfaceVectorField> tNhat(mesh.Sf()/mesh.magSf());

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

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

        // Non-orthogonal contribution:  (k/|Sf|) . (grad phi)_f
        tmp<GeometricField<cmptType, fvsPatchField, surfaceMesh>> tCorr
        (
            mesh.nonOrthCorrectionVectors() & tGradF()
        );
        tCorr.ref().setOriented();
        GeometricField<cmptType, fvsPatchField, surfaceMesh>& corr = tCorr.ref();

        // Add non-conjunctional (skew) contribution:
        //   nHat . ( Hessian_f . m )
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

        // Keep the tmp alive until we're done with vfCmpt.
        tmp<GeometricField<cmptType, fvPatchField, volMesh>> tvfCmpt
        (
            vf.component(cmpt)
        );
        const GeometricField<cmptType, fvPatchField, volMesh>& vfCmpt =
            tvfCmpt();

        // ---------------- Internal faces ----------------
        // Skip faces with |phi_N - phi_P| ~ 0.  Limit each component
        // increment so deltaCoeffs never goes negative; a final pass
        // below reverts faces where the correction would invalidate
        // the laplacian diagonal (GAMG divide-by-zero on D).
        forAll(owner, facei)
        {
            const scalar dPhi =
                vfCmpt[neighbour[facei]] - vfCmpt[owner[facei]];

            if (mag(dPhi) > SMALL)
            {
                deltaCoeffs[facei] += limitDeltaIncrement
                (
                    deltaCoeffs[facei],
                    corr[facei]/dPhi
                );
            }
        }

        // ---------------- Coupled patches ----------------
        // For cyclic / processor patches the neighbour value has to be
        // fetched via patchNeighbourField().  Non-coupled patches don't
        // need any correction because nonOrthCorrectionVectors is zero
        // there.
        surfaceScalarField::Boundary& deltaBf =
            deltaCoeffs.boundaryFieldRef();

        forAll(deltaBf, patchi)
        {
            fvsPatchScalarField& pDelta = deltaBf[patchi];

            if (!pDelta.coupled())
            {
                continue;
            }

            const fvPatchField<cmptType>& pvfCmpt =
                vfCmpt.boundaryField()[patchi];

            const Field<cmptType> pvfInt(pvfCmpt.patchInternalField());
            const Field<cmptType> pvfNbr(pvfCmpt.patchNeighbourField());
            const fvsPatchField<cmptType>& pCorr =
                corr.boundaryField()[patchi];

            forAll(pDelta, patchFacei)
            {
                const scalar dPhi = pvfNbr[patchFacei] - pvfInt[patchFacei];

                if (mag(dPhi) > SMALL)
                {
                    pDelta[patchFacei] += limitDeltaIncrement
                    (
                        pDelta[patchFacei],
                        pCorr[patchFacei]/dPhi
                    );
                }
            }
        }
    }

    // Final stabilisation on internal faces and coupled patches.
    const surfaceScalarField& deltaOrth = mesh.nonOrthDeltaCoeffs();

    forAll(deltaCoeffs, facei)
    {
        deltaCoeffs[facei] = stabiliseDeltaCoeff
        (
            deltaCoeffs[facei],
            deltaOrth[facei]
        );
    }

    surfaceScalarField::Boundary& deltaBf =
        deltaCoeffs.boundaryFieldRef();

    forAll(deltaBf, patchi)
    {
        fvsPatchScalarField& pDelta = deltaBf[patchi];
        const fvsPatchScalarField& pOrth =
            deltaOrth.boundaryField()[patchi];

        forAll(pDelta, patchFacei)
        {
            pDelta[patchFacei] = stabiliseDeltaCoeff
            (
                pDelta[patchFacei],
                pOrth[patchFacei]
            );
        }
    }

    return tdeltaCoeffs;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
