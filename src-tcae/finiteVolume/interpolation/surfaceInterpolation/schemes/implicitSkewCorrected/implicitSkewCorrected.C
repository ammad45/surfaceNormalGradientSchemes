/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Ported to OpenFOAM-dev (TCAE 26.03) from work
     \\/     M anipulation  | on ESI OpenFOAM-v2412.
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

#include "implicitSkewCorrected.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "linear.H"
#include "gaussGrad.H"
#include "gradScheme.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::surfaceScalarField>
Foam::implicitSkewCorrected<Type>::weights
(
    const VolField<Type>& vf
) const
{
    const fvMesh& mesh = this->mesh();

    // Start from the base scheme's weights.  This is what
    // fvm::div(phi, vf) would use if we did nothing.
    tmp<surfaceScalarField> tw
    (
        surfaceScalarField::New
        (
            "implicitSkewCorrectedWeights(" + vf.name() + ')',
            tScheme_().weights(vf)
        )
    );

    // Nothing to do if the mesh is not skew.
    const skewCorrectionVectors& scv =
        skewCorrectionVectors::New(mesh);

    if (!scv.skew())
    {
        return tw;
    }

    surfaceScalarField& w = tw.ref();

    // Length-dimensioned skew vector m = C_f - f'.
    const dimensionedScalar unitLength("unitLength", dimLength, 1.0);
    const surfaceVectorField m(unitLength*scv());

    typedef typename pTraits<Type>::cmptType cmptType;
    typedef typename outerProduct<vector, cmptType>::type CmptGradType;

    // Accumulate the numerator  Σ_i (m . grad vf_i) . (vf_N_i - vf_P_i)
    // and denominator Σ_i (vf_N_i - vf_P_i)^2 so that the scalar
    // implicit weight adjustment is  numer / denom  =
    //   scalar case:  m . grad phi / (phi_N - phi_P)              (paper Eq. 14)
    //   vector case:  ((m . grad U) . (U_N - U_P)) / |U_N - U_P|^2
    scalarField numer(mesh.nInternalFaces(), Zero);
    scalarField denom(mesh.nInternalFaces(), Zero);

    // Per-patch boundary accumulators for coupled patches
    // (one scalarField per patch, sized to the patch face count)
    PtrList<scalarField> numerBf(mesh.boundary().size());
    PtrList<scalarField> denomBf(mesh.boundary().size());
    forAll(mesh.boundary(), patchi)
    {
        const label pSize = mesh.boundary()[patchi].size();
        numerBf.set(patchi, new scalarField(pSize, Zero));
        denomBf.set(patchi, new scalarField(pSize, Zero));
    }

    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; ++cmpt)
    {
        // Cell-centred gradient of one component with the user-selected
        // grad scheme (paper §5.1 suggests leastSquares).
        tmp<VolField<CmptGradType>> tGradVol
        (
            fv::gradScheme<cmptType>::New
            (
                mesh,
                mesh.schemes().grad("grad(" + vf.name() + ')')
            )().grad(vf.component(cmpt))
        );

        // Interpolate the cell gradient to the face (linear).
        tmp<SurfaceField<CmptGradType>> tGradF
        (
            linear<CmptGradType>(mesh).interpolate(tGradVol())
        );

        // m . (grad vf_i)_{f'}
        const SurfaceField<cmptType> skewCorr(m & tGradF());

        tGradVol.clear();
        tGradF.clear();

        // Component field (owner/neighbour and coupled boundary values)
        tmp<VolField<cmptType>> tvfCmpt(vf.component(cmpt));
        const VolField<cmptType>& vfCmpt = tvfCmpt();

        // ---------------- Internal faces ----------------
        const labelUList& owner = mesh.owner();
        const labelUList& neighbour = mesh.neighbour();

        forAll(owner, facei)
        {
            const scalar dvf =
                vfCmpt[neighbour[facei]] - vfCmpt[owner[facei]];

            numer[facei] += skewCorr[facei]*dvf;
            denom[facei] += dvf*dvf;
        }

        // ---------------- Coupled patches ----------------
        forAll(mesh.boundary(), patchi)
        {
            const fvsPatchField<cmptType>& pSkew =
                skewCorr.boundaryField()[patchi];

            if (!pSkew.coupled())
            {
                continue;
            }

            const Field<cmptType> pvfInt
            (
                vfCmpt.boundaryField()[patchi].patchInternalField()
            );
            const Field<cmptType> pvfNbr
            (
                vfCmpt.boundaryField()[patchi].patchNeighbourField()
            );

            scalarField& pNumer = numerBf[patchi];
            scalarField& pDenom = denomBf[patchi];

            forAll(pNumer, patchFacei)
            {
                const scalar dvf = pvfNbr[patchFacei] - pvfInt[patchFacei];

                pNumer[patchFacei] += pSkew[patchFacei]*dvf;
                pDenom[patchFacei] += dvf*dvf;
            }
        }
    }

    // Combine numer / denom into the implicit weight adjustment
    // (skip faces where the field is nearly uniform: the correction
    // is meaningless there and dividing by a tiny denominator would
    // destroy the matrix).
    scalarField& wI = w.primitiveFieldRef();

    forAll(wI, facei)
    {
        if (denom[facei] > SMALL)
        {
            wI[facei] += numer[facei]/denom[facei];
        }
    }

    surfaceScalarField::Boundary& wBf = w.boundaryFieldRef();

    forAll(wBf, patchi)
    {
        fvsPatchScalarField& pw = wBf[patchi];

        if (!pw.coupled())
        {
            continue;
        }

        const scalarField& pNumer = numerBf[patchi];
        const scalarField& pDenom = denomBf[patchi];

        forAll(pw, patchFacei)
        {
            if (pDenom[patchFacei] > SMALL)
            {
                pw[patchFacei] += pNumer[patchFacei]/pDenom[patchFacei];
            }
        }
    }

    return tw;
}


// ************************************************************************* //
