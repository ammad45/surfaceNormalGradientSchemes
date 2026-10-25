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

#include "newCorrectedSnGrad.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "linear.H"
#include "fvcGrad.H"
#include "gaussGrad.H"
#include "skewCorrectionVectors.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::SurfaceField<Type>>
Foam::fv::newCorrectedSnGrad<Type>::fullGradCorrection
(
    const VolField<Type>& vf
) const
{
    const fvMesh& mesh = this->mesh();

    tmp<SurfaceField<Type>> tssf
    (
        SurfaceField<Type>::New
        (
            "snGradCorr(" + vf.name() + ')',
            mesh,
            vf.dimensions()*mesh.nonOrthDeltaCoeffs().dimensions()
        )
    );
    SurfaceField<Type>& ssf = tssf.ref();
    ssf = dimensioned<Type>(ssf.dimensions(), Zero);

    typedef typename pTraits<Type>::cmptType cmptType;
    typedef typename outerProduct<vector, cmptType>::type CmptGradType;
    typedef typename outerProduct<vector, CmptGradType>::type HessianType;

    const skewCorrectionVectors& scv =
        skewCorrectionVectors::New(mesh);

    const bool useSkew = scv.skew();

    tmp<surfaceVectorField> tSkewVec;
    if (useSkew)
    {
        const dimensionedScalar unitLength("unitLength", dimLength, 1.0);
        tSkewVec = tmp<surfaceVectorField>(unitLength*scv());
    }

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

        // Non-orthogonal correction
        tmp<SurfaceField<cmptType>> tCorr
        (
            mesh.nonOrthCorrectionVectors() & tGradF()
        );

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

            tCorr.ref() += tNhat() & (tHf() & tSkewVec());
        }

        ssf.replace(cmpt, tCorr());
    }

    return tssf;
}


template<class Type>
Foam::tmp<Foam::SurfaceField<Type>>
Foam::fv::newCorrectedSnGrad<Type>::correction
(
    const VolField<Type>& vf
) const
{
    const fvMesh& mesh = this->mesh();

    tmp<SurfaceField<Type>> tssf
    (
        SurfaceField<Type>::New
        (
            "snGradCorr(" + vf.name() + ')',
            mesh,
            vf.dimensions()*mesh.nonOrthDeltaCoeffs().dimensions()
        )
    );
    SurfaceField<Type>& ssf = tssf.ref();

    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; ++cmpt)
    {
        ssf.replace
        (
            cmpt,
            newCorrectedSnGrad<typename pTraits<Type>::cmptType>(mesh)
           .fullGradCorrection(vf.component(cmpt))
        );
    }

    return tssf;
}


// ************************************************************************* //
