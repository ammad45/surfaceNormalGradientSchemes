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

#include "newCorrectedSnGrad.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "linear.H"
#include "fvcGrad.H"
#include "gaussGrad.H"
#include "fvc.H"
#include "skewCorrectionVectors.H"


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::fv::newCorrectedSnGrad<Type>::fullGradCorrection
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    const fvMesh& mesh = this->mesh();

    auto tssf = tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>::New
    (
        IOobject
        (
            "snGradCorr("+vf.name()+')',
            vf.instance(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        vf.dimensions()*mesh.nonOrthDeltaCoeffs().dimensions()
    );
    auto& ssf = tssf.ref();

    ssf.setOriented();
    ssf = dimensioned<Type>(ssf.dimensions(), Zero);

    typedef typename pTraits<Type>::cmptType cmptType;
    typedef typename outerProduct<vector, cmptType>::type CmptGradType;
    typedef typename outerProduct<vector, CmptGradType>::type HessianType;

    const skewCorrectionVectors& scv = skewCorrectionVectors::New(mesh);
    const bool useSkew = scv.skew();

    tmp<surfaceVectorField> tSkewVec;
    if (useSkew)
    {
        const dimensionedScalar unitLength("unitLength", dimLength, 1.0);
        tSkewVec = tmp<surfaceVectorField>(unitLength*scv());
    }

    const surfaceVectorField& nonOrthCorrVecs = mesh.nonOrthCorrectionVectors();
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

        // Non-orthogonal correction
        tmp<GeometricField<cmptType, fvsPatchField, surfaceMesh>> tCorr
        (
            nonOrthCorrVecs & tGradF()
        );
        tCorr.ref().setOriented();

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

            tCorr.ref() += tSkew();
        }

        ssf.replace(cmpt, tCorr());
    }

    return tssf;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::fv::newCorrectedSnGrad<Type>::correction
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    const fvMesh& mesh = this->mesh();

    auto tssf = tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>::New
    (
        IOobject
        (
            "snGradCorr("+vf.name()+')',
            vf.instance(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        vf.dimensions()*mesh.nonOrthDeltaCoeffs().dimensions()
    );
    auto& ssf = tssf.ref();
    ssf.setOriented();

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
