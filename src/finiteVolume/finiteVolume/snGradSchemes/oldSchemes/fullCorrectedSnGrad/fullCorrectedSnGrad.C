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

#include "fullCorrectedSnGrad.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "linear.H"
#include "fvcGrad.H"
#include "gaussGrad.H"
#include "skewCorrectionVectors.H"


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::fv::fullCorrectedSnGrad<Type>::fullGradCorrection
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

    typedef typename outerProduct<vector, typename pTraits<Type>::cmptType>::type CmptGradType;

    const skewCorrectionVectors& scv = skewCorrectionVectors::New(mesh);
    tmp<surfaceVectorField> tScv(scv());
    const surfaceVectorField& nonOrthCorrVecs = mesh.nonOrthCorrectionVectors();

    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; ++cmpt)
    {
        //I have calculated the gradients components-wise to also include vectors and tensors
        tmp<GeometricField<CmptGradType, fvPatchField, volMesh>> tGradVol
        (
            gradScheme<typename pTraits<Type>::cmptType>::New
            (
                mesh,
                mesh.gradScheme("grad(" + vf.name() + ')')
            )().grad(vf.component(cmpt))
        );

        
        tmp<GeometricField<CmptGradType, fvsPatchField, surfaceMesh>> tGradF
        (
            linear<CmptGradType>(mesh).interpolate(tGradVol())
        );

        //This is what normally is calculated in correctedSnGrad.C  
        tmp<GeometricField<typename pTraits<Type>::cmptType, fvsPatchField, surfaceMesh>> tNonOrth
        (
            nonOrthCorrVecs & tGradF()
        );

        //This is the additional skew correction that is not normally calculated in correctedSnGrad.C
        //This still needs to refined as I have to find a way to calculate grad of grad (Future work)
        tmp<GeometricField<CmptGradType, fvPatchField, volMesh>> tGradVolGauss
        (
            fv::gaussGrad<typename pTraits<Type>::cmptType>(mesh).grad(vf.component(cmpt))
        );
        tmp<GeometricField<CmptGradType, fvsPatchField, surfaceMesh>> tGradFGauss
        (
            linear<CmptGradType>(mesh).interpolate(tGradVolGauss())
        );
        tmp<GeometricField<typename pTraits<Type>::cmptType, fvsPatchField, surfaceMesh>> tSkew
        (
            tScv() & tGradFGauss()
        );

        // Ensure consistent orientation before combining
        tNonOrth.ref().setOriented();
        tSkew.ref().setOriented();

        ssf.replace(cmpt, tNonOrth() + tSkew());
    }

    return tssf;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::fv::fullCorrectedSnGrad<Type>::correction
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
            fullCorrectedSnGrad<typename pTraits<Type>::cmptType>(mesh)
           .fullGradCorrection(vf.component(cmpt))
        );
    }

    return tssf;
}




// ************************************************************************* //
