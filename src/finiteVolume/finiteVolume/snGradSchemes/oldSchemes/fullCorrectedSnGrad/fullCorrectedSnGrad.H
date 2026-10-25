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

Class
    Foam::fv::correctedSnGrad

Group
    grpFvSnGradSchemes

Description
    Surface gradient scheme with full explicit non-orthogonal and skew correction.

Usage
    Minimal example by using \c system/fvSchemes:
    \verbatim
    snGradSchemes
    {
        snGrad(<term>)       fullCorrected;
    }
    \endverbatim

SourceFiles
    fullCorrectedSnGrad.C

\*---------------------------------------------------------------------------*/

#ifndef fullCorrectedSnGrad_H
#define fullCorrectedSnGrad_H

#include "snGradScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                        Class correctedSnGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class fullCorrectedSnGrad
:
    public snGradScheme<Type>
{
    // Private Member Functions

        //- No copy assignment
        void operator=(const fullCorrectedSnGrad&) = delete;
       
	
	//tmp<GeometricField<typename outerProduct<vector, Type>::type, fvPatchField, volMesh>>
       // computeVolumeGradient
       // (
       //     const GeometricField<Type, fvPatchField, volMesh>& vf
       // ) const;

public:

    //- Runtime type information
    TypeName("fullCorrected");


    // Constructors

        //- Construct from mesh
        fullCorrectedSnGrad(const fvMesh& mesh)
        :
            snGradScheme<Type>(mesh)
        {}

        //- Construct from mesh and data stream
        fullCorrectedSnGrad(const fvMesh& mesh, Istream&)
        :
            snGradScheme<Type>(mesh)
        {}


    //- Destructor
    virtual ~fullCorrectedSnGrad() = default;


    // Member Functions

        //- Return the interpolation weighting factors for the given field
        virtual tmp<surfaceScalarField> deltaCoeffs
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const
        {
            return this->mesh().nonOrthDeltaCoeffs();
        }

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const noexcept
        {
            return true;
        }

        //- Return the explicit correction to the correctedSnGrad
        //- for the given field using the gradient of the field
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        fullGradCorrection
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        //- Return the explicit correction to the correctedSnGrad
        //- for the given field using the gradients of the field components
        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        correction(const GeometricField<Type, fvPatchField, volMesh>&) const;
};


// * * * * * * * * Template Member Function Specialisations  * * * * * * * * //

template<>
tmp<surfaceScalarField> fullCorrectedSnGrad<scalar>::correction
(
    const volScalarField& vsf
) const;


template<>
tmp<surfaceVectorField> fullCorrectedSnGrad<vector>::correction
(
    const volVectorField& vvf
) const;


template<>
tmp<surfaceTensorField> fullCorrectedSnGrad<tensor>::correction
(
    const volTensorField& vtf
) const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fullCorrectedSnGrad.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
