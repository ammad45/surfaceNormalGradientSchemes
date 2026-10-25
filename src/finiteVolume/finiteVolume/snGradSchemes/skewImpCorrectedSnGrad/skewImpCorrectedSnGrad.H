/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 Zeljko Tukovic, FSB Zagreb.
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
    Foam::fv::skewImpCorrectedSnGrad

Group
    grpFvSnGradSchemes

Description
    Surface gradient scheme with skewness and full
    explicit non-orthogonal corrections.

Usage
    Minimal example by using \c system/fvSchemes:
    \verbatim
    snGradSchemes
    {
        snGrad(<term>)       skewCorrected;
    }
    \endverbatim

SourceFiles
    skewImpCorrectedSnGrad.C

\*---------------------------------------------------------------------------*/

#ifndef skewImpCorrectedSnGrad_H
#define skewImpCorrectedSnGrad_H

#include "snGradScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                     Class skewImpCorrectedSnGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class skewImpCorrectedSnGrad
:
    public snGradScheme<Type>
{
    // Private Member Functions
        //mutable tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tnonOrthCorr_;
        //- No copy assignment
        void operator=(const skewImpCorrectedSnGrad&) = delete;


public:

    //- Runtime type information
    TypeName("skewImpCorrected");


    // Constructors

        //- Construct from mesh
          //- Construct from mesh
        skewImpCorrectedSnGrad(const fvMesh& mesh)
        :
            snGradScheme<Type>(mesh)   
        {}


        //- Construct from mesh and data stream
        skewImpCorrectedSnGrad(const fvMesh& mesh, Istream&)
        :
            snGradScheme<Type>(mesh)
        {}



    //- Destructor
    virtual ~skewImpCorrectedSnGrad() = default;


    // Member Functions

        //- Return the interpolation weighting factors for the given field
        virtual tmp<surfaceScalarField> deltaCoeffs
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf 
        ) const;
        //{
        //    return this->mesh().nonOrthDeltaCoeffs();
        //}

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const noexcept
        {
            return false;
        }

        //- Return the explicit correction to the skewImpCorrectedSnGrad
        //- for the given field using the gradient of the field
        //tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        //fullGradCorrection
        //(
        //    const GeometricField<Type, fvPatchField, volMesh>&
        //) const;

        //- Return the explicit correction to the skewImpCorrectedSnGrad
        //- for the given field using the gradients of the field components
        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        correction(const GeometricField<Type, fvPatchField, volMesh>&) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "skewImpCorrectedSnGrad.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
