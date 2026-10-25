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
    Foam::fv::implicitStabCorrectedSnGrad

Group
    grpFvSnGradSchemes

Description
    Stabilised, positivity-preserving variant of implicitCorrectedSnGrad.

    The magnitude of the non-orthogonal and non-conjunctional corrections
    is added to the standard nonOrthDeltaCoeffs.  Unlike
    \c implicitCorrectedSnGrad the correction is added as a magnitude
    (never divided by (phi_N - phi_P)), so the effective coefficient is
    always at least the orthogonal one and the resulting Laplacian
    matrix remains diagonally dominant.  Coupled patches are also
    corrected.

Usage
    Minimal example by using \c system/fvSchemes:
    \verbatim
    snGradSchemes
    {
        snGrad(<term>)       implicitStabCorrected;
    }
    \endverbatim

SourceFiles
    implicitStabCorrectedSnGrad.C

\*---------------------------------------------------------------------------*/

#ifndef implicitStabCorrectedSnGrad_H
#define implicitStabCorrectedSnGrad_H

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
class implicitStabCorrectedSnGrad
:
    public snGradScheme<Type>
{

    //tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tnonOrthCorr_;
    // Private Member Functions

        //- No copy assignment
        void operator=(const implicitStabCorrectedSnGrad&) = delete;
       
	
	//tmp<GeometricField<typename outerProduct<vector, Type>::type, fvPatchField, volMesh>>
       // computeVolumeGradient
       // (
       //     const GeometricField<Type, fvPatchField, volMesh>& vf
       // ) const;

    public:

    //- Runtime type information
    TypeName("implicitStabCorrected");


    // Constructors

        //- Construct from mesh
        implicitStabCorrectedSnGrad(const fvMesh& mesh)
        :
            snGradScheme<Type>(mesh)
        {}

        //- Construct from mesh and data stream
        implicitStabCorrectedSnGrad(const fvMesh& mesh, Istream&)
        :
            snGradScheme<Type>(mesh)
        {}


    //- Destructor
    virtual ~implicitStabCorrectedSnGrad() = default;


    // Member Functions

        //- Return the interpolation weighting factors for the given field
        virtual tmp<surfaceScalarField> deltaCoeffs
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const;
      // {
      //     return this->mesh().nonOrthDeltaCoeffs();
      // }

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const noexcept
        {
            return false;
        }

        //- Return the explicit correction to the correctedSnGrad
        //- for the given field using the gradient of the field
       // tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
       // fullGradCorrection
       // (
       //     const GeometricField<Type, fvPatchField, volMesh>&
       // ) const;

        //- Return the explicit correction to the correctedSnGrad
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
    #include "implicitStabCorrectedSnGrad.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
