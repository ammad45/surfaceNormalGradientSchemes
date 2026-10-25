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

Class
    Foam::fv::newCorrectedSnGrad

Group
    grpFvSnGradSchemes

Description
    Explicit surface-normal gradient scheme with both non-orthogonal
    ("k.grad phi") and non-conjunctional (skew, "nHat.Hessian.m")
    corrections applied as a deferred-correction source term.

    Equivalent to the built-in \c corrected scheme extended with the
    additional skewness correction from arXiv:1705.11063.

Usage
    In \c system/fvSchemes:
    \verbatim
    snGradSchemes
    {
        snGrad(<term>)   newCorrected;
    }
    \endverbatim

SourceFiles
    newCorrectedSnGrad.C

\*---------------------------------------------------------------------------*/

#ifndef newCorrectedSnGrad_H
#define newCorrectedSnGrad_H

#include "snGradScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                    Class newCorrectedSnGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class newCorrectedSnGrad
:
    public snGradScheme<Type>
{
    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const newCorrectedSnGrad&) = delete;


public:

    //- Runtime type information
    TypeName("newCorrected");


    // Constructors

        //- Construct from mesh
        newCorrectedSnGrad(const fvMesh& mesh)
        :
            snGradScheme<Type>(mesh)
        {}

        //- Construct from mesh and data stream
        newCorrectedSnGrad(const fvMesh& mesh, Istream&)
        :
            snGradScheme<Type>(mesh)
        {}


    //- Destructor
    virtual ~newCorrectedSnGrad()
    {}


    // Member Functions

        //- Return the interpolation weighting factors for the given
        //  field: use the standard over-relaxed non-orthogonal delta.
        virtual tmp<surfaceScalarField> deltaCoeffs
        (
            const VolField<Type>&
        ) const
        {
            return this->mesh().nonOrthDeltaCoeffs();
        }

        //- This scheme uses an explicit correction
        virtual bool corrected() const
        {
            return true;
        }

        //- Return the explicit correction (non-orthogonal + skew) for
        //  the whole VolField (built up component-wise).
        virtual tmp<SurfaceField<Type>>
        correction(const VolField<Type>&) const;

        //- Return the correction for a single component. Called via
        //  the recursion below.
        tmp<SurfaceField<Type>>
        fullGradCorrection(const VolField<Type>&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "newCorrectedSnGrad.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
