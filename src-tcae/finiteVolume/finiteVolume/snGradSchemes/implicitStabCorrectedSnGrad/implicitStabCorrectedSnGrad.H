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
    Foam::fv::implicitStabCorrectedSnGrad

Group
    grpFvSnGradSchemes

Description
    Stabilised, positivity-preserving variant of implicitCorrectedSnGrad.

    The correction magnitude of the non-orthogonal and non-conjunctional
    corrections is added to the standard nonOrthDeltaCoeffs.  Unlike
    \c implicitCorrectedSnGrad the correction is added as a magnitude
    (never divided by (phi_N - phi_P)), so the effective coefficient is
    always at least the orthogonal one and the resulting Laplacian
    matrix remains diagonally dominant.  This trades some accuracy for
    robustness on strongly distorted meshes where (phi_N - phi_P) can be
    very small.

    Reference:  Hill et al., "Boundedness-preserving implicit correction
    of mesh-induced errors for VoF based heat and mass transfer",
    arXiv:1705.11063.

Usage
    In \c system/fvSchemes:
    \verbatim
    snGradSchemes
    {
        snGrad(<term>)    implicitStabCorrected;
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
                Class implicitStabCorrectedSnGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class implicitStabCorrectedSnGrad
:
    public snGradScheme<Type>
{
    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const implicitStabCorrectedSnGrad&) = delete;


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
    virtual ~implicitStabCorrectedSnGrad()
    {}


    // Member Functions

        //- Return the (implicitly stabilised) delta coefficients.
        virtual tmp<surfaceScalarField> deltaCoeffs
        (
            const VolField<Type>& vf
        ) const;

        //- No explicit correction (correction is implicit via delta).
        virtual bool corrected() const
        {
            return false;
        }

        //- Return the (null) explicit correction to snGrad
        virtual tmp<SurfaceField<Type>>
        correction(const VolField<Type>&) const
        {
            return tmp<SurfaceField<Type>>(nullptr);
        }
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
