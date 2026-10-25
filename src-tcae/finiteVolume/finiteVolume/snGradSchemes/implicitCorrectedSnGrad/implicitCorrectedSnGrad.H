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
    Foam::fv::implicitCorrectedSnGrad

Group
    grpFvSnGradSchemes

Description
    Surface-normal gradient scheme with a fully implicit correction for
    both mesh non-orthogonality and non-conjunctionality (skewness). The
    correction is folded into the face delta coefficient so that the
    Laplacian operator becomes bounded and diagonally dominant without
    an explicit deferred-correction source term.

    The construction follows the semi-implicit / NO/NC scheme described in
    Hill et al., "Boundedness-preserving implicit correction of
    mesh-induced errors for VoF based heat and mass transfer",
    arXiv:1705.11063, equations (17)-(19).

    Given the standard corrected surface-normal gradient decomposition

        snGrad = deltaOrth * (phi_N - phi_P)
               + kHat . (grad phi)_f'                     (non-orthogonal)
               + nHat . ( Hessian . m )                   (non-conjunctional)

    with kHat = nonOrthCorrectionVectors, nHat = Sf/|Sf| and m the skew
    vector Cf - f', the effective delta coefficient is

        deltaEff = deltaOrth
                 + [ kHat.(grad phi)_f' + nHat.(Hessian.m) ] / (phi_N - phi_P).

    Corrections are skipped when |phi_N - phi_P| is negligible.  Each
    component increment is limited so deltaEff never turns negative;
    faces where the correction would (nearly) cancel the orthogonal
    coefficient fall back to deltaOrth so the laplacian matrix stays
    diagonally valid for GAMG.

Usage
    Minimal example, in \c system/fvSchemes:
    \verbatim
    snGradSchemes
    {
        snGrad(<term>)   implicitCorrected;
    }
    \endverbatim

SourceFiles
    implicitCorrectedSnGrad.C

\*---------------------------------------------------------------------------*/

#ifndef implicitCorrectedSnGrad_H
#define implicitCorrectedSnGrad_H

#include "snGradScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                   Class implicitCorrectedSnGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class implicitCorrectedSnGrad
:
    public snGradScheme<Type>
{
    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const implicitCorrectedSnGrad&) = delete;


public:

    //- Runtime type information
    TypeName("implicitCorrected");


    // Constructors

        //- Construct from mesh
        implicitCorrectedSnGrad(const fvMesh& mesh)
        :
            snGradScheme<Type>(mesh)
        {}

        //- Construct from mesh and data stream
        implicitCorrectedSnGrad(const fvMesh& mesh, Istream&)
        :
            snGradScheme<Type>(mesh)
        {}


    //- Destructor
    virtual ~implicitCorrectedSnGrad()
    {}


    // Member Functions

        //- Return the (implicitly corrected) delta coefficients for the
        //  given field. The non-orthogonal and non-conjunctional
        //  corrections are folded into the returned coefficient.
        virtual tmp<surfaceScalarField> deltaCoeffs
        (
            const VolField<Type>& vf
        ) const;

        //- No explicit correction (the correction is implicit via
        //  deltaCoeffs above).
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
    #include "implicitCorrectedSnGrad.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
