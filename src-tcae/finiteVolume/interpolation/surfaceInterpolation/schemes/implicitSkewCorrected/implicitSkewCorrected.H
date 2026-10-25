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

Class
    Foam::implicitSkewCorrected

Description
    Face-interpolation scheme with an *implicit* non-conjunctionality
    (skew) correction, in the sense of Hill et al.,
    arXiv:1705.11063, Eq. 14.  The correction

        (grad phi)_{f'} . m           (paper Eq. 3)

    with  m = C_f − f'  (the skew vector, i.e. the offset of the
    face centre from the intersection of the P-N line with the face) is
    folded into the interpolation weight of a wrapped base scheme
    rather than being added as an explicit deferred-correction source
    term.

    For a **scalar** field phi with a base weight w_base the effective
    (implicit) interpolation weight is

        w_eff = w_base + ( m . (grad phi)_{f'} ) / (phi_N - phi_P)

    which is exactly the paper's Eq. 14 (without the TVD limiter, which
    is not meaningful for velocity / turbulence fields).

    For a **vector / tensor** field the correction is a Type-valued
    quantity while OpenFOAM's implicit convection machinery
    (gaussConvectionScheme) only consumes a **single scalar weight per
    face**.  The natural scalar reduction is the projection of the
    correction onto the direction of the neighbour difference:

        w_eff = w_base
              + ( (m . grad U)_{f'} . (U_N - U_P) ) / |U_N - U_P|^2

    which is exact when U_N - U_P is aligned with the correction and
    reduces to the scalar formula above component-wise.  The tangential
    component of the correction (orthogonal to U_N - U_P) is dropped;
    this is a genuine approximation for vectors, but it keeps the
    scheme *fully implicit* (no explicit source term is introduced).

    Boundedness / stability treatment (identical to the diffusion
    scheme implicitCorrectedSnGrad):

      * If |vf_N - vf_P| < SMALL for a face, the correction is
        skipped for that face and the base weight is used unchanged –
        prevents division blow-up in near-uniform regions and inside
        stagnation zones.
      * The base scheme's own explicit correction()/corrected() is
        passed through, so this scheme can be layered on top of e.g.
        linearUpwind and gain second-order upwind reconstruction
        without duplicating that work.

    Suggested combinations for turbomachinery:
    \verbatim
    divSchemes
    {
        default              none;
        div(phi,U)           Gauss implicitSkewCorrected linearUpwindV grad(U);
        div(phi,k)           Gauss implicitSkewCorrected linearUpwind  grad(k);
        div(phi,omega)       Gauss implicitSkewCorrected linearUpwind  grad(omega);
    }
    gradSchemes
    {
        default              leastSquares;
        grad(U)              leastSquares;
    }
    \endverbatim

SourceFiles
    implicitSkewCorrected.C

\*---------------------------------------------------------------------------*/

#ifndef implicitSkewCorrected_H
#define implicitSkewCorrected_H

#include "surfaceInterpolationScheme.H"
#include "skewCorrectionVectors.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class implicitSkewCorrected Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class implicitSkewCorrected
:
    public surfaceInterpolationScheme<Type>
{
    // Private Data

        //- Wrapped base interpolation scheme (linear, linearUpwind, ...)
        tmp<surfaceInterpolationScheme<Type>> tScheme_;


public:

    //- Runtime type information
    TypeName("implicitSkewCorrected");


    // Constructors

        //- Construct from mesh and Istream
        implicitSkewCorrected
        (
            const fvMesh& mesh,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            tScheme_
            (
                surfaceInterpolationScheme<Type>::New(mesh, is)
            )
        {}

        //- Construct from mesh, faceFlux and Istream
        implicitSkewCorrected
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            tScheme_
            (
                surfaceInterpolationScheme<Type>::New(mesh, faceFlux, is)
            )
        {}

        //- Disallow default bitwise copy construction
        implicitSkewCorrected(const implicitSkewCorrected&) = delete;


    // Member Functions

        //- Return the (implicit, skew-corrected) interpolation weights
        //  for the given field.
        virtual tmp<surfaceScalarField> weights
        (
            const VolField<Type>& vf
        ) const;

        //- Passes-through the base scheme's `corrected()`.
        //  If the wrapped base scheme requires an explicit correction
        //  (e.g. linearUpwind's gradient reconstruction), we keep it.
        virtual bool corrected() const
        {
            return tScheme_().corrected();
        }

        //- Passes-through the base scheme's explicit correction.
        virtual tmp<SurfaceField<Type>>
        correction(const VolField<Type>& vf) const
        {
            return tScheme_().correction(vf);
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const implicitSkewCorrected&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "implicitSkewCorrected.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
