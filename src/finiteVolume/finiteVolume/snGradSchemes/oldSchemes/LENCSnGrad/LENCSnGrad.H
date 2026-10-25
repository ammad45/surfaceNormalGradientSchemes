/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
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
    Foam::fv::LENCSnGrad

SourceFiles
    LENCSnGrad.C

Authors
    Simon Hill         <simon.hill@linde.com>
    Daniel Deising     <deising@mma.tu-darmstadt.de>
    All rights reserved.

Description
    Simple central-difference scalar-type snGrad scheme with implicit 
    limited non-orthogonal and non-conjunctional correction.

    You may refer to this software as :
        S.Hill, D.Deising, T.Acher, H.Klein, D.Bothe, H.Marschall, 2018.
        "Boundedness-preserving implicit correction of mesh-induced 
         errors for VOF based heat and mass transfer", 
        J. Comput. Phys. vol. 352, pp. 285-300.

    This code has been developed by :
        Simon Hill (main developer).
    
    Method Development and Intellectual Property :
        Simon Hill       <simon.hill@linde.com>
        Daniel Deising   <deising@mma.tu-darmstadt.de>
        Holger Marschall <marschall@mma.tu-darmstadt.de>
        Dieter Bothe     <bothe@mma.tu-darmstadt.de>

        Linde Engineering AG, Pullach/Germany
        Mathematical Modeling and Analysis
        Center of Smart Interfaces
        Technische Universitaet Darmstadt
       
    If you use this software for your scientific work or your publications,
    please don't forget to acknowledge explicitly the use of it.

\*---------------------------------------------------------------------------*/

#ifndef LENCSnGrad_H
#define LENCSnGrad_H

#include "snGradScheme.H"
// new
#include "linear.H"
#include "fvcGrad.H"
#include "skewCorrectionVectors.H"
#include "dimensionedType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                 Class LENCSnGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class LENCSnGrad
:
    public snGradScheme<Type>
{
    // Private Member Data
        
        //- Explicit Non-Orthogonal correction field
        // Has to be mutable that the const function deltaCoeffs
        // is allowed to change this variable
        mutable tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tnonOrthCorr_;
        
    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const LENCSnGrad&);


public:

    //- Runtime type information
    TypeName("LENC");


    // Constructors

        //- Construct from mesh
        LENCSnGrad(const fvMesh& mesh)
        :
            snGradScheme<Type>(mesh),
            tnonOrthCorr_
            (
                new surfaceScalarField
                (
                    IOobject
                    (
                        "nonOrthCorr",
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    mesh,
                    dimensionSet(0, -1, 0, 0, 0, 0, 0)
                )
            )
        {}


        //- Construct from mesh and data stream
        LENCSnGrad(const fvMesh& mesh, Istream&)
        :
            snGradScheme<Type>(mesh),
            tnonOrthCorr_
            (
                new surfaceScalarField
                (
                    IOobject
                    (
                        "nonOrthCorr",
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    mesh,
                    dimensionSet(0, -1, 0, 0, 0, 0, 0)
                )
            )
        {}


    //- Destructor
    virtual ~LENCSnGrad();

    // Member Functions

        //- Return the interpolation weighting factors for the given field
        virtual tmp<surfaceScalarField> deltaCoeffs
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const;

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const noexcept 
        {
            return true;
        }

        //- Return the explicit correction to the LENCSnGrad
        //  for the given field using the gradients of the field components
        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
        correction(const GeometricField<Type, fvPatchField, volMesh>&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "LENCSnGrad.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
