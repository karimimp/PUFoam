/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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
    Foam::rheologyPU

Description

SourceFiles
    rheologyPU.C

\*---------------------------------------------------------------------------*/

#ifndef rheologyPU_H
#define rheologyPU_H

#include "Time.H"
#include "PUgeneric.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class rheologyPU;
Istream& operator>>(Istream&, rheologyPU&);
Ostream& operator<<(Ostream&, const rheologyPU&);


/*---------------------------------------------------------------------------*\
                         Class rheologyPU Declaration
\*---------------------------------------------------------------------------*/

class rheologyPU
:
    public PUgeneric
{
    // Protected data
protected:
    volScalarField &XW;
    volScalarField &XOH;
    volVectorField &U;
    scalar XOH_Gel;

public:

    // Constructors

        //- Construct from components
        rheologyPU
        (
            fvMesh &mesh_, Time &runTime_, volScalarField &alpha2_,
            volScalarField &TS_, volScalarField &XW_, volScalarField &XOH_,
            volVectorField &U_, scalar XOH_Gel_
        );

    //- Destructor
    ~rheologyPU();

    // Member Functions
    scalar setGelPoint(scalar gelPoint);
    volScalarField waterLike();
    volScalarField CastroMacosko
                   (
                     scalar initOH, scalar initNCO,scalar initW
                   );
    volScalarField BirdCarreau
                   (
                        volScalarField &mu0, volScalarField &muinf,
                        volScalarField &muFoam
                   );

    // Member Operators

        void operator=(const rheologyPU&);


    // Friend Functions

    // Friend Operators

    // IOstream Operators

        friend Istream& operator>>(Istream&, rheologyPU&);
        friend Ostream& operator<<(Ostream&, const rheologyPU&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
