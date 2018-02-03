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
    Foam::blowingAgents

Description

SourceFiles
    blowingAgents.C

\*---------------------------------------------------------------------------*/

#ifndef blowingAgents_H
#define blowingAgents_H

#include "Time.H"
#include "PUgeneric.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class blowingAgents;
Istream& operator>>(Istream&, blowingAgents&);
Ostream& operator<<(Ostream&, const blowingAgents&);


/*---------------------------------------------------------------------------*\
                         Class blowingAgents Declaration
\*---------------------------------------------------------------------------*/

class blowingAgents
:
    public PUgeneric
{

public:

    // Constructors

    //- Construct from components
    blowingAgents(fvMesh &mesh_, Time &runTime_, volScalarField &alpha2_,
    volScalarField &TS_);

    //- Destructor
    ~blowingAgents();


    // Member Functions
    scalar n_pentaneInGas(blowingAgents &ba, scalar Tc, scalar To);
    scalar n_pentaneInLiquid
    (
        scalar molecualrMassBA, scalar liquidDensity, scalar p, scalar R,
        scalar TS
    );
    scalar R11InGas
    (
        blowingAgents &ba, scalar molecualrMassBA, scalar molecularMassNCO,
        scalar modelConstant, scalar Tc, scalar To
    );
    scalar CO2InGas
   (
        scalar molecualrMass, scalar liquidDensity, scalar p, scalar R,
        scalar TS
   );
   scalar CO2InLiquid
   (
        blowingAgents &ba, scalar initW, scalar molecualrMass,
        scalar liquidDensity, scalar XWc, scalar XWo
   );
   scalar CO2InGasNoPBE
   (
        scalar initW, scalar molecualrMass, scalar liquidDensity, scalar XW
   );

    // Member Operators

        void operator=(const blowingAgents&);


    // Friend Functions

    // Friend Operators

    // IOstream Operators

        friend Istream& operator>>(Istream&, blowingAgents&);
        friend Ostream& operator<<(Ostream&, const blowingAgents&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
