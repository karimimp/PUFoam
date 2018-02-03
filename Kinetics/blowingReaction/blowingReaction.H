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
    Foam::blowingReaction

Description

SourceFiles
    blowingReaction.C

\*---------------------------------------------------------------------------*/

#ifndef blowingReaction_H
#define blowingReaction_H

// #include ".H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class blowingReaction;
Istream& operator>>(Istream&, blowingReaction&);
Ostream& operator<<(Ostream&, const blowingReaction&);


/*---------------------------------------------------------------------------*\
                         Class blowingReaction Declaration
\*---------------------------------------------------------------------------*/

class blowingReaction
{
    // protected data
protected:
    const fvMesh& mesh;
    volScalarField &XW;
    volScalarField &TS;

public:


    // Constructors

        //- Construct from components
        blowingReaction(fvMesh &mesh_, volScalarField &XW_, volScalarField &TS_);


    //- Destructor
    ~blowingReaction();


    // Member Functions
    bool isBounded(scalar &X);
    scalar QKinW (scalar& AW, scalar& EW, scalar& temp, scalar& Lliq,
                  scalar& rhoPoly, scalar& rhoBL);

    // Member Operators

        void operator=(const blowingReaction&);

    // IOstream Operators

        friend Istream& operator>>(Istream&, blowingReaction&);
        friend Ostream& operator<<(Ostream&, const blowingReaction&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
