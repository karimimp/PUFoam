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
    Foam::gellingReaction

Description

SourceFiles
    gellingReactionI.H
    gellingReaction.C
    gellingReactionIO.C

\*---------------------------------------------------------------------------*/

#ifndef gellingReaction_H
#define gellingReaction_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class gellingReaction;
Istream& operator>>(Istream&, gellingReaction&);
Ostream& operator<<(Ostream&, const gellingReaction&);


/*---------------------------------------------------------------------------*\
                         Class gellingReaction Declaration
\*---------------------------------------------------------------------------*/

class gellingReaction
{
    // protected data
protected:
    const fvMesh& mesh;
    volScalarField &XOH;
    volScalarField &TS;

public:

    // Constructors

        //- Construct from components
        gellingReaction(fvMesh &mesh_, volScalarField &XOH_, 
            volScalarField &TS_);

    //- Destructor
    ~gellingReaction();


    // Member Functions
    scalar gellingSourceOldTime
           (
              scalar initOH, scalar initNCO, scalar initW, scalar XW, scalar XOH
           );

    // Member Operators

        void operator=(const gellingReaction&);


    // Friend Functions

    // Friend Operators

    // IOstream Operators

        friend Istream& operator>>(Istream&, gellingReaction&);
        friend Ostream& operator<<(Ostream&, const gellingReaction&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
