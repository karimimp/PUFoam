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

\*---------------------------------------------------------------------------*/
#include "fvCFD.H"
#include "blowingReaction.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blowingReaction::blowingReaction
(
    fvMesh &mesh_, volScalarField &XW_, volScalarField &TS_
)
:
    mesh(mesh_),
    XW(XW_),
    TS(TS_)
{};

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::blowingReaction::~blowingReaction()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
bool Foam::blowingReaction::isBounded(scalar &X)
{
    bool answer;
    if (X < scalar(1.0) && X > scalar(0.0)){
        answer = true;
    }
    else
    {
        answer = false;
    }
    return (answer);
}
scalar Foam::blowingReaction::QKinW
(
    scalar& AW, scalar& EW, scalar& temp, scalar& Lliq, scalar& rhoPoly,
    scalar& rhoBL
)
{
    // QKinW - arrhenius term times by dilution term
    // @param - similar to arrehenius and BAdilution functions
    scalar R = 8.3145;
    return (
                AW*(Foam::exp(-EW/(R*temp)))*
                (1/(1+Lliq*(rhoPoly/(Foam::max(rhoBL,ROOTVSMALL)))))
           );
}

// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

void Foam::blowingReaction::operator=(const blowingReaction& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("Foam::blowingReaction::operator=(const Foam::blowingReaction&)")
            << "Attempted assignment to self"
            << abort(FatalError);
    }
}

// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


// ************************************************************************* //
