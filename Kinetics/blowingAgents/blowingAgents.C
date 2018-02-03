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
#include "blowingAgents.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::blowingAgents::blowingAgents
(
    fvMesh &mesh_, Time &runTime_, volScalarField &alpha2_, volScalarField &TS_
)
:
    PUgeneric(mesh_, runTime_, alpha2_, TS_)
{};


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::blowingAgents::~blowingAgents()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
scalar Foam::blowingAgents::n_pentaneInGas
(
    blowingAgents &ba, scalar Tc, scalar To
)
{
    return (ba.ddT_LliqMax(Tc)*ba.dQdt(Tc, To));

}
scalar Foam::blowingAgents::n_pentaneInLiquid
(
    scalar molecualrMassBA, scalar liquidDensity, scalar p, scalar R,
    scalar TS
)
{
    return (
             (molecualrMassBA/1000.0)*(1.0/liquidDensity)*
             (p/(max((R*TS),ROOTVSMALL)))
           );
}

scalar Foam::blowingAgents::R11InGas
(
    blowingAgents &ba, scalar molecualrMassBA, scalar molecularMassNCO, 
    scalar modelConstant, scalar Tc, scalar To
)
{
    return (
               (-(molecualrMassBA/molecularMassNCO)*
               (1.0/(Foam::pow((1.0 - ba.xBL(Tc, modelConstant)),2)))*
               (modelConstant))*ba.dQdt(Tc, To)
            );
}

scalar Foam::blowingAgents::CO2InGas
(
    scalar molecualrMass, scalar liquidDensity, scalar p, scalar R, scalar TS
)
{
    return(
               (molecualrMass/1000.0)*(1.0/liquidDensity)*
               (p/(max((R*TS),ROOTVSMALL)))
          );
}

scalar Foam::blowingAgents::CO2InLiquid
(
    blowingAgents &ba, scalar initW, scalar molecualrMass, scalar liquidDensity,
    scalar XWc,
    scalar XWo
)
{
    return(
            ba.dQdt(XWc, XWo)*initW*(molecualrMass/scalar(1000.0))*
            (scalar(1.0)/liquidDensity)
          );
}

scalar Foam::blowingAgents::CO2InGasNoPBE
(
    scalar initW, scalar molecualrMass, scalar liquidDensity, scalar XW
)
{
    return(
            initW*XW*(molecualrMass/scalar(1000.0))*
            (scalar(1.0)/liquidDensity)
          );
}



// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

void Foam::blowingAgents::operator=(const blowingAgents& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("Foam::blowingAgents::operator=(const Foam::blowingAgents&)")
            << "Attempted assignment to self"
            << abort(FatalError);
    }
}

// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


// ************************************************************************* //
