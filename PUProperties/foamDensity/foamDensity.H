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
    Foam::foamDensity

Description

SourceFiles
    foamDensity.C

\*---------------------------------------------------------------------------*/

#ifndef foamDensity_H
#define foamDensity_H

#include "Time.H"
#include "PUgeneric.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class foamDensity;
Istream& operator>>(Istream&, foamDensity&);
Ostream& operator<<(Ostream&, const foamDensity&);


/*---------------------------------------------------------------------------*\
                         Class foamDensity Declaration
\*---------------------------------------------------------------------------*/

class foamDensity
:
    public PUgeneric
{
    // Protected data
    volScalarField& mOne;
    volScalarField& wCO2_g;
    volScalarField& wBA_g;
    volScalarField& wBA_l;
    volScalarField& p;
    volScalarField& XW;
    scalar L0;
    scalar rhoPoly;

public:

    // Constructors

        //- Construct from components
        foamDensity
        (
            fvMesh &mesh_, Time &runTime_, volScalarField &alpha2_,
            volScalarField &mOne_, volScalarField &wCO2_g_,
            volScalarField &wBA_g_, volScalarField &wBA_l_,
            volScalarField &p_, volScalarField &TS_, volScalarField &XW_,
            scalar L0_, scalar rhoPoly_
        );

    //- Destructor
    ~foamDensity();


    // Member Functions
    scalar LiquidMixtureDensity(scalar rhoPU);
    scalar InitialLiquidBlowingAgent(scalar initL);
    volScalarField rhoFoam(scalar M_CO2, scalar M_B);
    volScalarField rhoFoamNoPBE
    (
        scalar M_CO2, scalar M_B, scalar surfaceTension,
        scalar M_liq, scalar CW_0, scalar rhoBL
    );


    // Member Operators

        void operator=(const foamDensity&);


    // Friend Functions

    // Friend Operators

    // IOstream Operators

        friend Istream& operator>>(Istream&, foamDensity&);
        friend Ostream& operator<<(Ostream&, const foamDensity&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
