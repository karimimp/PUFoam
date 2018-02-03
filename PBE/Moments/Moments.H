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
    Foam::Moments

Description

SourceFiles
    Moments.C

\*---------------------------------------------------------------------------*/

#ifndef Moments_H
#define Moments_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class Moments;
Istream& operator>>(Istream&, Moments&);
Ostream& operator<<(Ostream&, const Moments&);


/*---------------------------------------------------------------------------*\
                         Class Moments Declaration
\*---------------------------------------------------------------------------*/

class Moments
{

public:

    // Constructors

        //- Construct null
        Moments();

    //- Destructor
    ~Moments();


    // Member Functions
    void PDA(double *we, double *vi, double *mom, int &n);
    scalar cc1Value(scalar L_l_val, scalar LMAX);
    void growthSource
    (
        scalar *sgBA, scalar *sgCO2, scalar *we, scalar *vi, int &nNodes,
        int *mOrder, scalar &CO2_l_val, scalar &L_l_val,
        scalar &wCO2_Max, scalar &LMAX
    );
    scalar growthRateConst();
    scalar growthRateCO2(scalar &CO2_l, scalar &wCO2_Max);
    scalar growthRateBA(scalar &L_l, scalar &LMAX);
    void coalescenceSource
    (
        double *sc, double *we, double *vi, int &nNodes, int *mOrder
    );
    scalar coalescenceKernel();


    // Member Operators

        void operator=(const Moments&);


    // Friend Functions

    // Friend Operators

    // IOstream Operators

        friend Istream& operator>>(Istream&, Moments&);
        friend Ostream& operator<<(Ostream&, const Moments&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
