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
    Foam::PUgeneric

Description

SourceFiles
    PUgeneric.C

\*---------------------------------------------------------------------------*/

#ifndef PUgeneric_H
#define PUgeneric_H

#include "Time.H"
#include "PUgeneric.H"
    
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class PUgeneric;
Istream& operator>>(Istream&, PUgeneric&);
Ostream& operator<<(Ostream&, const PUgeneric&);


/*---------------------------------------------------------------------------*\
                         Class PUgeneric Declaration
\*---------------------------------------------------------------------------*/

class PUgeneric
{
    // Protected data

protected:
    const fvMesh& mesh;
    Time &runTime;
    volScalarField &alpha2;
    volScalarField &TS;

public:


    // Constructors

        //- Construct from components
        PUgeneric (const fvMesh &mesh_, Time &runTime_, volScalarField &alpha2_,
        volScalarField &TS_);


    //- Destructor
    ~PUgeneric();


    // Member Functions

        scalar LliqMax (scalar& tmptr);
        scalar dQdt(scalar Qc, scalar Qo);
        scalar ddT_LliqMax (double& tmptr);
        scalar xBL (scalar &T, scalar &dxdT);
        scalar wBL_D(scalar &xBL, scalar &M_B, scalar &M_NCO);
        scalar henryCoefficient(scalar &T);
        scalar bubbleRadius (const scalar m0, const scalar m1);
        scalar partialPressureCO2
        (
            scalar &M_CO2, scalar &M_B, scalar &surfaceTension,
            scalar &wCO2_g, double &wBA_g, double &p, double &R
        );
        scalar wCO2Max
        (
            scalar &M_CO2, scalar &M_liq,
            scalar &pCO2, scalar &henryCoeff
        );
        scalar creamTemperature(scalar &xBL0, scalar &dxdT);
        scalar arrhenius (scalar& AOH, scalar& EOH, scalar& tempt);
        scalar BAdilution(scalar& L_l, scalar& rhoPoly, scalar& rhoBL);
        scalar exothermicGelling
        (
            scalar& DH_OH, scalar& COH_0, scalar& dXOHdt, scalar& rhoPoly,
            scalar& C_TOT
        );
        scalar exothermicBlowing
        (
            scalar& DH_W, scalar& CW_0, scalar& dWdt, scalar& rhoPoly,
            scalar& C_TOT
        );
        scalar endothermicEvaporation
        (
            scalar& latenth, scalar& dLdt, scalar& rhoPoly, scalar& C_TOT
        );
        scalar R11EvaporationRate
        (
            scalar M_B, scalar M_NCO, scalar xBL_val, scalar dxdT, 
            scalar latenth
        );
        scalar thermalDiffusivityGas(scalar &T);
        scalar npentaneThermalconductivityHighrho(scalar &rho);
        scalar npentaneThermalconductivityLowrho(scalar &rho);

    // Member Operators

        void operator=(const PUgeneric& );


    // Friend Functions

    // Friend Operators

    // IOstream Operators

        friend Istream& operator>>(Istream&, PUgeneric&);
        friend Ostream& operator<<(Ostream&, const PUgeneric&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
