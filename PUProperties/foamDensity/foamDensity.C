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
#include "foamDensity.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::foamDensity::foamDensity
(
    fvMesh &mesh_, Time &runTime_, volScalarField &alpha2_,
    volScalarField &mOne_, volScalarField &wCO2_g_,
    volScalarField &wBA_g_, volScalarField &wBA_l_,
    volScalarField &p_, volScalarField &TS_, volScalarField &XW_,
    scalar L0_, scalar rhoPoly_
)
:
    PUgeneric(mesh_, runTime_, alpha2_, TS_),
    mOne(mOne_),
    wCO2_g(wCO2_g_),
    wBA_g(wBA_g_),
    wBA_l(wBA_l_),
    p(p_),
    XW(XW_)
{
    L0 = InitialLiquidBlowingAgent(L0_);
    rhoPoly = LiquidMixtureDensity(rhoPoly_);
};



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::foamDensity::~foamDensity()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
scalar Foam::foamDensity::LiquidMixtureDensity(scalar rhoPU)
{
    rhoPoly = rhoPU;
    return rhoPoly;
}

scalar Foam::foamDensity::InitialLiquidBlowingAgent(scalar initL)
{
    L0 = initL;
    return L0;
}

volScalarField Foam::foamDensity::rhoFoam(scalar M_CO2, scalar M_B)
{
    scalar R = 8.3145;

    volScalarField rho_bubble_
    (
        IOobject
        (
            "rho_bubble_",
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("rho_bubble_", dimensionSet(1,-1,-2,1,0,0,0), 0.0)
    );
    volScalarField rho_foam_
    (
        IOobject
        (
            "rho_foam_",
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("rho_foam_", dimensionSet(1,-3,0,0,0,0,0), SMALL)
    );

    forAll(mesh.C(), celli)
    {
        rho_bubble_[celli] =
            (p[celli]/(R*TS[celli]))*
            ((wCO2_g[celli]*M_CO2+wBA_g[celli]*M_B)/
            Foam::max((scalar(1000.0)*(wCO2_g[celli]+wBA_g[celli])),ROOTVSMALL));

        rho_foam_[celli] =
           (
            (rho_bubble_[celli]*(mOne[celli]/(scalar(1.0) + mOne[celli]))
            + ((scalar(1.0) + L0)*(rhoPoly - L0*rhoPoly)*(scalar(1.0) -
            (mOne[celli]/(scalar(1.0) + mOne[celli])))))
           );
    }

    return (rho_foam_);
}

volScalarField Foam::foamDensity::rhoFoamNoPBE
(
    scalar M_CO2, scalar M_B, scalar surfaceTension,
    scalar M_liq, scalar CW_0, scalar rhoBL
)
{
    scalar R = 8.3145;
    volScalarField rho_foam_
    (
        IOobject
        (
            "rho_foam_",
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("rho_foam_", dimensionSet(1,-3,0,0,0,0,0), SMALL)
    );
    volScalarField CO_2
    (
        IOobject
        (
            "CO_2",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("CO_2", dimless, 0.0)
    );
    forAll(mesh.C(), celli)
    {
        scalar henry_coeff = henryCoefficient(TS[celli]);
        scalar bubble_radius = bubbleRadius(0.0, 0.0);
        scalar partialPressure_CO2 =
            partialPressureCO2
            (
                M_CO2, M_B, surfaceTension, wCO2_g[celli],
                wBA_g[celli], p[celli], bubble_radius
            );
        scalar wCO2_Max =
            wCO2Max(M_CO2, M_liq, partialPressure_CO2, henry_coeff);

        CO_2[celli] =
            (
                ((CW_0*XW[celli]*M_CO2)/(scalar(1000.0)*(rhoPoly))) - wCO2_Max
            );

        if (CO_2[celli] < 0.0)
        {
            CO_2[celli] = 0.0;
        }
        rho_foam_[celli] =
            (
                scalar(1.0)+L0)/
                (((CO_2[celli]*scalar(1000.0)*R*TS[celli])/(p[celli]*M_CO2))
             +  ((wBA_g[celli]*scalar(1000.0)*R*TS[celli])/(p[celli]*M_B))
             +  ((wBA_l[celli])/(rhoBL))
             +  (scalar(1.0)/(rhoPoly))
             );
    }
    return (rho_foam_);
}
// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

void Foam::foamDensity::operator=(const foamDensity& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("Foam::foamDensity::operator=(const Foam::foamDensity&)")
            << "Attempted assignment to self"
            << abort(FatalError);
    }
}

// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


// ************************************************************************* //
