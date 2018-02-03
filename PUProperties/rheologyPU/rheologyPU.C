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
#include "rheologyPU.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::rheologyPU::rheologyPU
(
    fvMesh &mesh_, Time &runTime_, volScalarField &alpha2_,
    volScalarField &TS_, volScalarField &XW_, volScalarField &XOH_,
    volVectorField &U_, scalar XOH_Gel_
)
:
    PUgeneric(mesh_, runTime_, alpha2_, TS_),
    XW(XW_),
    XOH(XOH_),
    U(U_)
{
    XOH_Gel = setGelPoint(XOH_Gel_);
};



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::rheologyPU::~rheologyPU()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
scalar Foam::rheologyPU::setGelPoint(scalar gelPoint)
{
    XOH_Gel = gelPoint;
    return XOH_Gel;
}

volScalarField Foam::rheologyPU::waterLike()
{
    volScalarField muFoamCorr_
    (
        IOobject
        (
        "muFoamCorr_",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("muFoamCorr_", dimensionSet(1,-1,-1,0,0,0,0), 1e-5)
   );
   forAll(mesh.C(), celli)
   {
        if (alpha2[celli] > 0.5)
        {
            muFoamCorr_[celli] = 1.0e-3;
        }
   }
    return (muFoamCorr_);
}

volScalarField Foam::rheologyPU::CastroMacosko
(
    scalar initOH, scalar initNCO,scalar initW
)
{
    // model constants
    scalar muinf = 10.3*1e-8;
    scalar XNCOGel = 1.0;
    scalar EnuR = 4970.0;
    scalar nua = 1.5;
    scalar nub = 1.0;
    scalar nuc = 0.0;
    volScalarField muFoamCorr_
    (
        IOobject
        (
        "muFoamCorr_",
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("muFoamCorr_",
            dimensionSet(1,-1,-1,0,0,0,0), 1e-5)
   );
    volScalarField XNCO_
    (
        IOobject
        (
        "XNCO_",
         mesh,
         IOobject::READ_IF_PRESENT,
         IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("XNCO_", dimless, SMALL)
    );
    forAll(mesh.C(), celli)
    {
        if (alpha2[celli] > 0.5)
        {
            XNCO_[celli] =
                (initNCO - (initNCO - XOH[celli]*initOH -
                2.0*XW[celli]*initW))/initNCO;
            muFoamCorr_[celli] =
                (
                    muinf*Foam::exp(EnuR/(TS[celli]))*
                    Foam::pow((XNCOGel/max((XNCOGel -
                    XNCO_[celli]),ROOTVSMALL)),(nua +
                    nub*XNCO_[celli] + nuc*XNCO_[celli]*
                    XNCO_[celli]))
                );
        }
    }
    return (muFoamCorr_);
}

volScalarField Foam::rheologyPU::BirdCarreau
(
    volScalarField &mu0, volScalarField &muinf,
    volScalarField &muFoam
)
{
    // model constants
    scalar Amu = 0.0387;
    scalar EmuR = 10000.0;
    scalar a = 1.5;
    scalar b = 1.0;
    scalar c = 0.0;
    scalar d = 0.001;
    scalar mu0l = 0.195;
    scalar muinfl = 0.266;
    scalar coeffalpha = 2.0;
    scalar lambdaBC = 11.35;
    scalar coeffn = 0.2;
    const scalar R = 8.3145;

    volScalarField shearRate = Foam::mag(fvc::grad(U));

    volScalarField muFoamCorr_
    (
        IOobject
        (
        "muFoamCorr_",
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("muFoamCorr_",
            dimensionSet(1,-1,-1,0,0,0,0), 1e-5)
   );

    forAll(mesh.C(), celli)
    {
        if(alpha2[celli] > 0.5)
        {
            mu0[celli] =
                (Foam::log(XOH[celli] + d) - Foam::log(d) +
                Foam::pow((XOH_Gel/max((XOH_Gel - XOH[celli])
                ,0.00001)),((a + b*XOH[celli] + c*XOH[celli]
                *XOH[celli]))))*mu0l;
            muinf[celli] =
                (Foam::log(XOH[celli] + d) - Foam::log(d) +
                Foam::pow((XOH_Gel/max((XOH_Gel-XOH[celli]),
                0.00001)),((a + b*XOH[celli] + c*XOH[celli]
                *XOH[celli]))))*muinfl;
            muFoam[celli] =
                ((muinf[celli] + (mag((mu0[celli] -
                muinf[celli]))*Foam::pow((1.0 +
                Foam::pow((shearRate[celli]*lambdaBC),
                coeffalpha)),((coeffn-1)/coeffalpha)))));
            muFoamCorr_[celli] =
                (muFoam[celli]*
                (Amu*Foam::exp(EmuR/R/TS[celli])));
        }
        else
        {
            muFoamCorr_[celli] = 1.0e-5;
        }
    }
    return (muFoamCorr_);
}
// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

void Foam::rheologyPU::operator=(const rheologyPU& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("Foam::rheologyPU::operator=(const Foam::rheologyPU&)")
            << "Attempted assignment to self"
            << abort(FatalError);
    }
}

// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


// ************************************************************************* //
