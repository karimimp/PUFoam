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
#include "Moments.H"

extern "C"{void dsteqr_(char &, int *, double *, double *, double *, int *, double *, int *); }

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Moments::Moments()
{};


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Moments::~Moments()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::Moments::PDA(double *we, double *vi, double *mom, int &n)
{
    // Construct P Matrix
    int i,j;

    double p[2*n+1][2*n+1];

    for(i = 0; i < 2*n+1; i++)
    {
        for(j = 0;j < 2*n+1; j++)
        {
            p[i][j] = 0.0;
        }
    }

    double norm_mom[2*n];

    // Normalize the moments
    for(i = 0; i < 2*n; i++)
    {
        norm_mom[i] = mom[i]/(max(mom[0],1.0e-10));
    }

    // First column of P matrix
    p[0][0] = 1.0;
    p[0][1] = 1.0;

    // Second column of P matrix
    for(i = 1; i < 2*n; i++)
    {
        p[i][1] = Foam::pow(-1,double (i))*norm_mom[i];
    }

    // Recursion method for calculation of P
    for(j = 2; j < 2*n+1; j++)
    {
        for(i = 0; i < 2*n+2-j; i++)
            {
                p[i][j] = p[0][j-1]*p[i+1][j-2] - p[0][j-2]*p[i+1][j-1];
            }
    }

    // Computing zeta
    double zeta[2*n];
    for( i = 0; i < 2*n; i++)
    {
        zeta[i] = 0.0;
    }

    for(i = 1; i < 2*n; i++)
    {
        if(p[0][i]*p[0][i-1]>0)
        {
            zeta[i] = p[0][i+1]/(p[0][i]*p[0][i-1]);
        }
        else
        {
            zeta[i] = 0.0;
        }
    }

    // Coefficients of Jacobi matrix
    double aa[n];
    double bb[n-1], cc[n-1];

    for(i = 0; i < n-1; i++)
    {
        bb[i] = cc[i] = 0.0;
    }

    for(i = 1;i <= n; i++)
    {
        aa[i-1]=zeta[2*i-1]+zeta[2*i-2];
    }

    for(i = 1; i <= n-1; i++)
    {
        bb[i-1]=zeta[2*i]*zeta[2*i-1];
    }

    for(i = 1;i <= n-1; i++)
    {
        cc[i-1]= Foam::sqrt(fabs(bb[i-1]));
    }


    // Eigenvalues and Eigenvectors of Jacobi matrix
    double evec[n][n];
    for(i = 0; i < n;i++)
    {
        for(j = 0;j < n; j++)
        {
            evec[i][j] = 0.0;
        }
    }


    double work[2*n-2];
    int info;
    char choice='I';

    dsteqr_(choice,&n,aa,cc,&evec[0][0],&n,work,&info);

    for(i = 0; i < n; i++)
    {
        vi[i] = aa[i];

        if(vi[i] < 0.0)
        {
            vi[i] = 0.0;
        }
    }

    for(j = 0; j < n; j++)
    {
        we[j] = mom[0]*Foam::pow(evec[j][0],2);

        if(we[j] < 0)
        {
            we[j] = 0.0;
        }
    }
} // End of PDA

scalar Foam::Moments::cc1Value(scalar L_l_val, scalar LMAX)
{
    scalar cc1;
    if (L_l_val > LMAX)
    {
        cc1 = 1.0;
    }
    else
    {
        cc1 = 0.0;
    }
    return cc1;
}

void Foam::Moments::growthSource
(
    scalar *sgBA, scalar *sgCO2, scalar *we, scalar *vi, int &nNodes,
    int *mOrder, scalar &CO2_l_val, scalar &L_l_val,
    scalar &wCO2_Max, scalar &LMAX
)
{
    // growthSource - source term for the bubble growth
    // @param - double *sgBA - source of growth due to blowing agent
    // @param - double *sgCO2 - source of growth due to CO2
    // @param - double *we - weights of quadrature
    // @param - duoble *vi - nodes of quadrature
    // @param - int &nNodes - number of nodes
    // @param - int *mOrder - order of moments
    // @param - double &CO2_l_val - weight fraction of CO2 in liquid
    // @param - double &L_l_val - weight fraction of liquid blowing agent in liquid
    // @param - double &tmp_val - temperature
    // @param - double &wCO2_Max - weight fraction of maximum allowable CO2 in liquid
    // @param - double &cc1_val - constant
    // @param - double &LMAX - weight fraction of maximum allowable blowing agent in liquid
    int i;
    int counter = 0;
    scalar k;
    scalar c_L, c_CO2;

    c_CO2 = growthRateCO2(CO2_l_val, wCO2_Max);     // growth rate constant due to CO2

    c_L = growthRateBA(L_l_val, LMAX);          // growth rate constant due to blowing agent

    while(counter < 2*nNodes)
    {
        sgBA[counter] = 0.0;
        sgCO2[counter] = 0.0;

        k = static_cast<double>(mOrder[counter]);

        if(counter == 0)
        {
            sgBA[counter] = 0.0;
            sgCO2[counter] = 0.0;
        }
        else if(counter == 1)
        {
            for(i = 0; i < nNodes; i++)
            {
                if(vi[i] > 0.0)
                {
                    sgBA[counter] += c_L*we[i];
                    sgCO2[counter] += c_CO2*we[i];
                }
                else
                {
                    sgBA[counter] = sgBA[counter];
                    sgCO2[counter] = sgCO2[counter];
                }
            }
        }
        else
        {
            for(i = 0; i < nNodes; i++)
            {
                if(vi[i] > 0.0)
                {
                    sgBA[counter] += k*c_L*we[i]*Foam::pow(vi[i], k - 1);
                    sgCO2[counter] += k*c_CO2*we[i]*Foam::pow(vi[i], k - 1);
                }
                else
                {
                    sgBA[counter] = sgBA[counter];
                    sgCO2[counter] = sgCO2[counter];
                }
            }
        }
        counter++;
    } // end of while
}

scalar Foam::Moments::growthRateConst()
{
// growthRateConst - constant growth rate
// @param -
    return 0.0955;
}

scalar Foam::Moments::growthRateCO2(scalar &CO2_l, scalar &wCO2_Max)
{
// growthRateCO2 - growth rate due to CO2
// @param - scalar &CO2_l - weight fraction of CO2 in liquid
// @param - scalar &wCO2_Max - weight fraction of maximum allowable CO2 in liquid
    scalar G0 = 1.0e-14; // m3s-1
    if (CO2_l > wCO2_Max)
    {
        return (G0*(CO2_l - wCO2_Max)/max(wCO2_Max,4.4e-4));
    }
    else
    {
        return (0.0);
    }
}

scalar Foam::Moments::growthRateBA(scalar &L_l, scalar &LMAX)
{
// growthRateBA - growth rate due to blowing agent
// @param - scalar &L_l - weight fraction of liquid blowing agent in liquid
// @param - scalar &LMAX - weight fraction of maximum allowable blowing agent in liquid
    scalar G0 = 1.0e-14;

    if (LMAX != 0.0 && L_l > LMAX)
    {
        return (G0*(L_l - LMAX)/LMAX);
    }
    else
    {
        return (0.0);
    }
}

void Foam::Moments::coalescenceSource
(
    double *sc, double *we, double *vi, int &nNodes, int *mOrder
)
{
    // coalescenceSource - source term for the bubbles coalescence
    // @param - double *sc - source of coalescence
    // @param - double *we - weights of quadrature
    // @param - duoble *vi - nodes of quadrature
    // @param - int &nNodes - number of nodes
    // @param - int *mOrder - order of moments

    int counter = 0;
    double k;
    int i, j;
    double coalescenceRate;

    // coalescence kernel beta0*(v(i) + v(j))
    // if counter = 0 the inside braket term = -1
    // if counter = 1 the inside braket = 0
    // if counter > 1 check vi[i]*vi[j] != 0 then do the summation. otherwise = zero
    // Multiply by 0.5 when adding to sg

    while(counter < 2*nNodes)
    {
        sc[counter] = 0.0;
        k = static_cast<double>(mOrder[counter]);

        coalescenceRate = coalescenceKernel();

        if(counter == 0)
        {
            for(i = 0; i < nNodes; i++)
            {
                for(j = 0; j < nNodes; j++)
                {
                    if (vi[i]*vi[j] != 0)
                    {
                        sc[counter] += coalescenceRate*(vi[i]+vi[j])*we[i]*we[j]*(-1.0);
                    }
                    else
                    {
                        sc[counter] = 0.0;
                    }
                }
            }
        }
        else if(counter == 1)
        {
            sc[counter] = 0.0;
        }
        else
        {
            for(i = 0; i < nNodes; i++)
            {
                for(j = 0; j < nNodes; j++)
                {
                    if(vi[i]*vi[j] != 0)
                    {
                        sc[counter] +=
                            (
                                coalescenceRate*(vi[i] + vi[j])*we[i]*we[j]
                                *(Foam::pow((vi[i] + vi[j]),k) -
                                Foam::pow(vi[i],k) - Foam::pow(vi[j],k))
                            );
                    }
                    else
                    {
                        sc[counter]  = 0.0;
                    }
                }
            }
        }
        counter++;
    } // End of while loop
}

scalar Foam::Moments::coalescenceKernel()
{
    // coalescenceKernel - constant kernel for coalescence
    return 0.0;
}
// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

void Foam::Moments::operator=(const Moments& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn("Foam::Moments::operator=(const Foam::Moments&)")
            << "Attempted assignment to self"
            << abort(FatalError);
    }
}

// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


// ************************************************************************* //
