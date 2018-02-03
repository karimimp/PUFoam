/* Description:
Source terms for blowing agent and CO2 in both liquid and gas
*/
if (PBESwitch)
{
    // Kinetics dissolvedBA(mesh, runTime, alpha2, TS);
    blowingAgents sourceTerm(mesh, runTime, alpha2, TS);
    PUgeneric PUProperties(mesh, runTime, alpha2, TS);

    forAll(mesh.C(), celli)
    {
        if (blowingAgent == "n-pentane")
        {
            if (alpha2[celli] > 0.5)
            {
                if (PUProperties.LliqMax(TS[celli]) < L0)
                {
                    wBA_gSource[celli] = (-cc1[celli]*
                                          sourceTerm.n_pentaneInGas(sourceTerm,
                                            TS[celli], TS.oldTime()[celli])
                                         );
                    wBA_lSource[celli] =
                        sourceTerm.n_pentaneInLiquid
                        (M_B, rhoPoly, p[celli], RR,TS[celli])*(-g1_BA[celli]);
                }
                else
                {
                    wBA_gSource[celli] = 0.0;
                    wBA_lSource[celli] = 0.0;
                }
            }
            else
            {
                wBA_lSource[celli] = 0.0;
                wBA_gSource[celli] = 0.0;
            }
        }

        if (blowingAgent == "R-11")
        {
            scalar xBL_val, wBL_Max;
            xBL_val = PUProperties.xBL(TS[celli],dxdT);
            wBL_Max = PUProperties.wBL_D(xBL_val, M_B, M_NCO);

            if (wBL_Max > 1.0e-4 && alpha2[celli] > 0.5)
            {
                wBA_gSource[celli] = (
                                      cc1[celli]*
                                      sourceTerm.R11InGas
                                      (
                                        sourceTerm, M_B, M_NCO,
                                        dxdT, TS[celli], TS.oldTime()[celli])
                                     );
                wBA_lSource[celli] = (
                                      -cc1[celli]*
                                      sourceTerm.R11InGas
                                      (
                                        sourceTerm, M_B, M_NCO,
                                        dxdT, TS[celli], TS.oldTime()[celli])
                                     );
            }
            else
            {
                wBA_gSource[celli] = 0.0;
                wBA_lSource[celli] = 0.0;
            }
        }

        if (CW_0 == 0.0)
        {
            wCO2_lSource[celli] = 0.0;
        }
        else
        {
            if ( XOH[celli] < XOH_Gel )
            {
                wCO2_lSource[celli] =
                (
                    sourceTerm.CO2InLiquid(sourceTerm, M_CO2, CW_0, rhoPoly,
                    XW[celli], XW.oldTime()[celli]) - g1_CO2[celli]*
                    sourceTerm.CO2InGas(M_CO2, rhoPoly, p[celli], RR, TS[celli])
                );
            }

        }
        scalar henry_Coefficient, bubble_radius, partialPressure_CO2, wCO2_Max;
        henry_Coefficient = PUProperties.henryCoefficient(TS[celli]);
        bubble_radius = PUProperties.bubbleRadius(mZero[celli], mOne[celli]);
        partialPressure_CO2 =
            PUProperties.partialPressureCO2
                (
                    M_CO2, M_B, surfaceTension, wCO2_g[celli], wBA_g[celli],
                    p[celli], bubble_radius
                );
        wCO2_Max =
            PUProperties.wCO2Max
                (M_CO2, M_liq, partialPressure_CO2, henry_Coefficient);

        if(wCO2_l[celli] < wCO2_Max || CW_0 == 0.0 || XOH[celli] > XOH_Gel)
        {
            wCO2_gSource[celli] = 0.0;
        }
        else
        {
            wCO2_gSource[celli] =
                (
                    g1_CO2[celli]*
                    sourceTerm.CO2InGas(M_CO2, rhoPoly, p[celli], RR, TS[celli])
                );
        }
    }
}
else
{
    blowingAgents sourceTerm(mesh, runTime, alpha2, TS);
    // Kinetics dissolvedBA(mesh, runTime, alpha2, TS);
    PUgeneric PUProperties(mesh, runTime, alpha2, TS);

    forAll(mesh.C(), celli)
    {
        if (blowingAgent == "n-pentane")
        {
            if (alpha2[celli] > 0.5)
            {
                if (L0 < PUProperties.LliqMax(TS[celli]))
                {
                    wBA_gSource[celli] = scalar(-1.0)*
                                         sourceTerm.n_pentaneInGas(sourceTerm,
                                            TS[celli], TS.oldTime()[celli]);
                    wBA_lSource[celli]  = sourceTerm.n_pentaneInGas(sourceTerm,
                                            TS[celli], TS.oldTime()[celli]);
                }
                else
                {
                    wBA_gSource[celli]  = 0.0;
                    wBA_lSource[celli]  = 0.0;
                }
            }
            else
            {
                wBA_lSource[celli]  = 0.0;
                wBA_gSource[celli]  = 0.0;
            }
        }

        if (blowingAgent == "R-11")
        {
            scalar xBL0;
            xBL0 = ((M_NCO/M_B)*L0)/(1.0 + (M_NCO/M_B)*L0);
            creamT[celli] = PUProperties.creamTemperature(xBL0, dxdT);
            if (alpha2[celli] > 0.5)
            {
                if
                (
                    TS[celli] > creamT[celli]
                    && wBA_l[celli] > 1.0e-8
                    && TS[celli] < 383.0
                )
                {
                    wBA_gSource[celli]  =
                        sourceTerm.R11InGas
                        (
                            sourceTerm, M_B, M_NCO,
                            dxdT, TS[celli], TS.oldTime()[celli]
                        );
                    wBA_lSource[celli]  =
                        scalar(-1.0)*
                        sourceTerm.R11InGas
                        (
                            sourceTerm, M_B, M_NCO,
                            dxdT, TS[celli], TS.oldTime()[celli]
                        );
                }
                else
                {
                    wBA_gSource[celli]  = 0.0;
                    wBA_lSource[celli]  = 0.0;
                }
            }
            else
            {
                wBA_lSource[celli]  = 0.0;
                wBA_gSource[celli]  = 0.0;
            }

        }

        wCO2_lSource[celli] = scalar(0.0);

        scalar henry_Coefficient, bubble_radius, partialPressure_CO2, wCO2_Max;
        henry_Coefficient = PUProperties.henryCoefficient(TS[celli]);
        bubble_radius = PUProperties.bubbleRadius(mZero[celli], mOne[celli]);
        partialPressure_CO2 =
            PUProperties.partialPressureCO2
                (
                    M_CO2, M_B, surfaceTension, wCO2_g[celli], wBA_g[celli],
                    p[celli], bubble_radius
                );
        wCO2_Max =
            PUProperties.wCO2Max
                (M_CO2, M_liq, partialPressure_CO2, henry_Coefficient);

        if (alpha2[celli] > 0.5 && XOH[celli] <= XOH_Gel && CW_0 != 0.0)
        {
                        blowingAgents sourceTerm(mesh, runTime, alpha2, TS);

            wCO2_gSource[celli] =
                (
                sourceTerm.CO2InGasNoPBE(CW_0, M_CO2, rhoPoly, XW[celli])
                - wCO2_Max
                );
        }
        else
        {
            wCO2_gSource[celli] = 0.0;
        }

        if (wCO2_gSource[celli] < 0.0)
        {
            wCO2_gSource[celli] = scalar(0.0);
        }
    }
}