/* Source terms for moments transport equations */

if (PBESwitch)
{
    Moments moments;
    if (PBEMethod == "QMOM" && nNodes == 2)
    {
        static label nodes = nNodes;
        label mOrder[2*nodes];
        for (label i = 0; i < 2*nNodes; i++){mOrder[i] = i;}

        scalar mom[2*nodes], we[nodes], vi[nodes], sgBA[2*nodes],
               sgCO2[2*nodes], sc[2*nodes], sources[2*nodes];
        scalar LMAX, xBL_val;

        forAll(mesh.C(), celli)
        {
            // Using old moments to compute the new weights and nodes
            if (alpha2[celli] > 0.5)
            {
                mom[0] = mZero[celli];
                mom[1] = mOne[celli];
                mom[2] = mTwo[celli];
                mom[3] = mThree[celli];

                // calling the functions, PDA, growthSource, coalescenceSource
                moments.PDA(we, vi, mom, nodes);

                weight0[celli] = we[0];
                weight1[celli] = we[1];
                node0[celli] = vi[0];
                node1[celli] = vi[1];

                // Required variables for wCO2_Max
                scalar henry_Coefficient, bubble_radius,
                       partialPressure_CO2, wCO2_Max;
                henry_Coefficient = PUProperties.henryCoefficient(TS[celli]);
                bubble_radius = PUProperties.bubbleRadius(mom[0], mom[1]);
                partialPressure_CO2 =
                    PUProperties.partialPressureCO2
                    (
                        M_CO2, M_B, surfaceTension,
                        wCO2_g[celli], wBA_g[celli], p[celli], bubble_radius
                    );
                wCO2_Max =
                    PUProperties.wCO2Max
                    (
                        M_CO2, M_liq, partialPressure_CO2, henry_Coefficient
                    );

                if (blowingAgent == "n-pentane")
                {
                    LMAX = PUProperties.LliqMax (TS[celli]);
                }
                if (blowingAgent == "R-11")
                {
                    xBL_val = PUProperties.xBL(TS[celli], dxdT);
                    LMAX = PUProperties.wBL_D(xBL_val, M_B, M_NCO);
                }
                if (blowingAgent == "no")
                {
                    LMAX = 0.0;
                }
                moments.growthSource(sgBA, sgCO2, we, vi, nodes, mOrder,
                        wCO2_l[celli], wBA_l[celli], wCO2_Max, LMAX);

                g1_CO2[celli] = sgCO2[1];
                g1_BA[celli] = sgBA[1];
                cc1[celli] = moments.cc1Value(wBA_l[celli], LMAX);

                if (XOH[celli] > XOH_Gel)
                {
                    for(label j = 0; j < 2*nodes; j++)
                    {
                        sc[j] = 0.0;
                    }
                }
                else
                {
                    moments.coalescenceSource(sc, we, vi, nodes, mOrder);
                }

                for(label j = 0; j < 2*nodes; j++)
                {
                    sources[j] = sgBA[j] + sgCO2[j] + 0.5*sc[j];
                }
                mZeroSource[celli] = sources[0];
                mOneSource[celli] = sources[1];
                mTwoSource[celli] = sources[2];
                mThreeSource[celli] = sources[3];
            }
            else
            {
                mZeroSource[celli] = ROOTVSMALL;
                mOneSource[celli] = ROOTVSMALL;
                mTwoSource[celli] = ROOTVSMALL;
                mThreeSource[celli] = ROOTVSMALL;
            }
        } //end of forAll
    }

    if ( PBEMethod == "QMOM" && nNodes == 3 )
    {
        static label nodes = nNodes;
        label mOrder[2*nodes];
        for (label i=0;i<2*nNodes;i++){mOrder[i] = i;}
        scalar mom[2*nodes], we[nodes], vi[nodes], sgBA[2*nodes],
                sgCO2[2*nodes], sc[2*nodes], sources[2*nodes];

        scalar LMAX, xBL_val;

        forAll(mesh.C(), celli)
        {
            if(alpha2[celli] > 0.5)
            {
                mom[0] = mZero[celli];
                mom[1] = mOne[celli];
                mom[2] = mTwo[celli];
                mom[3] = mThree[celli];
                mom[4] = mFour[celli];
                mom[5] = mFive[celli];

                // calling the functions, PDA, growthSource, coalescenceSource
                moments.PDA(we, vi, mom, nodes);
                weight0[celli] = we[0];
                weight1[celli] = we[1];
                weight2[celli] = we[2];
                node0[celli] = vi[0];
                node1[celli] = vi[1];
                node2[celli] = vi[2];

                // Required variables for wCO2_Max
                scalar henry_Coefficient, bubble_radius,
                       partialPressure_CO2, wCO2_Max;
                henry_Coefficient = PUProperties.henryCoefficient(TS[celli]);
                bubble_radius = PUProperties.bubbleRadius(mom[0], mom[1]);
                partialPressure_CO2 =
                    PUProperties.partialPressureCO2
                    (
                        M_CO2, M_B, surfaceTension, wCO2_g[celli],
                        wBA_g[celli], p[celli], bubble_radius
                    );
                wCO2_Max =
                    PUProperties.wCO2Max
                    (
                        M_CO2, M_liq, partialPressure_CO2, henry_Coefficient
                    );

                if (blowingAgent == "n-pentane")
                {
                    LMAX = PUProperties.LliqMax (TS[celli]);
                }
                if (blowingAgent == "R-11")
                {
                    xBL_val = PUProperties.xBL(TS[celli], dxdT);
                    LMAX = PUProperties.wBL_D(xBL_val, M_B, M_NCO);
                }
                if (blowingAgent == "no")
                {
                    LMAX = 0.0;
                }
                moments.growthSource(sgBA, sgCO2, we, vi, nodes, mOrder,
                                wCO2_l[celli], wBA_l[celli], wCO2_Max, LMAX);
                g1_CO2[celli] = sgCO2[1];
                g1_BA[celli] = sgBA[1];
                cc1[celli] = moments.cc1Value(wBA_l[celli], LMAX);

                if (XOH[celli] > XOH_Gel)
                {
                    for(int j = 0; j < 2*nodes; j++)
                    {
                        sc[j] = 0.0;
                    }
                }
                else
                {
                    moments.coalescenceSource(sc, we, vi, nodes, mOrder);
                }

                for(int j = 0; j < 2*nodes; j++)
                {
                    sources[j] = sgBA[j] + sgCO2[j] + 0.5*sc[j];
                }
                mZeroSource[celli] = sources[0];
                mOneSource[celli] = sources[1];
                mTwoSource[celli] = sources[2];
                mThreeSource[celli] = sources[3];
                mFourSource[celli] = sources[4];
                mFiveSource[celli] = sources[5];
            }
            else
            {
                mZeroSource[celli] = ROOTVSMALL;
                mOneSource[celli] = ROOTVSMALL;
                mTwoSource[celli] = ROOTVSMALL;
                mThreeSource[celli] = ROOTVSMALL;
                mFourSource[celli] = ROOTVSMALL;
                mFiveSource[celli] = ROOTVSMALL;
            }
        } //end of forAll
    }
}