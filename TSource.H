    scalar dXOHdt, dWdt, dLdt, C_TOT, xBL_val;

    forAll(mesh.C(), celli)
    {
        dXOHdt = PUProperties.dQdt(XOH[celli], XOH.oldTime()[celli]);
        dWdt = PUProperties.dQdt(XW[celli], XW.oldTime()[celli]);
        if (PBESwitch)
        {
            dLdt = PUProperties.dQdt(wBA_l[celli], wBA_l.oldTime().oldTime()[celli]);
        }
        else
        {
            if (PUProperties.LliqMax(TS[celli]) < L0)
            {
                dLdt = -1.0*PUProperties.ddT_LliqMax(TS[celli])*
                        PUProperties.dQdt(TS[celli], TS.oldTime()[celli]);
            }
            else
            {
                dLdt = 0.0;
            }
        }

        if (blowingAgent == "n-pentane")
        {
            C_TOT = C_Poly;
            if (alpha2[celli] > 0.5 && XOH[celli] < XOH_Gel)
            {
                TSSource[celli] =
                    PUProperties.exothermicGelling
                    (
                        DH_OH, COH_0, dXOHdt, rhoPoly, C_TOT
                    )
                  + PUProperties.exothermicBlowing
                    (
                        DH_W, CW_0, dWdt, rhoPoly, C_TOT
                    )
                  + PUProperties.endothermicEvaporation
                    (
                        latenth, dLdt, rhoPoly, C_TOT
                    );
            }
            else if (alpha2[celli] > 0.5 && XOH[celli] >= XOH_Gel)
            {
                TSSource[celli] =
                    PUProperties.endothermicEvaporation
                    (
                        latenth, dLdt, rhoPoly, C_TOT
                    );
            }
            else
            {
                TSSource[celli] = ROOTVSMALL;
            }
        }

        if (blowingAgent == "R-11")
        {
            if (wBA_g[celli] > 0.0)
            {
                xBL_val = PUProperties.xBL(TS[celli], dxdT);

                C_TOT = C_Poly + PUProperties.R11EvaporationRate
                                (
                                    M_B, M_NCO, xBL_val, dxdT, latenth
                                );

                if (alpha2[celli] > 0.5 && XOH[celli] < XOH_Gel)
                {
                    TSSource[celli] =
                        PUProperties.exothermicGelling
                        (
                            DH_OH, COH_0, dXOHdt, rhoPoly, C_TOT
                        )
                      + PUProperties.exothermicBlowing
                        (
                            DH_W, CW_0, dWdt, rhoPoly, C_TOT
                        );
                }
                else if (alpha2[celli] > 0.5 && XOH[celli] >= XOH_Gel)
                {
                    TSSource[celli] =
                        PUProperties.R11EvaporationRate
                        (
                            M_B, M_NCO, xBL_val, dxdT, latenth
                        );
                }
                else
                {
                    TSSource[celli] = ROOTVSMALL;
                }
            }
            else
            {
                C_TOT   = C_Poly;
                if (alpha2[celli] > 0.5 && XOH[celli] < XOH_Gel)
                {
                    TSSource[celli] =
                        PUProperties.exothermicGelling
                        (
                            DH_OH, COH_0, dXOHdt, rhoPoly, C_TOT
                        )
                      + PUProperties.exothermicBlowing
                        (
                            DH_W, CW_0, dWdt, rhoPoly, C_TOT
                        );
                }
                else
                {
                    TSSource[celli]  = ROOTVSMALL;
                }
            }
        }
        if (blowingAgent == "no")
        {
            C_TOT       = C_Poly;
            if (alpha2[celli] > 0.5 && XOH[celli] < XOH_Gel)
            {
                TSSource[celli] =
                    PUProperties.exothermicGelling
                    (
                        DH_OH, COH_0, dXOHdt, rhoPoly, C_TOT
                    )
                  + PUProperties.exothermicBlowing
                    (
                        DH_W, CW_0, dWdt, rhoPoly, C_TOT
                    );
            }
            else
            {
                TSSource[celli] = ROOTVSMALL;
            }
        }
        if (blowingAgent == "n-pentane")
        {
            if (alpha2[celli] > 0.5)
            {
                if(rho_foam[celli] > 48)
                {
                    thermalConductivity[celli] =
                        PUProperties.npentaneThermalconductivityHighrho
                        (rho_foam[celli]);
                }
                else
                {
                    thermalConductivity[celli] =
                        PUProperties.npentaneThermalconductivityLowrho
                        (rho_foam[celli]);
                }
                thermalDiff_foam[celli] =
                    (thermalConductivity[celli]/(rho_foam[celli]*C_Poly));
            }
            else
            {
                thermalDiff_foam[celli] =
                    PUProperties.thermalDiffusivityGas(TS[celli]);
            }

            thermalDiff_gas[celli] = PUProperties.thermalDiffusivityGas(TS[celli]);
        }
        else
        {
            if (alpha2[celli] > 0.5)
            {
                thermalConductivity[celli] = 0.03;
                thermalDiff_foam[celli] =
                    (thermalConductivity[celli]/(rho_foam[celli]*C_Poly));
            }
            else
            {
                thermalDiff_foam[celli] =
                    PUProperties.thermalDiffusivityGas(TS[celli]);
            }
            thermalDiff_gas[celli] = PUProperties.thermalDiffusivityGas(TS[celli]);
        }
        thermalDiffusivity[celli] = alpha1[celli]*thermalDiff_gas[celli] +
                                    alpha2[celli]*thermalDiff_foam[celli];
}