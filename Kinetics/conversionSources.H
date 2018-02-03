
    // Kinetics kinetics(mesh, runTime, alpha2, TS);
    blowingReaction sourceTermBlowing(mesh, XW, TS);
    gellingReaction sourceTermGelling(mesh, XOH, TS);
    PUgeneric PUProperties(mesh, runTime, alpha2, TS);
    forAll(mesh.C(), celli)
    {
        // XW Source term
        if (XW[celli] < 0.0)
        {
            XW[celli] = ROOTVSMALL;
        }
        if (XW[celli] > 1.0)
        {
            XW[celli] = ROOTVSMALL;
        }

        if
        (
            (sourceTermBlowing.QKinW(AW, EW, TS[celli], L0, rhoPoly, rhoBL) < 0.0)
            || CW_0 == 0.0
        )
        {
            XWSource[celli] = ROOTVSMALL;
        }

        if (XOH[celli] < XOH_Gel)
        {
            XWSource[celli] =
                sourceTermBlowing.QKinW(AW, EW, TS[celli], L0, rhoPoly, rhoBL);
        }

        // XOH Source term
        if (XW[celli] < 0.0)
        {
            XW[celli] = ROOTVSMALL;
        }
        else if (XW[celli] > 1.0)
        {
            XW[celli] = ROOTVSMALL;
        }
        else
        {
            XW[celli] = XW[celli];
        }

        if (XOH[celli] < 0.0)
        {
            XOH[celli] = ROOTVSMALL;
        }
        if (XOH[celli] > 1.0)
        {
            XOH[celli] = ROOTVSMALL;
        }

        if (XOH[celli] < XOH_Gel)
        {
            if
            (
                (PUProperties.arrhenius(AOH, EOH, TS.oldTime()[celli])*
                sourceTermGelling.gellingSourceOldTime(COH_0, CNCO_0, CW_0,
                                                       XW[celli], XOH[celli])*
                PUProperties.BAdilution(L0, rhoPoly, rhoBL)) < 0.0
            )
            {
                ScNonXOH[celli] = ROOTVSMALL;
                Coeff_1stOrd[celli] = ROOTVSMALL;
                SuXOH[celli] = ROOTVSMALL;
                SpXOH[celli] = ROOTVSMALL;
            }
            else
            {
                ScNonXOH[celli] = PUProperties.BAdilution(L0, rhoPoly, rhoBL)*
                                  (PUProperties.arrhenius(AOH, EOH, TS[celli])*CNCO_0
                                  - 2.0*PUProperties.arrhenius(AOH, EOH, TS[celli])*
                                  CW_0*XW[celli]
                                  );
                Coeff_1stOrd[celli] = PUProperties.BAdilution(L0, rhoPoly, rhoBL)*
                                      (PUProperties.arrhenius(AOH, EOH, TS[celli])*
                                       COH_0 +
                                       PUProperties.arrhenius(AOH, EOH, TS[celli])*
                                       CNCO_0 -
                                       2.0*PUProperties.arrhenius(AOH, EOH, TS[celli])*
                                       CW_0*XW[celli]
                                      );
                SuXOH[celli] = scalar(-1.0)*
                               PUProperties.BAdilution(L0, rhoPoly, rhoBL)*
                               PUProperties.arrhenius(AOH, EOH, TS.oldTime()[celli])*
                               COH_0*Foam::pow(XOH[celli],2.0);
                SpXOH[celli] = scalar(2.0)*
                               PUProperties.BAdilution(L0, rhoPoly, rhoBL)*
                               PUProperties.arrhenius(AOH, EOH, TS.oldTime()[celli])*
                               COH_0*XOH[celli];
            }
        }
    }