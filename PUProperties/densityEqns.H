dimensionedScalar oneTwo
(
    "oneTwo",
    dimDensity,
    1.0
);

rho_gas = 1.2*oneTwo;

if (PBESwitch)
{
    foamDensity foam_density
                (
                    mesh, runTime, alpha2, mOne, wCO2_g, wBA_g, wBA_l, p,
                    TS, XW, L0, rhoPoly
                );

    forAll(mesh.C(), celli)
    {
        if (wCO2_g[celli] < 0.0)
        {
            wCO2_g[celli] = ROOTVSMALL;
        }
        if (wBA_g[celli] < 0.0)
        {
            wBA_g[celli] = ROOTVSMALL;
        }
        if (TS[celli] < 0.0)
        {
            TS[celli] = 298.0;
        }
    }

    rho_foam = foam_density.rhoFoam(M_CO2, M_B);
}
else
{
    foamDensity foam_density
                (
                    mesh, runTime, alpha2, mOne, wCO2_g, wBA_g, wBA_l, p, 
                    TS, XW, L0, rhoPoly
                );

    rho_foam =
        foam_density.rhoFoamNoPBE
        (
            M_CO2, M_B, surfaceTension, M_liq, CW_0, rhoBL
        );
}