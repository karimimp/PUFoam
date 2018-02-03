/*
Checking if the cell values of wBA_l, wBA_g, wCO2_l, and wCO2_g are bounded.
*/
    forAll(mesh.C(), celli)
    {
        if (wBA_l[celli] > L0)
        {
            wBA_l[celli] = L0;
        }
        if (wBA_l[celli] < 0.0)
        {
            wBA_l[celli] = ROOTVSMALL;
        }
        if (wBA_g[celli] < 0.0)
        {
            wBA_g[celli] = ROOTVSMALL;
        }
        if (wCO2_g[celli] < 0.0)
        {
            wCO2_g[celli] = ROOTVSMALL;
        }
    }
