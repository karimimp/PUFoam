rheologyPU foamViscosity(mesh, runTime, alpha2, TS, XW, XOH, U, XOH_Gel);
// Constant viscosity (Newtonian)
if (viscosityModel == "constant")
{
	muFoamCorr = foamViscosity.waterLike();
}

// Castro-Macosko model (Newtonian)
if (viscosityModel == "castro-macosko")
{
	muFoamCorr = foamViscosity.CastroMacosko(COH_0, CNCO_0, CW_0);
}

// Modified Bird-Carreau model (non-Newtonian)
if (viscosityModel == "bird-carreau")
{
	muFoamCorr = foamViscosity.BirdCarreau(mu0, muinf, muFoam);
}

muMixture 	= alpha1*muAir + alpha2*muFoamCorr;
