volTensorField gradU = fvc::grad(U);
  fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(rhoPhi, U) 
     	- fvm::laplacian(muMixture, U) 
    );

    UEqn.relax();
    //*momentumPredictor : flag to indicate to solve for momentum

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    fvc::interpolate(interface.sigmaK())*fvc::snGrad(alpha1)
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );

        K = 0.5*magSqr(U);
    }