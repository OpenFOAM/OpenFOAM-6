{
    volScalarField tempT = T; 
    volScalarField& he = thermo.he();

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + mvConvection->fvmDiv(phi, he)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      + (
            he.name() == "e"
          ? fvc::div
            (
                fvc::absolute(phi/fvc::interpolate(rho), U),
                p,
                "div(phiv,p)"
            )
          : -dpdt
        )
      - fvm::laplacian(turbulence->alphaEff(), he)
     ==
        Qdot
      + fvOptions(rho, he)
    );
    
    forAll(Y, k)
    {
        EEqn -= fvc::laplacian(turbulence->alphaEff()*hsi[k], Y[k]);
        EEqn -= fvc::div(J[k], hsi[k], "div(Ji,hsi)");
    }

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);

    thermo.correct();

    Info << "min/max(T) = "
         << min(T).value() << ", " << max(T).value() << endl;
    // update oldT for CoDAC
    //oldT = tempT;
}


