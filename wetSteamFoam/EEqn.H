{
    volScalarField& he = steam.he();

    surfaceScalarField mshPhi(
        IOobject("meshPhi", mesh.time().timeName(), mesh),
        mesh, dimensionedScalar("0", dimVolume/dimTime, 0.0));

    MRF.makeAbsolute(mshPhi);

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + fvm::div(phi, he)
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
      + fvc::div(mshPhi,p)
      - fvm::laplacian(turbulence->alphaEff(), he)
     ==
      - fvc::div(turbulence->devRhoReff() & U)
      + fvOptions(rho, he)
    );

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);

    steam.correct();
}
