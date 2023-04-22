#if (WM_PROJECT_VERSION==6 || OPENFOAM >= 1912)
CorrectPhi
(
    U,
    phi,
    p,
    rho,
    psi,
    dimensionedScalar("rAUf", dimTime, 1),
    divrhoU(),
    pimple //,
    //true
);
#else
CorrectPhi
(
    U,
    phi,
    p,
    rho,
    psi,
    dimensionedScalar("rAUf", dimTime, 1),
    divrhoU(),
    pimple,
    true
);
#endif
