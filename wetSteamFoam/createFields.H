#include "createRDeltaT.H"

Info<< "Reading thermophysical properties\n" << endl;

WetSteam::wetSteamSystem steam(mesh);

fluidThermo& thermo = steam.gas();
thermo.validate(args.executable(), "h");


volScalarField& p = steam.p();


volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    steam.rho()
);


Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


#include "compressibleCreatePhi.H"

pressureControl pressureControl(p, rho, pimple.dict(), false);

mesh.setFluxRequired(p.name());


Info<< "Creating turbulence model\n" << endl;
autoPtr<WetSteam::turbulenceModel> turbulence
(
    WetSteam::turbulenceModel::New
    (
        rho,
        U,
        phi,
        steam
    )
);


Info<< "Creating condensation model\n" << endl;
autoPtr<WetSteam::condensationModel> condensation
(
    WetSteam::condensationModel::New(rho, U, phi, steam)
);


Info<< "Creating field dpdt\n" << endl;
volScalarField dpdt
(
    IOobject
    (
        "dpdt",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("dpdt", p.dimensions()/dimTime, 0)
);

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));

#include "createMRF.H"
#include "createFvOptions.H"
