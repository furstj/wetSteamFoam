/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "condensationMonodispersion.H"
#include "addToRunTimeSelectionTable.H"
#include "fvm.H"
#include "physicoChemicalConstants.H"
#include "fundamentalConstants.H"
#include "mathematicalConstants.H"
#include "turbulentWetSteamModel.H"

namespace Foam
{
namespace WetSteam
{

defineTypeNameAndDebug(condensationMonodispersion, 0);
addToRunTimeSelectionTable(condensationModel, condensationMonodispersion, params);

condensationMonodispersion::condensationMonodispersion(
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi,
    wetSteamSystem& steam
):
    condensationModel(rho, U, phi, steam),

    Q0_(
        IOobject
        (
            "Q0",
            U.mesh().time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),

    mDotGL_(
        IOobject
        (
            "mDotGL",
            U.mesh().time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("zero", dimDensity/dimTime, 0.0)
    ),

    m1_(
        "m1",
        dimMass,
        steam.gas().lookupOrDefault<scalar>("molecularMass", 2.99046e-26)
    ),

    beta_(
	"beta",
	dimless,
	steam.gas().lookupOrDefault<scalar>("surfaceTensionCorrection", 1)
    ),

    Sct_(
        "Sct",
        dimless,
        steam.gas().lookupOrDefault<scalar>("SchmidtNumber", 0.9)
    )
{
}


void condensationMonodispersion::correct()
{
    const fvMesh& mesh = steam().mesh();
    volScalarField& w = steam().w();
    const volScalarField Ts = steam().TSat();
    const volScalarField ps = steam().pSat();
    const volScalarField T = steam().T();
    const volScalarField p = steam().p();
    const volScalarField sigma = steam().liquid().sigma();
    const volScalarField rho_l = steam().liquid().rho();
    const volScalarField rho_g = steam().gas().rho();
    const volScalarField L = steam().L();

    const WetSteam::turbulenceModel& turbModel =
        mesh.lookupObject<WetSteam::turbulenceModel>
        (
            IOobject::groupName
            (
                turbulenceModel::propertiesName,
                p.group()
            )
        );

    const scalar pi = constant::mathematical::pi;
    const scalar kB = constant::physicoChemical::k.value();
    const scalar M  = constant::physicoChemical::NA.value()*m1_.value();
    const scalar Rg = constant::physicoChemical::R.value()/M;
    
    const scalar rMin = 1e-12;
    const dimensionedScalar wMin("wMin", dimless, 1.e-6);

    // Droplet nucleation rate [m^-3 s^-1]
    volScalarField J(
        IOobject("J", mesh.time()),
        mesh,
        dimensionedScalar("zero", dimless/dimTime/dimVolume, 0.0)
    );

    // Droplet growth rate [kg/m3/s]
    volScalarField mDot(
        IOobject("mDot", mesh.time()),
        mesh,
        dimensionedScalar("zero", dimDensity/dimTime, 0.0)
    );

    // Droplet critical radius [m]
    volScalarField rc(
        IOobject("rc", mesh.time()),
        mesh,
        dimensionedScalar("zero", dimLength, 0.0)
    );

    forAll(J, i)
    {

        // Knudsen number
        scalar Kn = 0;
        
        // Average droplet radius
        scalar r = 0;
        if (w[i] > wMin.value() && Q0_[i] > 0)
        {
            r = pow(3*w[i]/(4*Q0_[i]*pi*rho_l[i]), 1.0/3.0);
            scalar eta = 1.823e-6*sqrt(T[i])/(1 + 673/T[i]);
            Kn = eta*sqrt(2*pi*Rg*T[i])/(4*r*p[i]);
        }

        scalar lambda_g = 7.341e-3 - 1.013e-5*T[i] + 1.801e-7*sqr(T[i]) - 9.1e-11*pow3(T[i]);
        scalar rho = rho_g[i]/(1 - w[i]);
        
        if (T[i] >= Ts[i])   // Evaporation
        {
            if (r <= rMin) // Complete evaporation
            {
                w[i] = 0;
                Q0_[i] = 0;
		J[i] = 0;
		mDot[i] = 0;
            }
            else // No new droplets, only evaporation
            {
                J[i] = 0;
                scalar rDot = lambda_g*(Ts[i] - T[i])/(L[i]*rho_l[i]*(1 + 3.18*Kn))/r;
                mDot[i] = 4*pi*rho*sqr(r)*rDot*rho_l[i]*Q0_[i];
            }
        }
        else  // Condensation (T < Ts)
        {
#ifdef HALAMA
            rc[i] = 2*sigma[i]/(L[i]*rho_l[i]*log(Ts[i]/T[i]));
#else
            scalar S = p[i]/ps[i];
            rc[i] = 2*sigma[i]/(rho_l[i]*Rg*T[i]*log(S));
#endif            
            J[i] = sqrt(2*sigma[i]/(pi*pow3(m1_.value())))*sqr(rho_g[i])/rho_l[i]*
                exp(-beta_.value()*4*pi*sqr(rc[i])*sigma[i]/(3*kB*T[i]));

           if (r > rc[i])
           {
 #ifdef HALAMA
               scalar theta = Ts[i]/T[i] - 1;
               scalar dTbyLog = T[i]*(1 + theta/2 - sqr(theta)/12 + pow3(theta)/24 
               - 19*pow4(theta)/720 + 3*pow5(theta)/160);
               scalar rDot = lambda_g/(L[i]*rho_l[i]*(1 + 3.18*Kn))*
                   ((Ts[i] - T[i])/r - 2*sigma[i]/(L[i]*rho_l[i]*sqr(r))*dTbyLog);
#else
               scalar rDot = lambda_g*(Ts[i] - T[i])/(L[i]*rho_l[i]*(1 + 3.18*Kn)) *
                   (r - rc[i])/sqr(r);
#endif
               mDot[i] = 4*pi*rho*sqr(r)*rDot*rho_l[i]*Q0_[i];
           }

        }

    }

    mDotGL_ = 4.0/3.0*pi*pow3(rc)*J*rho_l + mDot;

    volScalarField Dt("Dt", rho_*turbModel.nut()/Sct_);
    
    {
        fvScalarMatrix Q0Eqn
        (
            fvm::ddt(rho_, Q0_) + fvm::div(phi_, Q0_) - fvm::laplacian(Dt, Q0_)
            ==
            J
        );

	Q0Eqn.relax();
        Q0Eqn.solve();
    }

    {
        fvScalarMatrix wEqn
        (
            fvm::ddt(rho_, w) + fvm::div(phi_, w) - fvm::laplacian(Dt, w)
            ==
            4.0/3.0*pi*pow3(rc)*J*rho_l - fvm::SuSp(-mDot/(w+wMin), w)
        );
    
        wEqn.relax();
        wEqn.solve();
    }
    

    // Explicitly constrain w & Q0 in dry steam region
    forAll(w, i)
    {
        if (w[i] < wMin.value() && J[i] <= 0)
	{
	    w[i] = 0;
	    Q0_[i] = 0;
	}
	w[i] = max(w[i],0);
	Q0_[i] = max(Q0_[i],0);
    }
}


}
}

// ************************************************************************* //
