/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "wetSteamSystem.H"

namespace Foam
{
namespace WetSteam
{



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

wetSteamSystem::wetSteamSystem
(
    const fvMesh& mesh,
    const word& phaseName
)
    :
    
    mesh_(mesh),
    
    gas_(fluidThermo::New(mesh)),

    gasProps_(gasProperties::New(gas_())),
    
    liquid_(*this, gas_()),
    
    w_(
        IOobject
        (
            "w",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    
    he_(
        IOobject
        (
            "h",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimEnergy/dimMass,
        replaceFixedEnergyTypes(gas_->he().boundaryField().types())
    ),
   
    psi_(
	"psi_s",
        gas_->psi()/(1-w_)
    ),

    saturation_(saturationCurve::New(gas_()))
{
    volScalarField Q = w_*L();
    setField(he_, gas().he() - Q);
}


tmp<volScalarField> wetSteamSystem::rho() const
{
    return gas().rho()/(1.0 - w_); 
}

tmp<scalarField> wetSteamSystem::rho(label patchi) const
{
    const scalarField& wp = w().boundaryField()[patchi];
    return gas().rho(patchi)/(1.0 - wp); 
}

tmp<volScalarField> wetSteamSystem::mu() const
{
    tmp<volScalarField> mug = gas().mu();
    if (mug()[0] > 0)
    {
        return (1 - w())*mug + w()*liquid().mu();
    }
    else
    {
        return mug;
    }
}

tmp<scalarField> wetSteamSystem::mu(const label patchi) const
{
    tmp<scalarField> mug = gas().mu(patchi)().clone();
    if (mug().size() > 0 && mug()[0] > 0)
    {
        const scalarField& pp = p().boundaryField()[patchi];
        const scalarField& Tp = T().boundaryField()[patchi];
        const scalarField& wp = w().boundaryField()[patchi];
        scalarField& mup = mug.ref(); 
        forAll(mup, i)
        {
            mup[i] = (1 - wp[i])*mup[i] + wp[i]*liquid().properties().mu(pp[i],Tp[i]);
        }
        return mug;
    }
    else
    {
        return mug;
    }
}


tmp<volScalarField> wetSteamSystem::nu() const
{
    return mu()/rho();
}

tmp<scalarField> wetSteamSystem::nu(label patchi) const
{
    return mu(patchi)/rho(patchi);
}


tmp<volScalarField> wetSteamSystem::alpha() const
{
    tmp<volScalarField> alphag = gas().alpha();
    if (alphag()[0] > 0)
    {
        return (1 - w())*alphag + w()*liquid().alpha();
    }
    else
    {
        return alphag;
    }
}

tmp<scalarField> wetSteamSystem::alpha(const label patchi) const
{
    tmp<scalarField> alphag = gas().alpha(patchi).clone();
    if (alphag()[0] > 0)
    {
        const scalarField& pp = p().boundaryField()[patchi];
        const scalarField& Tp = T().boundaryField()[patchi];
        const scalarField& wp = w().boundaryField()[patchi];
        scalarField& alphap = alphag.ref();
        forAll(alphap, i)
        {
            alphap[i] = (1 - wp[i])*alphap[i] + wp[i]*liquid().properties().alphah(pp[i],Tp[i]);
        }
        return alphag;
    }
    else
    {
        return alphag;
    }
}


tmp<volScalarField> wetSteamSystem::alphahe() const
{
    return alpha();
}


tmp<scalarField> wetSteamSystem::alphahe(const label patchi) const
{
    return alpha(patchi);
}


tmp<volScalarField> wetSteamSystem::alphaEff(const volScalarField& alphat) const
{
    return alpha() + alphat;
}

tmp<scalarField> wetSteamSystem::alphaEff(const scalarField& alphat, const label patchi) const
{
    return alpha(patchi) + alphat;
}

// TODO
tmp<volScalarField> wetSteamSystem::kappa() const
{
    notImplemented(__FUNCTION__);
    return gas().kappa();
}

tmp<scalarField> wetSteamSystem::kappa(const label patchi) const
{
    notImplemented(__FUNCTION__);
    return gas().kappa(patchi);
}

tmp<volScalarField> wetSteamSystem::kappaEff(const volScalarField& alphat) const
{
    notImplemented(__FUNCTION__);
    return kappa();
}

tmp<scalarField> wetSteamSystem::kappaEff(const scalarField& alphat, const label patchi) const
{
    notImplemented(__FUNCTION__);
    return kappa(patchi);
}

// L = hg(p,T) - hl(p,Ts) = hg(p,T) - hg(p,Ts) + hg(p,Ts) - hl(p,Ts) = 
//   ~ hg(p,T) - hg(p,Ts) + 1/rhog(p,Ts)*dps/dt ~ Cp*(T - Ts) + 1/rhog*dpsdt
//
#define L_FROM_SATCURVE
tmp<volScalarField> wetSteamSystem::L() const
{
#ifdef L_FROM_SATCURVE
    volScalarField Ts = TSat();

    const std::function<scalar(scalar,scalar)> f = [&](scalar pp, scalar TT)
        {return gasProps_->rho(pp,TT);};
    volScalarField rhos = applyFunction2(f, gas().p(), Ts, "rhos", dimDensity);

    return tmp<volScalarField>
        (
	    //new volScalarField("L", saturation_->dpsdT(Ts)*Ts/rhos)
            new volScalarField("L", gas().Cp()*(T() - Ts) + saturation_->dpsdT(Ts)*Ts/rhos)
        );
#else
    const std::function<scalar(scalar)> f = [&](scalar TT)
        {return 461.52*(-2.7246e-2*sqr(TT) + 2*1.6853e-5*pow(TT,3) + 2.4576*TT + 6094.4642);};

    return applyFunction1(f, T(), "L", dimEnergy/dimMass);
#endif
}

tmp<volScalarField> wetSteamSystem::S() const
{
    volScalarField Ts = TSat();
    volScalarField Ls = L();

    const std::function<scalar(scalar,scalar)> s = [&](scalar pp, scalar TT)
        {return gasProps_->S(pp,TT);};
    volScalarField Sg =
        applyFunction2(s, gas().p(), gas().T(), "Sg", dimEnergy/dimMass/dimTemperature);

    return tmp<volScalarField>
        (
            new volScalarField("S", Sg - w_*Ls/Ts)
        );
}

void wetSteamSystem::correct()
{
    volScalarField Q = w_*L();
    
    setField(gas().he(), he_ + Q);
    gas().he().correctBoundaryConditions();
    gas().correct();
    setField(he_, gas().he() - Q);
    setField(psi_, gas().psi()/(1-w_));
}


void wetSteamSystem::setField(volScalarField& dest, const volScalarField& src)
{
    dest = src;
    forAll(dest.boundaryField(), patchi)
    {
        scalarField& dp = dest.boundaryFieldRef()[patchi];
        const scalarField& sp = src.boundaryField()[patchi];
        forAll(dp, f)
        {
            dp[f] = sp[f];
        }
    }
}

wordList wetSteamSystem::replaceFixedEnergyTypes(const wordList& bcTypes) const
{
    wordList result;
    for (const auto& bcType : bcTypes)
    {
        result.append( (bcType=="fixedEnergy" ? "fixedValue" : bcType) );
    }
    return result;
}

// ************************************************************************* //

wetSteamSystem::liquidThermo::liquidThermo(const wetSteamSystem& steam, const dictionary& dict)
    :
    steam_(steam),
    liquidProps_(
        liquidProperties::New(dict.lookupOrDefault<word>("liquidProperties", "H2O"))
    )
{
}

const liquidProperties& wetSteamSystem::liquidThermo::properties() const
{
    return liquidProps_();
}

tmp<volScalarField> wetSteamSystem::liquidThermo::rho() const
{
    return rho(steam_.p(), steam_.T());
}

tmp<volScalarField> wetSteamSystem::liquidThermo::rho(
    const volScalarField& p, const volScalarField& T
) const
{
    const std::function<scalar(scalar,scalar)> f = [&](scalar pp, scalar TT)
        {return liquidProps_->rho(pp,TT);};
    return applyFunction2(f, p, T, "rho_l", dimDensity);
}


tmp<volScalarField> wetSteamSystem::liquidThermo::sigma() const
{
    return sigma(steam_.p(), steam_.T());
}

tmp<volScalarField> wetSteamSystem::liquidThermo::sigma(
    const volScalarField& p, const volScalarField& T
) const
{
    const std::function<scalar(scalar,scalar)> f = [&](scalar pp, scalar TT)
        {return liquidProps_->sigma(pp,TT);};
    return applyFunction2(f, p, T, "sigma", dimForce/dimLength);
}


tmp<volScalarField> wetSteamSystem::liquidThermo::mu() const
{
    const std::function<scalar(scalar,scalar)> f = [&](scalar pp, scalar TT)
        {return liquidProps_->mu(pp,TT);};
    return applyFunction2(f, steam_.p(), steam_.T(), "mu", dimPressure*dimTime);
}

tmp<volScalarField> wetSteamSystem::liquidThermo::alpha() const
{
    const std::function<scalar(scalar,scalar)> f = [&](scalar pp, scalar TT)
        {return liquidProps_->alphah(pp,TT);};
    return applyFunction2(f, steam_.p(), steam_.T(), "alpha", dimPressure*dimTime);
}


// ************************************************************************* //

tmp<volScalarField> applyFunction1
(
    const std::function<scalar(scalar)> f,
    const volScalarField& src,
    const word name,
    const dimensionSet dim
) 
{
    const fvMesh& mesh = src.mesh();
    tmp<volScalarField> tres
        (
            new volScalarField
            (
                IOobject
                (
                    name,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh,
                dim
            )
        );
    volScalarField& r = tres.ref();

    forAll(r, celli)
    {
        r[celli] = f(src[celli]);
    }
    
    forAll(src.boundaryField(), patchi)
    {
        scalarField& rp = r.boundaryFieldRef()[patchi];
        const scalarField& sp = src.boundaryField()[patchi];
        forAll(rp, facei)
        {
            rp[facei] = f(sp[facei]);
        }
    }
    
    
    return tres;
}


tmp<volScalarField> applyFunction2
(
    const std::function<scalar(scalar,scalar)> f,
    const volScalarField& src1,
    const volScalarField& src2,
    const word name,
    const dimensionSet dim
) 
{
    const fvMesh& mesh = src1.mesh();
    tmp<volScalarField> tres
        (
            new volScalarField
            (
                IOobject
                (
                    name,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh,
                dim
            )
        );
    volScalarField& r = tres.ref();

    forAll(r, celli)
    {
        r[celli] = f(src1[celli], src2[celli]);
    }
    
    forAll(src1.boundaryField(), patchi)
    {
        scalarField& rp = r.boundaryFieldRef()[patchi];
        const scalarField& s1p = src1.boundaryField()[patchi];
        const scalarField& s2p = src2.boundaryField()[patchi];
        forAll(rp, facei)
        {
            rp[facei] = f(s1p[facei], s2p[facei]);
        }
    }
    
    
    return tres;
}

// ************************************************************************* //

}
}
