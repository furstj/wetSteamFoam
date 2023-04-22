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

Class
    Foam::wetSteamThermo

Description

SourceFiles
    wetSteamThermo.C

\*---------------------------------------------------------------------------*/

#ifndef wetSteamThermo_H
#define wetSteamThermo_H

#include "fluidThermo.H"
#include "liquidProperties.H"
#include "saturationCurve.H"
#include "gasProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace WetSteam
{

/*---------------------------------------------------------------------------* \
    Class wetSteamSystem declaration
\*---------------------------------------------------------------------------*/

class wetSteamSystem
{
protected:
    
    class liquidThermo
    {
        const wetSteamSystem& steam_;
        autoPtr<liquidProperties> liquidProps_;
        liquidThermo(const liquidThermo&) = delete;
        void operator=(const liquidThermo&) = delete;
        
    public:
        liquidThermo(const wetSteamSystem& steam, const dictionary& dict);

        const liquidProperties& properties() const;
        
        tmp<volScalarField> rho() const;
        tmp<volScalarField> rho(const volScalarField& p, const volScalarField& T) const;

        tmp<volScalarField> mu() const;
        tmp<volScalarField> alpha() const;

        tmp<volScalarField> sigma() const;
        tmp<volScalarField> sigma(const volScalarField& p, const volScalarField& T) const;
};
    
    // Protected data

    const fvMesh& mesh_;

    //- Thermodynamics of gaseous phase
    autoPtr<fluidThermo> gas_;

    //- Thermohysical properties of gaseous phase
    autoPtr<gasProperties> gasProps_;
    
    //- Thermophysical properties of liquid phase
    liquidThermo liquid_;
    
    //- Wetness (mass fraction of wet phase) [-]
    volScalarField w_;

    //- Mixture enthalpy [J/kg]
    volScalarField he_;

    //- Mixture compressibility [s^2/m^2]
    volScalarField psi_;

    //- Saturation curve
    autoPtr<saturationCurve> saturation_;

private:
    // Private Member Functions
    
    //- Construct as copy (not implemented)
    wetSteamSystem(const wetSteamSystem&) = delete;
    void operator=(const wetSteamSystem&) = delete;

    static void setField(volScalarField& dest, const volScalarField& src);

    //- Liquid phase density
    // scalar rho_l(scalar T) const;

    wordList replaceFixedEnergyTypes(const wordList& bcTypes) const;
    
public:
    
    // Constructors
    
    //- Construct from mesh and phase name
    wetSteamSystem
    (
        const fvMesh&,
        const word& phaseName = word::null
    );
    
    
    //- Destructor
    virtual ~wetSteamSystem()
        {};

    //- Access to gas phase
    fluidThermo& gas()
        {
            return gas_();
        }

    const fluidThermo& gas() const
        {
            return gas_();
        }

    const gasProperties& gasProps()
        {
	    return *gasProps_;
	}

    //- Access to liquid phase properties
    const liquidThermo& liquid() const
        {
            return liquid_;
        }
    
    const fvMesh& mesh() const
        {
            return mesh_;
        }
    
    //- Access to wetness
    volScalarField& w()
        {
            return w_;
        }

    const volScalarField& w() const
        {
            return w_;
        }

    volScalarField& psi()
	{
	    return psi_;
	}

    //- Access to mixture pressure [Pa]
    volScalarField& p()
        {
            return gas().p();
        }

    const volScalarField& p() const
        {
            return gas().p();
        }

    //- Access to mixture temperature [K]
    const volScalarField& T() const
        {
            return gas().T();
        }

    
    //- Access to mixture enthalpy [J/kg]
    volScalarField& he()
        {
            return he_;
        }

    //- Mixture density [kg/m3]
    tmp<volScalarField> rho() const;
    tmp<scalarField> rho(label patchi) const;

    //- Mixture dynamic viscosity [Pa s]
    tmp<volScalarField> mu() const;
    tmp<scalarField> mu(label patchi) const;

    //- Mixture kinematic viscosity [m2 s]
    tmp<volScalarField> nu() const;
    tmp<scalarField> nu(label patchi) const;

    //- Thermal diffusivity for enthalpy of mixture [kg/m/s]
    tmp<volScalarField> alpha() const;
    tmp<scalarField> alpha(label patchi) const;

    //- Thermal diffusivity for enthalpy of mixture [kg/m/s]
    tmp<volScalarField> alphahe() const;
    tmp<scalarField> alphahe(label patchi) const;

    tmp<volScalarField> alphaEff(const volScalarField& alphat) const;
    tmp<scalarField> alphaEff(const scalarField& alphat, const label patchi) const;
    
    //- Thermal diffusivity for temperature of mixture [J/m/s/K]
    tmp<volScalarField> kappa() const;
    tmp<scalarField> kappa(label patchi) const;

    tmp<volScalarField> kappaEff(const volScalarField& alphat) const;
    tmp<scalarField> kappaEff(const scalarField& alphat, const label patchi) const;

    //- Saturation pressure at current temperature
    tmp<volScalarField> pSat() const
        {
            return tmp<volScalarField>
                (
                    new volScalarField("pSat", saturation_->ps(gas().T()))
                );
        }
    
    //- Saturation temperature at current pressure
    tmp<volScalarField> TSat() const
        {
            return tmp<volScalarField>
                (
                    new volScalarField("TSat", saturation_->Ts(gas().p()))
                );
        }

    //- Latent heat 
    tmp<volScalarField> L() const;

    //- Mixture entropy
    tmp<volScalarField> S() const;
    
    // Calculate and correct gas phase properties
    void correct();
    
};


tmp<volScalarField> applyFunction1
(
    const std::function<scalar(scalar)> f,
    const volScalarField& src,
    const word name,
    const dimensionSet dim
);


tmp<volScalarField> applyFunction2
(
    const std::function<scalar(scalar,scalar)> f,
    const volScalarField& src1,
    const volScalarField& src2,
    const word name,
    const dimensionSet dim
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
}
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#ifdef NoRepository
//    #include "wetSteamSystem.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
