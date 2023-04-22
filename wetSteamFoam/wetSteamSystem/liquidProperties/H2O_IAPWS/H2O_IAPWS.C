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

#include "H2O_IAPWS.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(H2O_IAPWS, 0);
    addToRunTimeSelectionTable(liquidProperties, H2O_IAPWS,);
    addToRunTimeSelectionTable(liquidProperties, H2O_IAPWS, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::H2O_IAPWS::H2O_IAPWS()
:
    liquidProperties
    (
        18.015286,
        647.096,
        22.064e+6,
        1.0/322.0,
        0.229,
        273.16,
        611.655,
        373.124,
        6.1709e-30,
        0.3449,
        4.7813e+4
    ),
    rho_(),
    pv_(73.649, -7258.2, -7.3037, 4.1653e-06, 2),
    hl_(647.13, 2889425.47876769, 0.3199, -0.212, 0.25795, 0),
    Cp_
    (
        15341.1046350264,
       -116.019983347211,
        0.451013044684985,
       -0.000783569247849015,
        5.20127671384957e-07,
        0
    ),
    h_
    (
       -17957283.7993676,
        15341.1046350264,
       -58.0099916736053,
        0.150337681561662,
       -0.000195892311962254,
        1.04025534276991e-07
    ),
    Cpg_
    (
        1851.73466555648,
        1487.53816264224,
        2609.3,
        493.366638912018,
        1167.6
    ),
    B_
    (
       -0.0012789342214821,
        1.4909797391063,
       -1563696.91923397,
        1.85445462114904e+19,
       -7.68082153760755e+21
    ),
    mu_(-51.964, 3670.6, 5.7331, -5.3495e-29, 10),
    mug_(2.6986e-06, 0.498, 1257.7, -19570),
    kappa_(-0.4267, 0.0056903, -8.0065e-06, 1.815e-09, 0, 0),
    kappag_(6.977e-05, 1.1243, 844.9, -148850),         
    sigma_(),
    D_(15.0, 15.0, 18.015, 28)
{}


Foam::H2O_IAPWS::H2O_IAPWS
(
    const liquidProperties& l,
    const IAPWS_rhoSat& density,
    const NSRDSfunc1& vapourPressure,
    const NSRDSfunc6& heatOfVapourisation,
    const NSRDSfunc0& heatCapacity,
    const NSRDSfunc0& enthalpy,
    const NSRDSfunc7& idealGasHeatCapacity,
    const NSRDSfunc4& secondVirialCoeff,
    const NSRDSfunc1& dynamicViscosity,
    const NSRDSfunc2& vapourDynamicViscosity,
    const NSRDSfunc0& thermalConductivity,
    const NSRDSfunc2& vapourThermalConductivity,
    const IAPWS_sigma& surfaceTension,
    const APIdiffCoefFunc& vapourDiffussivity
)
:
    liquidProperties(l),
    rho_(density),
    pv_(vapourPressure),
    hl_(heatOfVapourisation),
    Cp_(heatCapacity),
    h_(enthalpy),
    Cpg_(idealGasHeatCapacity),
    B_(secondVirialCoeff),
    mu_(dynamicViscosity),
    mug_(vapourDynamicViscosity),
    kappa_(thermalConductivity),
    kappag_(vapourThermalConductivity),
    sigma_(surfaceTension),
    D_(vapourDiffussivity)
{}


Foam::H2O_IAPWS::H2O_IAPWS(const dictionary& dict)
:
    H2O_IAPWS()
{
    readIfPresent(*this, dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::H2O_IAPWS::writeData(Ostream& os) const
{
    liquidProperties::writeData(*this, os);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const H2O_IAPWS& l)
{
    l.writeData(os);
    return os;
}


// ************************************************************************* //
