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

#include "makeGasProperties.H"

#include "specie.H"
#include "perfectGas.H"
#include "IAPWSIF97metaGas.H"
#include "IAPWSIF97metaThermo.H"
#include "IAPWSIF97Transport.H"
#include "sensibleEnthalpy.H"

#include "thermo.H"

#include "constTransport.H"
#include "sutherlandTransport.H"

#include "hPolynomialThermo.H"
#include "polynomialTransport.H"


namespace Foam
{

makeGasProperties
(
    constTransport,
    sensibleEnthalpy,
    hPolynomialThermo,
    perfectGas,
    specie
);

makeGasProperties
(
    sutherlandTransport,
    sensibleEnthalpy,
    hPolynomialThermo,
    perfectGas,
    specie
);

makeGasProperties
(
    polynomialTransport,
    sensibleEnthalpy,
    hPolynomialThermo,
    perfectGas,
    specie
);


// IAPWS IF97
makeGasProperties
(
    constTransport,
    sensibleEnthalpy,
    IAPWSIF97metaThermo,
    IAPWSIF97metaGas,
    specie
);

makeGasProperties(
    IAPWSIF97Transport,
    sensibleEnthalpy,
    IAPWSIF97metaThermo,
    IAPWSIF97metaGas,
    specie
);

}

// ************************************************************************* //
