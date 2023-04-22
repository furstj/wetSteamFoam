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
    Foam::saturationIAPWS95

Description

SourceFiles
    saturationIAPWS95.C

\*---------------------------------------------------------------------------*/

#ifndef saturationIAPWS95_H
#define saturationIAPWS95_H

#include "saturationCurve.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace WetSteam
{

/*---------------------------------------------------------------------------*\
                         Class saturationCurve Declaration
\*---------------------------------------------------------------------------*/

class saturationIAPWS95
:
    public saturationCurve
{
private:

    // Critical values
    const double Tc = 647.096;
    const double pc = 22.064e6;

    // Saturation line pSat(T) (Wagner & Pruss, 2000)
    const double a1 = -7.85951783;
    const double a2 = 1.84408259;
    const double a3 = -11.7866497;
    const double a4 = 22.6807411;
    const double a5 = -15.9618719;
    const double a6 = 1.80122502;

public:

    TypeName("IAPWS95");

    // Constructor
    saturationIAPWS95(const dictionary& dict);
    
    // Member Functions

    virtual scalar ps(scalar T) const;
    virtual scalar Ts(scalar p) const;
    virtual scalar dpsdT(scalar T) const;

    virtual ~saturationIAPWS95() {};
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
