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

#include "saturationIAPWS95.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace WetSteam
{

defineTypeNameAndDebug(saturationIAPWS95, 0);
addToRunTimeSelectionTable(saturationCurve, saturationIAPWS95, dict);


saturationIAPWS95::saturationIAPWS95(const dictionary& dict):
    saturationCurve(dict)
{}


scalar saturationIAPWS95::dpsdT(scalar T) const
{
    double theta = 1 - T/Tc;
    double theta05 = sqrt(theta);
    double theta15 = theta*theta05;
    double theta20 = theta*theta;
    double theta25 = theta20*theta05;
    double theta30 = theta20*theta;
    double theta35 = theta30*theta05;
    double theta40 = theta20*theta20;
    double theta65 = theta40*theta25;
    double theta75 = theta40*theta35;

    double lnP = Tc/T*(a1*theta + a2*theta15 + a3*theta30 + a4*theta35 + a5*theta40 + a6*theta75);
    double p = pc*exp(lnP);
    
    return -p/T*(lnP + a1 + 1.5*a2*theta05 + 3*a3*theta20 +
    3.5*a4*theta25 + 4*a5*theta30 + 7.5*a6*theta65);
}


scalar saturationIAPWS95::ps(scalar T) const
{
    double theta = 1 - T/Tc;
    double theta05 = sqrt(theta);
    double theta15 = theta*theta05;
    double theta30 = theta*theta*theta;
    double theta35 = theta30*theta05;
    double theta40 = theta30*theta;
    double theta75 = theta40*theta35;
    
    return pc*exp(Tc/T*(
        a1*theta + a2*theta15 + a3*theta30 + a4*theta35 + a5*theta40 + a6*theta75));
}


scalar saturationIAPWS95::Ts(scalar p) const
{
    double T = 647.096;
    
    // Newton method for f = ps(T) - p = 0    
    int iter = 0;
    while (iter<10) {
        double f = ps(T) - p;
        if (fabs(f)<1e-4*p) break;
        T -= f / dpsdT(T);
        iter++;
    }

    return T;
}


}
}

