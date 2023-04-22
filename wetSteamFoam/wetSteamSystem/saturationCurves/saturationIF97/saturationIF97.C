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

#include "saturationIF97.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace WetSteam
{

defineTypeNameAndDebug(saturationIF97, 0);
addToRunTimeSelectionTable(saturationCurve, saturationIF97, dict);


saturationIF97::saturationIF97(const dictionary& dict):
    saturationCurve(dict)
{}


scalar saturationIF97::dpsdT(scalar T) const
{
    scalar psat = ps(T);
    scalar beta = pow(psat*1e-6, 0.25);
    scalar theta = T + n9/(T - n10);

    // Derivatives of implicit equation (29) with respect to beta and theta
    scalar dfdbeta = (2*beta + n3)*sqr(theta) + (2*beta*n1 + n4)*theta + 2*n2*beta + n5;
    scalar dfdtheta= (2*theta + n1)*sqr(beta) + (2*n3*theta + n4)*beta + 2*n6*theta + n7; 

    scalar dthetadT = (1 - n9/sqr(T - n10));
    scalar dbetadtheta = -dfdtheta/dfdbeta;
    scalar dpdbeta = 4*pow3(beta)*1e6;

    return dpdbeta*dbetadtheta*dthetadT;
}


scalar saturationIF97::ps(scalar T) const
{
    T = max(min(T, 647.096), 273.15);
    scalar theta = T + n9/(T - n10);
    scalar A = sqr(theta) + n1*theta + n2;
    scalar B = n3*sqr(theta) + n4*theta + n5;
    scalar C = n6*sqr(theta) + n7*theta + n8;
    
    scalar p = pow4(2*C/(-B + sqrt(sqr(B) - 4*A*C)))*1e6;
    return p;
}


scalar saturationIF97::Ts(scalar p) const
{
    p = max(min(p, 22.064e6), 611.213);
    
    scalar beta = pow(p*1e-6, 0.25);
    
    scalar E = sqr(beta) + n3*beta + n6;
    scalar F = n1*sqr(beta) + n4*beta + n7;
    scalar G = n2*sqr(beta) + n5*beta + n8;
    scalar D = 2*G/(-F - sqrt(sqr(F) - 4*E*G));

    scalar T = 0.5*(n10 + D - sqrt(sqr(n10 + D) - 4*(n9 + n10*D)));
    return T;
}


}
}

