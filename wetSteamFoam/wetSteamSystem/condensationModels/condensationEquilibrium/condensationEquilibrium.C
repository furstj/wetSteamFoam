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

#include "condensationEquilibrium.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace WetSteam
{

defineTypeNameAndDebug(condensationEquilibrium, 0);
addToRunTimeSelectionTable(condensationModel, condensationEquilibrium, params);

condensationEquilibrium::condensationEquilibrium(
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi,
    wetSteamSystem& steam
):
    condensationModel(rho, U, phi, steam)
{
}


void condensationEquilibrium::correct()
{
    volScalarField Ts = steam().TSat();
    const volScalarField& T = steam().T();
    volScalarField Cp = steam().gas().Cp();
    volScalarField L = steam().L();
    volScalarField& w = steam().w();
    
    forAll(T, i)
    {
        w[i] += (Ts[i] - T[i])*Cp[i]/L[i];
        w[i] = max(w[i], 0);
    }

}


}
}

// ************************************************************************* //
