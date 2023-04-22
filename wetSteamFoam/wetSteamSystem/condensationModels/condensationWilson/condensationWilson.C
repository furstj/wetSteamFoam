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

#include "condensationWilson.H"
#include "addToRunTimeSelectionTable.H"
#include "fvm.H"

namespace Foam
{
namespace WetSteam
{

defineTypeNameAndDebug(condensationWilson, 0);
addToRunTimeSelectionTable(condensationModel, condensationWilson, params);

condensationWilson::condensationWilson(
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi,
    wetSteamSystem& steam
):
    condensationModel(rho, U, phi, steam),

    WilsonLine_(steam.gas().lookupOrDefault("WilsonLine", 1.0))
{
    Info << "Wilson line at " << WilsonLine_*100 << "% steam quality." << endl << endl;
}


void condensationWilson::correct()
{
    volScalarField& w = steam().w();

    solve(
        fvm::ddt(rho_, w) + fvm::div(phi_, w)
    );

    volScalarField Ts = steam().TSat();
    const volScalarField& T = steam().T();
    volScalarField Cp = steam().gas().Cp();
    volScalarField L = steam().L();

    scalar wCrit = 1 - WilsonLine_;
    
    forAll(T, i)
    {
        scalar weq = w[i] + (1 - w[i])*Cp[i]*(Ts[i] - T[i])/L[i];
        if (weq>wCrit)
        {
            w[i] = weq;
        }
            
    }
    w.correctBoundaryConditions();
}


}
}

// ************************************************************************* //
