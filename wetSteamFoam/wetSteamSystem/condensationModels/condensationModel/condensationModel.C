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

#include "condensationModel.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace WetSteam
{

defineTypeNameAndDebug(condensationModel, 0);
defineRunTimeSelectionTable(condensationModel, params);


autoPtr<condensationModel>
condensationModel::New
(
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi,
    wetSteamSystem& steam
)
{
    word condensationModelName
    (
        steam.gas().lookup("condensationModel")
    );

    Info<< "Selecting condensation model "
        << condensationModelName << endl;

    auto cstrIter = paramsConstructorTablePtr_->cfind(condensationModelName);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown condesnationModel "
            << condensationModelName << endl << endl
            << "Valid condensationModels are : " << endl
            << paramsConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<condensationModel>(cstrIter()(rho, U, phi, steam));
}


condensationModel::condensationModel(
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi,
    wetSteamSystem& steam
):
    steam_(steam),
    mesh_(U.mesh()),
    time_(U.time()),
    rho_(rho),
    U_(U),
    phi_(phi)
{
}

}
}
// ************************************************************************* //
