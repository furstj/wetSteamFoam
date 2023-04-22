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

#include "saturationCurve.H"
#include "wetSteamSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace WetSteam
{

defineTypeNameAndDebug(saturationCurve, 0);
defineRunTimeSelectionTable(saturationCurve, dict);


autoPtr<saturationCurve>
saturationCurve::New
(
    const dictionary& dict
)
{
    word name
    (
        dict.lookup("saturationCurve")
    );

    Info<< "Selecting saturation curve "
        << name << endl;

    auto cstrIter = dictConstructorTablePtr_->cfind(name);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown saturationCurve "
            << name << endl << endl
            << "Valid saturation curves are : " << endl
            << dictConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<saturationCurve>(cstrIter()(dict));
}

tmp<volScalarField> saturationCurve::ps(const volScalarField& T) const
{
    const std::function<scalar(scalar)> f = [&](scalar x){return ps(x);};
    return applyFunction1(f, T, "pSat", dimPressure);
}

//- Saturation temperature at given pressure [K]
tmp<volScalarField> saturationCurve::Ts(const volScalarField& p) const
{
    const std::function<scalar(scalar)> f = [&](scalar x){return Ts(x);};
    return applyFunction1(f, p, "TSat", dimTemperature);
}

//- Derivative of saturation pressure by temperature [Pa/K]
tmp<volScalarField> saturationCurve::dpsdT(const volScalarField& T) const
{
    const std::function<scalar(scalar)> f = [&](scalar x){return dpsdT(x);};
    return applyFunction1(f, T, "dpsdT", dimPressure/dimTemperature);
}

}
}
