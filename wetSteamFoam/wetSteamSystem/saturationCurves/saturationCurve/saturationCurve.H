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
    Foam::saturationCurve

Description

SourceFiles
    saturationCurve.C

\*---------------------------------------------------------------------------*/

#ifndef saturationCurve_H
#define saturationCurve_H

#include "IOdictionary.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include <functional>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace WetSteam
{

/*---------------------------------------------------------------------------*\
                         Class saturationCurve Declaration
\*---------------------------------------------------------------------------*/

class saturationCurve
{
    
public:
    TypeName("saturationCurve");
    
    declareRunTimeSelectionTable
    (
        autoPtr,
        saturationCurve,
        dict,
        (
            const dictionary& dict
        ),
        (dict)
    );

    saturationCurve(const dictionary& dict) {};
    
    // Selectors

    //- Select null constructed
    static autoPtr<saturationCurve> New(const dictionary& dict);


    //- Destructor
    virtual ~saturationCurve() {};


    // Member Functions

    //- Saturation pressure at given temperature [Pa]
    virtual scalar ps(scalar T) const = 0;
    virtual tmp<volScalarField> ps(const volScalarField& T) const;
    
    //- Saturation temperature at given pressure [K]
    virtual scalar Ts(scalar p) const = 0;
    virtual tmp<volScalarField> Ts(const volScalarField& p) const;

    //- Derivative of saturation pressure by temperature [Pa/K]
    virtual scalar dpsdT(scalar T) const = 0;
    virtual tmp<volScalarField> dpsdT(const volScalarField& T) const;


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
