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
    Foam::saturationIF97

Description

SourceFiles
    saturationIF97.C

\*---------------------------------------------------------------------------*/

#ifndef saturationIF97_H
#define saturationIF97_H

#include "saturationCurve.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace WetSteam
{

/*---------------------------------------------------------------------------*\
                         Class saturationCurve Declaration
\*---------------------------------------------------------------------------*/

class saturationIF97
:
    public saturationCurve
{
private:

    // Critical values
    const scalar Tc = 647.096;
    const scalar pc = 22.064e6;

    const scalar n1 =  0.11670521452767E+04;
    const scalar n2 = -0.72421316703206E+06;
    const scalar n3 = -0.17073846940092E+02;
    const scalar n4 =  0.12020824702470E+05;
    const scalar n5 = -0.32325550322333E+07;
    const scalar n6 =  0.14915108613530E+02;
    const scalar n7 = -0.48232657361591E+04;
    const scalar n8 =  0.40511340542057E+06;
    const scalar n9 = -0.23855557567849E+00;
    const scalar n10=  0.65017534844798E+03;

public:

    TypeName("IF97");

    // Constructor
    saturationIF97(const dictionary& dict);
    
    // Member Functions

    virtual scalar ps(scalar T) const;
    virtual scalar Ts(scalar p) const;
    virtual scalar dpsdT(scalar T) const;

    virtual ~saturationIF97() {};
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
