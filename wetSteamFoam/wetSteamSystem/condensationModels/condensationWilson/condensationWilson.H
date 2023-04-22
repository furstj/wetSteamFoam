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
    Foam::condensationModel

Description

SourceFiles
    condensationModel.C

\*---------------------------------------------------------------------------*/

#ifndef condensationWilson_H
#define condensationWilson_H

#include "condensationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace WetSteam
{

/*---------------------------------------------------------------------------* \
                         Class condensationModel Declaration
\*---------------------------------------------------------------------------*/

class condensationWilson
:
    public condensationModel
{
private:
    // Private data

    //- Disallow default bitwise copy construct
    condensationWilson(const condensationWilson&) = delete;
    
    //- Disallow default bitwise assignment
    void operator=(const condensationWilson&) = delete;

    scalar WilsonLine_;
    
public:

    TypeName("Wilson");
        
    // Constructors

    //- Construct from components
    condensationWilson(
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& phi,
            wetSteamSystem& steam
    );
    
    //- Destructor
    virtual ~condensationWilson() {};
    

    // Member Functions

    //- Main update/correction routine
    virtual void correct();
    
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
