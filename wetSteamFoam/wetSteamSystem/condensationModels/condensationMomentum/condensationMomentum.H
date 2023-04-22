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

#ifndef condensationMomentum_H
#define condensationMomentum_H

#include "condensationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace WetSteam
{

/*---------------------------------------------------------------------------* \
                         Class condensationModel Declaration
\*---------------------------------------------------------------------------*/

class condensationMomentum
:
    public condensationModel
{
private:
    // Private data

    //- Disallow default bitwise copy construct
    condensationMomentum(const condensationMomentum&) = delete;
    
    //- Disallow default bitwise assignment
    void operator=(const condensationMomentum&) = delete;

protected:
    //- Number of droplets [-]
    volScalarField Q0_; 

    //- Sum of radii [m]
    volScalarField Q1_;

    //- Sum of radii squared [m2]
    volScalarField Q2_;

    //- Mass flux through gas-liquid interface [kg/(s m3)]
    volScalarField mDotGL_;

    //- Mass of single molecule
    const dimensionedScalar m1_;

    //- Surface tension correction coefficient
    const dimensionedScalar beta_;

    //- Turbulent schmidt number
    const dimensionedScalar Sct_;

    //- Use Kantrowitz correction factor
    Switch kantrowitz_;

    //- Use Courtney correction factor
    Switch courtney_;
    
public:

    TypeName("momentum");
        
    // Constructors

    //- Construct from components
    condensationMomentum(
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& phi,
            wetSteamSystem& steam
    );
    
    //- Destructor
    virtual ~condensationMomentum() {};
    

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
