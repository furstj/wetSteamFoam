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

#ifndef condensationModel_H
#define condensationModel_H

#include "IOdictionary.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "wetSteamSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace WetSteam
{

/*---------------------------------------------------------------------------*\
                         Class condensationModel Declaration
\*---------------------------------------------------------------------------*/

class condensationModel
{
protected:
    // Prprotected data

    //- Reference to wetSteamSystem
    wetSteamSystem& steam_;

    //- Refernce to the mesh
    const fvMesh& mesh_;

    //- Time
    const Time& time_;

    const volScalarField& rho_;
    const volVectorField& U_;
    const surfaceScalarField& phi_;
    
    // Private Member Functions

private:
    //- Disallow default bitwise copy construct
    condensationModel(const condensationModel&) = delete;

    //- Disallow default bitwise assignment
    void operator=(const condensationModel&) = delete;


public:
    TypeName("condensationModel");

    // Declare runtime constructor selection table

    declareRunTimeSelectionTable
    (
        autoPtr,
        condensationModel,
        params,
        (
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& phi,
            wetSteamSystem& steam
        ),
        (rho, U, phi, steam)
    );

    
    // Constructors

    //- Construct from components
    condensationModel(
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& phi,
            wetSteamSystem& steam
    );
       

    // Selectors

    
    //- Return a reference to the selected condensation model
    static autoPtr<condensationModel> New(
        const volScalarField& rho,
        const volVectorField& U,
        const surfaceScalarField& phi,
        wetSteamSystem& steam
    );
    

    //- Destructor
    virtual ~condensationModel() = default;
    

    // Member Functions

    //- Main update/correction routine
    virtual void correct() = 0;
    
    // Access
    const wetSteamSystem& steam() const
        {
            return steam_;
        }
    
    wetSteamSystem& steam()
        {
            return steam_;
        }

    const fvMesh& mesh() const
        {
            return mesh_;
        }
    // Check
    
    // Edit
    
    // Write
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

}
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
