/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::IAPWSIF97metaThermo

Description
    IAPWSIF97META tables based thermodynamics package templated
    into the equation of state.

SourceFiles
    IAPWSIF97metaThermoI.H
    IAPWSIF97metaThermo.C

\*---------------------------------------------------------------------------*/

#ifndef IAPWSIF97metaThermo_H
#define IAPWSIF97metaThermo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class EquationOfState> class IAPWSIF97metaThermo;

template<class EquationOfState>
inline IAPWSIF97metaThermo<EquationOfState> operator+
(
    const IAPWSIF97metaThermo<EquationOfState>&,
    const IAPWSIF97metaThermo<EquationOfState>&
);

template<class EquationOfState>
inline IAPWSIF97metaThermo<EquationOfState> operator*
(
    const scalar,
    const IAPWSIF97metaThermo<EquationOfState>&
);

template<class EquationOfState>
inline IAPWSIF97metaThermo<EquationOfState> operator==
(
    const IAPWSIF97metaThermo<EquationOfState>&,
    const IAPWSIF97metaThermo<EquationOfState>&
);

template<class EquationOfState>
Ostream& operator<<
(
    Ostream&,
    const IAPWSIF97metaThermo<EquationOfState>&
);


/*---------------------------------------------------------------------------*\
                           Class IAPWSIF97metaThermo Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState>
class IAPWSIF97metaThermo
:
    public EquationOfState
{

public:

    // Public data
    
private:
    
    // Private data
    
    
public:
    
    // Constructors
    
    //- Construct from components
    inline IAPWSIF97metaThermo
    (
        const EquationOfState& st
    );
    
    //- Construct from dictionary
    IAPWSIF97metaThermo(const dictionary& dict);
    
    //- Construct as a named copy
    inline IAPWSIF97metaThermo(const word&, const IAPWSIF97metaThermo&);
    
    
    // Member Functions
    
    //- Return the instantiated type name
    static word typeName()
        {
            return "IAPWSIF97meta<" + EquationOfState::typeName() + '>';
        }
    
    //- Limit the temperature to be in the range Tlow_ to Thigh_
    inline scalar limit(const scalar T) const;
    

    // Fundamental properties
    
    //- Heat capacity at constant pressure [J/(kg K)]
    inline scalar Cp(const scalar p, const scalar T) const;

    //- Heat capacity at constant volume [J/(kg K)]
    inline scalar Cv(const scalar p, const scalar T) const;

    //- Absolute Enthalpy [J/kg]
    inline scalar Ha(const scalar p, const scalar T) const;

    //- Sensible enthalpy [J/kg]
    inline scalar Hs(const scalar p, const scalar T) const;

    //- Chemical enthalpy [J/kg]
    inline scalar Hc() const;

    //- Entropy [J/(kg K)]
    inline scalar S(const scalar p, const scalar T) const;


    // Derivative term used for Jacobian

    //- Derivative of Gibbs free energy w.r.t. temperature
    inline scalar dGdT(const scalar p, const scalar T) const;
    
    //- Temperature derivative of heat capacity at constant pressure
    inline scalar dCpdT(const scalar p, const scalar T) const;


    // I-O
    
    //- Write to Ostream
    void write(Ostream& os) const;
    
    
    // Member operators
    
    inline void operator+=(const IAPWSIF97metaThermo&);
    
    
    // Friend operators
    
    friend IAPWSIF97metaThermo operator+ <EquationOfState>
    (
        const IAPWSIF97metaThermo&,
        const IAPWSIF97metaThermo&
    );
    
    friend IAPWSIF97metaThermo operator* <EquationOfState>
    (
        const scalar,
        const IAPWSIF97metaThermo&
    );
    
    friend IAPWSIF97metaThermo operator== <EquationOfState>
    (
        const IAPWSIF97metaThermo&,
        const IAPWSIF97metaThermo&
    );
    
    
    // Ostream Operator
    
    friend Ostream& operator<< <EquationOfState>
    (
        Ostream&,
        const IAPWSIF97metaThermo&
    );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "IAPWSIF97metaThermoI.H"

#ifdef NoRepository
    #include "IAPWSIF97metaThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
