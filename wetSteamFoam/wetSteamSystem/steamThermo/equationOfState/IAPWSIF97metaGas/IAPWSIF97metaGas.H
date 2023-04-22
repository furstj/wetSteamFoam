/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2018 OpenFOAM Foundation
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
    Foam::IAPWSIF97metaGas

Description
    IAPWSIF97metaGas gas equation of state.

SourceFiles
    IAPWSIF97metaGasI.H
    IAPWSIF97metaGas.C

\*---------------------------------------------------------------------------*/

#ifndef IAPWSIF97metaGas_H
#define IAPWSIF97metaGas_H

#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class IAPWSIF97metaGas;

template<class Specie>
inline IAPWSIF97metaGas<Specie> operator+
(
    const IAPWSIF97metaGas<Specie>&,
    const IAPWSIF97metaGas<Specie>&
);

template<class Specie>
inline IAPWSIF97metaGas<Specie> operator*
(
    const scalar,
    const IAPWSIF97metaGas<Specie>&
);

template<class Specie>
inline IAPWSIF97metaGas<Specie> operator==
(
    const IAPWSIF97metaGas<Specie>&,
    const IAPWSIF97metaGas<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const IAPWSIF97metaGas<Specie>&
);



/*---------------------------------------------------------------------------*\
                           Class IAPWSIF97metaGas Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class IAPWSIF97metaGas
:
    public Specie
{
    // Private data
    
    
    
public:
    
    // Constructors
    
    //- Construct from components
    inline IAPWSIF97metaGas(const Specie& sp);
    
    //- Construct from dictionary
    IAPWSIF97metaGas(const dictionary& dict);
    
    //- Construct as named copy
    inline IAPWSIF97metaGas(const word& name, const IAPWSIF97metaGas&);
    
    //- Construct and return a clone
    inline autoPtr<IAPWSIF97metaGas> clone() const;

    // Selector from dictionary
    inline static autoPtr<IAPWSIF97metaGas> New
    (
        const dictionary& dict
    );
    
    
    // Member functions
    
    //- Return the instantiated type name
    static word typeName()
        {
            return "IAPWSIF97metaGas<" + word(Specie::typeName_()) + '>';
        }
    
    // Fundamental properties
    
    
    //- Is the equation of state is incompressible i.e. rho != f(p)
    static const bool incompressible = false;
    
    //- Is the equation of state is isochoric i.e. rho = const
    static const bool isochoric = false;
    
    //- Return density [kg/m^3]
    inline scalar rho(scalar p, scalar T) const;

    //- Return enthalpy departure [J/kg]
    inline scalar H(const scalar p, const scalar T) const;

    //- Return Cp departure [J/(kg K]
    inline scalar Cp(scalar p, scalar T) const;
    
    //- Return Cv departure [J/(kg K]
    inline scalar Cv(scalar p, scalar T) const;
    
    //- Return entropy [J/(kg K)]
    inline scalar S(const scalar p, const scalar T) const;
    
    //- Return compressibility rho/p [s^2/m^2]
    inline scalar psi(scalar p, scalar T) const;
    
    //- Return compression factor []
    inline scalar Z(scalar p, scalar T) const;
    
    //- Return (Cp - Cv) [J/(kg K]
    inline scalar CpMCv(scalar p, scalar T) const;
    
    
    // IO
    
    //- Write to Ostream
    void write(Ostream& os) const;
    
    // Member operators
    
    inline void operator+=(const IAPWSIF97metaGas&);
    inline void operator*=(const scalar);
    
    
    // Friend operators
    
    friend IAPWSIF97metaGas operator+ <Specie>
    (
        const IAPWSIF97metaGas&,
        const IAPWSIF97metaGas&
    );
    
    friend IAPWSIF97metaGas operator* <Specie>
    (
        const scalar s,
        const IAPWSIF97metaGas&
    );
    
    friend IAPWSIF97metaGas operator== <Specie>
    (
        const IAPWSIF97metaGas&,
        const IAPWSIF97metaGas&
    );
    
    
    // Ostream Operator
    
    friend Ostream& operator<< <Specie>
    (
        Ostream&,
        const IAPWSIF97metaGas&
    );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "IAPWSIF97metaGasI.H"

#ifdef NoRepository
    #include "IAPWSIF97metaGas.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
