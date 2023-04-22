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
    Foam::IAPWSIF97Transport

Description
    Transport package using IAPWS IF97.

    Templated into a given thermodynamics package (needed for thermal
    conductivity).

SourceFiles
    IAPWSIF97TransportI.H
    IAPWSIF97Transport.C

\*---------------------------------------------------------------------------*/

#ifndef IAPWSIF97Transport_H
#define IAPWSIF97Transport_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo> class IAPWSIF97Transport;

template<class Thermo>
inline IAPWSIF97Transport<Thermo> operator+
(
    const IAPWSIF97Transport<Thermo>&,
    const IAPWSIF97Transport<Thermo>&
);

template<class Thermo>
inline IAPWSIF97Transport<Thermo> operator*
(
    const scalar,
    const IAPWSIF97Transport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const IAPWSIF97Transport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                    Class IAPWSIF97Transport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class IAPWSIF97Transport
:
    public Thermo
{
    // Private data

public:

    // Constructors

        //- Construct from components
        inline IAPWSIF97Transport
        (
            const Thermo& t
        );

        //- Construct as named copy
        inline IAPWSIF97Transport(const word&, const IAPWSIF97Transport&);

        //- Construct from dictionary
        IAPWSIF97Transport(const dictionary& dict);

        //- Construct from base thermo and dictionary
        IAPWSIF97Transport(const Thermo& t, const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<IAPWSIF97Transport> clone() const;

        // Selector from dictionary
        inline static autoPtr<IAPWSIF97Transport> New(const dictionary& dict);


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "IAPWSIF97<" + Thermo::typeName() + '>';
        }

        //- Dynamic viscosity [kg/ms]
        inline scalar mu(const scalar p, const scalar T) const;

        //- Thermal conductivity [W/mK]
        inline scalar kappa(const scalar p, const scalar T) const;

        //- Thermal diffusivity of enthalpy [kg/ms]
        inline scalar alphah(const scalar p, const scalar T) const;

        // Species diffusivity
        // inline scalar D(const scalar p, const scalar T) const;

        //- Write to Ostream
        void write(Ostream& os) const;


    // Member operators

        inline void operator=(const IAPWSIF97Transport&);

        inline void operator+=(const IAPWSIF97Transport&);

        inline void operator*=(const scalar);


    // Friend operators

        friend IAPWSIF97Transport operator+ <Thermo>
        (
            const IAPWSIF97Transport&,
            const IAPWSIF97Transport&
        );

        friend IAPWSIF97Transport operator* <Thermo>
        (
            const scalar,
            const IAPWSIF97Transport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const IAPWSIF97Transport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "IAPWSIF97TransportI.H"

#ifdef NoRepository
    #include "IAPWSIF97Transport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
