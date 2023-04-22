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

\*---------------------------------------------------------------------------*/

#include "IAPWSIF97metaThermo.H"
#include "specie.H"
#include "../../IF97/reg2meta.H"
#include <cassert>

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::IAPWSIF97metaThermo<EquationOfState>::IAPWSIF97metaThermo
(
    const EquationOfState& st
)
:
    EquationOfState(st)
{}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::IAPWSIF97metaThermo<EquationOfState>::IAPWSIF97metaThermo
(
    const word& name,
    const IAPWSIF97metaThermo& jt
)
:
    EquationOfState(name, jt)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::scalar Foam::IAPWSIF97metaThermo<EquationOfState>::limit
(
    const scalar T
) const
{
    return min(max(T, 273.15), 623.15);
}


template<class EquationOfState>
inline Foam::scalar Foam::IAPWSIF97metaThermo<EquationOfState>::Cp
(
    const scalar p,
    const scalar T
) const
{
    return IF97::reg2meta::cp0(p,T) + EquationOfState::Cp(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::IAPWSIF97metaThermo<EquationOfState>::Cv
(
    const scalar p,
        const scalar T
) const
{
    return IF97::reg2meta::cv0(p,T) + EquationOfState::Cv(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::IAPWSIF97metaThermo<EquationOfState>::Ha
(
    const scalar p,
    const scalar T
) const
{
    return IF97::reg2meta::h0(p,T) + EquationOfState::H(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::IAPWSIF97metaThermo<EquationOfState>::Hs
(
    const scalar p,
    const scalar T
) const
{
    return Ha(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::IAPWSIF97metaThermo<EquationOfState>::Hc() const
{
    return 0;
}


template<class EquationOfState>
inline Foam::scalar Foam::IAPWSIF97metaThermo<EquationOfState>::S
(
    const scalar p,
    const scalar T
) const
{
    return IF97::reg2meta::s0(p,T) + EquationOfState::S(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::IAPWSIF97metaThermo<EquationOfState>::dGdT
(
    const scalar p,
    const scalar T
) const
{
    return -S(p,T);
}


template<class EquationOfState>
inline Foam::scalar Foam::IAPWSIF97metaThermo<EquationOfState>::dCpdT
(
    const scalar p,
    const scalar T
) const
{
    return ::IF97::reg2meta::dcpdT(p,T);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline void Foam::IAPWSIF97metaThermo<EquationOfState>::operator+=
(
    const IAPWSIF97metaThermo<EquationOfState>& jt
)
{
    EquationOfState::operator+=(jt);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::IAPWSIF97metaThermo<EquationOfState> Foam::operator+
(
    const IAPWSIF97metaThermo<EquationOfState>& jt1,
    const IAPWSIF97metaThermo<EquationOfState>& jt2
)
{
    EquationOfState eofs = jt1;
    eofs += jt2;
    
    return IAPWSIF97metaThermo<EquationOfState>(eofs);
}


template<class EquationOfState>
inline Foam::IAPWSIF97metaThermo<EquationOfState> Foam::operator*
(
    const scalar s,
    const IAPWSIF97metaThermo<EquationOfState>& jt
)
{
    return IAPWSIF97metaThermo<EquationOfState>
    (
        s*static_cast<const EquationOfState&>(jt)
    );
}


template<class EquationOfState>
inline Foam::IAPWSIF97metaThermo<EquationOfState> Foam::operator==
(
    const IAPWSIF97metaThermo<EquationOfState>& jt1,
    const IAPWSIF97metaThermo<EquationOfState>& jt2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(jt1)
     == static_cast<const EquationOfState&>(jt2)
    );

    return IAPWSIF97metaThermo<EquationOfState>
    (
        eofs
    );
}


// ************************************************************************* //
