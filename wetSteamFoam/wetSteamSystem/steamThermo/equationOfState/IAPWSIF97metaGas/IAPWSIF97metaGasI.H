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

\*---------------------------------------------------------------------------*/

#include "IAPWSIF97metaGas.H"
#include "mathematicalConstants.H"
#include "../../IF97/IF97.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Specie>
inline Foam::IAPWSIF97metaGas<Specie>::IAPWSIF97metaGas(const Specie& sp)
    :
    Specie(sp)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::IAPWSIF97metaGas<Specie>::IAPWSIF97metaGas
(
    const word& name,
    const IAPWSIF97metaGas& pg
)
:
    Specie(name, pg)
{}


template<class Specie>
inline Foam::autoPtr<Foam::IAPWSIF97metaGas <Specie>>
Foam::IAPWSIF97metaGas<Specie>::clone() const
{
    return autoPtr<IAPWSIF97metaGas<Specie>>
    (
        new IAPWSIF97metaGas<Specie>(*this)
    );
}


template<class Specie>
inline Foam::autoPtr<Foam::IAPWSIF97metaGas<Specie>>
Foam::IAPWSIF97metaGas<Specie>::New
(
    const dictionary& dict
)
{
    return autoPtr<IAPWSIF97metaGas<Specie>>
    (
        new IAPWSIF97metaGas<Specie>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::scalar Foam::IAPWSIF97metaGas<Specie>::rho
(
    scalar p,
    scalar T
) const
{
    return 1/IF97::reg2meta::v(p,T);
}


template<class Specie>
inline Foam::scalar Foam::IAPWSIF97metaGas<Specie>::H(scalar p, scalar T) const
{
    return IF97::reg2meta::hr(p,T);
}


template<class Specie>
inline Foam::scalar Foam::IAPWSIF97metaGas<Specie>::Cp(scalar p, scalar T) const
{
    return IF97::reg2meta::cpr(p,T);
}

template<class Specie>
inline Foam::scalar Foam::IAPWSIF97metaGas<Specie>::Cv(scalar p, scalar T) const
{
    return IF97::reg2meta::cvr(p,T);
}


template<class Specie>
inline Foam::scalar Foam::IAPWSIF97metaGas<Specie>::S
(
    scalar p,
    scalar T
) const
{
    return IF97::reg2meta::sr(p,T);
}


template<class Specie>
inline Foam::scalar Foam::IAPWSIF97metaGas<Specie>::psi
(
    scalar p,
    scalar T
) const
{
    return 1.0/(p*IF97::reg2meta::v(p,T));
}


template<class Specie>
inline Foam::scalar Foam::IAPWSIF97metaGas<Specie>::Z
(
    scalar p,
    scalar T
) const
{
    return IF97::reg2meta::z(p,T);
}


template<class Specie>
inline Foam::scalar Foam::IAPWSIF97metaGas<Specie>::CpMCv
(
    scalar p,
    scalar T
) const
{
    return IF97::reg2meta::cp(p,T) - IF97::reg2meta::cv(p,T);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::IAPWSIF97metaGas<Specie>::operator+=
(
    const IAPWSIF97metaGas<Specie>& pg
)
{
    Specie::operator+=(pg);
}


template<class Specie>
inline void Foam::IAPWSIF97metaGas<Specie>::operator*=(const scalar s)
{
     Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


template<class Specie>
Foam::IAPWSIF97metaGas<Specie> Foam::operator+
(
    const IAPWSIF97metaGas<Specie>& pg1,
    const IAPWSIF97metaGas<Specie>& pg2
)
{
    return IAPWSIF97metaGas<Specie>
        (
            static_cast<const Specie&>(pg1)
            + static_cast<const Specie&>(pg2)
        );
}


template<class Specie>
Foam::IAPWSIF97metaGas<Specie> Foam::operator*
(
    const scalar s,
    const IAPWSIF97metaGas<Specie>& pg
)
{
    return IAPWSIF97metaGas<Specie>(s*static_cast<const Specie&>(pg));
}


template<class Specie>
Foam::IAPWSIF97metaGas<Specie> Foam::operator==
(
    const IAPWSIF97metaGas<Specie>& pg1,
    const IAPWSIF97metaGas<Specie>& pg2
)
{
    return IAPWSIF97metaGas<Specie>
        (
            static_cast<const Specie&>(pg1) == static_cast<const Specie&>(pg2)
        );
}


// ************************************************************************* //
