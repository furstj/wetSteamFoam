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

#include "specie.H"
#include "../../IF97/IF97.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::IAPWSIF97Transport<Thermo>::IAPWSIF97Transport
(
    const Thermo& t
)
:
    Thermo(t)
{}


template<class Thermo>
inline Foam::IAPWSIF97Transport<Thermo>::IAPWSIF97Transport
(
    const word& name,
    const IAPWSIF97Transport& st
)
:
    Thermo(name, st)
{}


template<class Thermo>
inline Foam::autoPtr<Foam::IAPWSIF97Transport<Thermo>>
Foam::IAPWSIF97Transport<Thermo>::clone() const
{
    return autoPtr<IAPWSIF97Transport<Thermo>>
    (
        new IAPWSIF97Transport<Thermo>(*this)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::IAPWSIF97Transport<Thermo>>
Foam::IAPWSIF97Transport<Thermo>::New
(
    const dictionary& dict
)
{
    return autoPtr<IAPWSIF97Transport<Thermo>>
    (
        new IAPWSIF97Transport<Thermo>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::IAPWSIF97Transport<Thermo>::mu
(
    const scalar p,
    const scalar T
) const
{
    scalar rho = 1./IF97::reg2meta::v(p,T);
    return IF97::mu(rho,T);
}


template<class Thermo>
inline Foam::scalar Foam::IAPWSIF97Transport<Thermo>::kappa
(
    const scalar p, const scalar T
) const
{
    scalar rho = 1./IF97::reg2meta::v(p,T);
    return IF97::lambda(rho,T);
}


template<class Thermo>
inline Foam::scalar Foam::IAPWSIF97Transport<Thermo>::alphah
(
    const scalar p,
    const scalar T
) const
{
    return kappa(p, T)/this->Cp(p, T);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline void Foam::IAPWSIF97Transport<Thermo>::operator=
(
    const IAPWSIF97Transport<Thermo>& st
)
{
    Thermo::operator=(st);
}


template<class Thermo>
inline void Foam::IAPWSIF97Transport<Thermo>::operator+=
(
    const IAPWSIF97Transport<Thermo>& st
)
{
    Thermo::operator+=(st);
}


template<class Thermo>
inline void Foam::IAPWSIF97Transport<Thermo>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::IAPWSIF97Transport<Thermo> Foam::operator+
(
    const IAPWSIF97Transport<Thermo>& st1,
    const IAPWSIF97Transport<Thermo>& st2
)
{
    return IAPWSIF97Transport<Thermo>
    (
        static_cast<const Thermo&>(st1) + static_cast<const Thermo&>(st2)
    );
}


template<class Thermo>
inline Foam::IAPWSIF97Transport<Thermo> Foam::operator*
(
    const scalar s,
    const IAPWSIF97Transport<Thermo>& st
)
{
    return IAPWSIF97Transport<Thermo>
    (
        s*static_cast<const Thermo&>(st)
    );
}


// ************************************************************************* //
