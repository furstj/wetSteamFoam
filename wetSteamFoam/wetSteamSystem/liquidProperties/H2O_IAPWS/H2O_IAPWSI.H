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

Foam::scalar Foam::H2O_IAPWS::IAPWS_rhoSat::f(scalar, scalar T) const
{
    const scalar Tc   = 647.096;
    const scalar rhoc = 322.0;
    const scalar tau = max(1.0 - T/Tc, 0.0);
    const scalar tau1  = pow(tau, 1.0/3.0);
    const scalar tau2  = sqr(tau1);
    const scalar tau5  = sqr(tau2)*tau1;
    const scalar tau11 = sqr(tau5)*tau1;
    const scalar tau16 = tau11*tau5;
    const scalar tau43 = sqr(tau16)*tau11;
    const scalar tau110= pow(tau11, 10);

    return rhoc*(1.0 + 1.99274064*tau1 + 1.09965342*tau2
    - 0.510839303*tau5 - 1.75493479*tau16 - 45.5170352*tau43
    - 6.74694450e5*tau110);
}

Foam::scalar Foam::H2O_IAPWS::IAPWS_sigma::f(scalar, scalar T) const
{
    const scalar Tc   = 647.096;
    const scalar tau = max(1.0 - T/Tc, 0.0);

    return 235.8e-3*pow(tau, 1.256)*(1 - 0.625*tau);
}


inline Foam::scalar Foam::H2O_IAPWS::rho(scalar p, scalar T) const
{
    return rho_.f(p, T);
}

inline Foam::scalar Foam::H2O_IAPWS::pv(scalar p, scalar T) const
{
    return pv_.f(p, T);
}


inline Foam::scalar Foam::H2O_IAPWS::hl(scalar p, scalar T) const
{
    return hl_.f(p, T);
}


inline Foam::scalar Foam::H2O_IAPWS::Cp(scalar p, scalar T) const
{
    return Cp_.f(p, T);
}


inline Foam::scalar Foam::H2O_IAPWS::h(scalar p, scalar T) const
{
    return h_.f(p, T);
}


inline Foam::scalar Foam::H2O_IAPWS::Cpg(scalar p, scalar T) const
{
    return Cpg_.f(p, T);
}


inline Foam::scalar Foam::H2O_IAPWS::B(scalar p, scalar T) const
{
    return B_.f(p, T);
}


inline Foam::scalar Foam::H2O_IAPWS::mu(scalar p, scalar T) const
{
    return mu_.f(p, T);
}


inline Foam::scalar Foam::H2O_IAPWS::mug(scalar p, scalar T) const
{
    return mug_.f(p, T);
}


inline Foam::scalar Foam::H2O_IAPWS::kappa(scalar p, scalar T) const
{
    return kappa_.f(p, T);
}


inline Foam::scalar Foam::H2O_IAPWS::kappag(scalar p, scalar T) const
{
    return kappag_.f(p, T);
}


inline Foam::scalar Foam::H2O_IAPWS::sigma(scalar p, scalar T) const
{
    return sigma_.f(p, T);
}


inline Foam::scalar Foam::H2O_IAPWS::D(scalar p, scalar T) const
{
    return D_.f(p, T);
}


inline Foam::scalar Foam::H2O_IAPWS::D(scalar p, scalar T, scalar Wb) const
{
    return D_.f(p, T, Wb);
}


// ************************************************************************* //
