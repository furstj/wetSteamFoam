#ifndef IF97_H
#define IF97_H

#include <cmath>

namespace IF97
{

constexpr double R       = 461.526;   // [J/kg/K]
constexpr double Tc      = 647.096;   // [K]
constexpr double pc      = 22.064e6;  // [Pa]
constexpr double rhoc    = 322.0;     // [kg/m3]

#include "reg2meta.H"

#include "viscosity.H"
#include "conductivity.H"

};


#endif
