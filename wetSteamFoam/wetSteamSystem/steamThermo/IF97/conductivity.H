#ifndef IF97_conductivity_H
#define IF97_conductivity_H

#include <cmath>

// Interpolation of lambda as in https://srd.nist.gov/JPCRD/jpcrd256.pdf
//
// Calculation omits the lambda2 correction for region near critical parameters

// Thermal conductivity in [W/mK]
inline double lambda(double rho, double T)
{
    const double tau = T/647.3;
    const double delta = rho/317.7;
    
    constexpr double a[] = {0.0102811, 0.0299621, 0.0156146, -0.00422464};
    double lambda0 = 0;
    for (int i=0; i<4; i++)
        lambda0 += a[i]*pow(tau,i);
    lambda0 *= sqrt(tau);

    constexpr double b0 = -0.397070;
    constexpr double b1 = 0.400302;
    constexpr double b2 = 1.060000;
    constexpr double B1 = -0.171587;
    constexpr double B2 = 2.392190;

    double lambda1 = b0 + b1*delta + b2*exp(B1*(delta + B2)*(delta + B2));
    return lambda0 + lambda1;
}

#endif
