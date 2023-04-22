#ifndef IF97_viscosity_H
#define IF97_viscosity_H

#include <cmath>

// Dynamic viscosity
inline double mu(double rho, double T)
{
    constexpr double muStar = 1.e-6;
    
    double delta = rho/rhoc;
    double tau = Tc/T;

    double s = 0;
    constexpr double H0[] = {1.67752, 2.20462, 0.6366564, -0.241605};
    for (int i=0; i<4; i++)
        s += H0[i]*pow(tau,i);
    double mu0 = 100/sqrt(tau)/s;

    
    constexpr double H1[6][7] = {
        { 5.20094E-1, 2.22531E-1, -2.81378E-1, 1.61913E-1, -3.25372E-2, 0.0,         0.0},
        { 8.50895E-2, 9.99115E-1, -9.06851E-1, 2.57399E-1,  0.0,        0.0,         0.0},
        {-1.08374,    1.88797,    -7.72479E-1, 0.0,         0.0,        0.0,         0.0},
        {-2.89555E-1, 1.26613,    -4.89837E-1, 0.0,         6.98452E-2, 0.0,        -4.35673E-3},
        { 0.0,        0.0,        -2.57040E-1, 0.0,         0.0,        8.72102E-3,  0.0},
        { 0.0,        1.20573E-1,  0.0,        0.0,         0.0,        0.0,        -5.93264E-4}
    };

    double s1 = 0;
    for (int i=0; i<6; i++)
    {
        double tau1 = pow(tau - 1, i);
        for (int j=0; j<7; j++)
        {
            if (H1[i][j] != 0.0)
            {
                s1 += H1[i][j]*tau1*pow(delta - 1, j);
            }
        }
    };
    double mu1 = exp(delta*s1);
    
    return muStar*mu0*mu1;
        
}

#endif
