#ifndef IF97_reg2meta_H
#define IF97_reg2meta_H

namespace reg2meta
{

constexpr double pStar   = 1e6;       // [Pa]
constexpr double TStar   = 540;       // [K]

constexpr double n0[] = {
    -0.96937268393049E+01,  0.10087275970006E+02, -0.56087911283020E-02,
    +0.71452738081455E-01, -0.40710498223928E+00,  0.14240819171444E+01,
    -0.43839511319450E+01, -0.28408632460772E+00,  0.21268463753307E-01
};

constexpr int J0[] = {0, 1, -5, -4, -3, -2, -1, 2, 3};

constexpr double nr[] = {
    -0.73362260186506E-02, -0.88223831943146E-01, -0.72334555213245E-01,
    -0.40813178534455E-02,  0.20097803380207E-02, -0.53045921898642E-01,
    -0.76190409086970E-02, -0.63498037657313E-02, -0.86043093028588E-01,
    +0.75321581522770E-02, -0.79238375446139E-02, -0.22888160778447E-03,
    -0.26456501482810E-02 }; 

constexpr int Ir[] = {1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5, 5};

constexpr int Jr[] = {0, 2, 5, 11, 1, 7, 16, 4, 16, 7, 10, 9, 10};

// Ideal part of Gibbs energy
inline double gamma0(double p, double T)
{
    double pi = p/pStar;
    double tau = TStar/T;
    double gamma = log(pi);
    for (auto i=0; i<9; i++)
        gamma += n0[i]*pow(tau, J0[i]);
    return gamma;
}

inline double gamma0_pi(double p, double T)
{
    return pStar/p;
}

inline double gamma0_pipi(double p, double T)
{
    return -pStar*pStar/(p*p);
}

inline double gamma0_tau(double p, double T)
{
    double tau = TStar/T;
    double gamma_tau = 0;
    for (auto i=0; i<9; i++)
        gamma_tau += n0[i]*J0[i]*pow(tau, J0[i]-1);
    return gamma_tau;
}

inline double gamma0_tautau(double p, double T)
{
    double tau = TStar/T;
    double gamma_tautau = 0;
    for (auto i=0; i<9; i++)
        gamma_tautau += n0[i]*J0[i]*(J0[i]-1)*pow(tau, J0[i]-2);
    return gamma_tautau;
}

inline double gamma0_tautautau(double p, double T)
{
    double tau = TStar/T;
    double gamma_tautautau = 0;
    for (auto i=0; i<9; i++)
        gamma_tautautau += n0[i]*J0[i]*(J0[i]-1)*(J0[i]-2)*pow(tau, J0[i]-3);
    return gamma_tautautau;
}

// Residual part of Gibbs energy
inline double gammar(double p, double T)
{
    double pi = p/pStar;
    double tau5 = TStar/T - 0.5;
    double gamma = 0;
    for (auto i=0; i<13; i++)
        gamma += nr[i]*pow(pi,Ir[i])*pow(tau5, Jr[i]);
    return gamma;
}

inline double gammar_pi(double p, double T)
{
    double pi = p/pStar;
    double tau5 = TStar/T - 0.5;
    double gamma = 0;
    for (auto i=0; i<13; i++)
        gamma += nr[i]*Ir[i]*pow(pi,Ir[i]-1)*pow(tau5, Jr[i]);
    return gamma;
}

inline double gammar_pipi(double p, double T)
{
    double pi = p/pStar;
    double tau5 = TStar/T - 0.5;
    double gamma = 0;
    for (auto i=0; i<13; i++)
        gamma += nr[i]*Ir[i]*(Ir[i]-1)*pow(pi,Ir[i]-2)*pow(tau5, Jr[i]);
    return gamma;
}

inline double gammar_tau(double p, double T)
{
    double pi = p/pStar;
    double tau5 = TStar/T - 0.5;
    double gamma = 0;
    for (auto i=0; i<13; i++)
        gamma += nr[i]*pow(pi,Ir[i])*Jr[i]*pow(tau5, Jr[i]-1);
    return gamma;
}

inline double gammar_tautau(double p, double T)
{
    double pi = p/pStar;
    double tau5 = TStar/T - 0.5;
    double gamma = 0;
    for (auto i=0; i<13; i++)
        gamma += nr[i]*pow(pi,Ir[i])*Jr[i]*(Jr[i]-1)*pow(tau5, Jr[i]-2);
    return gamma;
}

inline double gammar_pitau(double p, double T)
{
    double pi = p/pStar;
    double tau5 = TStar/T - 0.5;
    double gamma = 0;
    for (auto i=0; i<13; i++)
        gamma += nr[i]*Ir[i]*pow(pi,Ir[i]-1)*Jr[i]*pow(tau5, Jr[i]-1);
    return gamma;
}

inline double gammar_tautautau(double p, double T)
{
    double pi = p/pStar;
    double tau5 = TStar/T - 0.5;
    double gamma = 0;
    for (auto i=0; i<13; i++)
        gamma += nr[i]*pow(pi,Ir[i])*Jr[i]*(Jr[i]-1)*(Jr[i]-2)*pow(tau5, Jr[i]-3);
    return gamma;
}

//======================================================================


// Specific volume [m3/kg]
inline double v0(double p, double T)
{
    return R*T/p;
}


inline double vr(double p, double T)
{
    return R*T/pStar*gammar_pi(p,T);
}

inline double v(double p, double T)
{
    return v0(p,T) + vr(p,T);
}

// Compression factor [-] (z=pv/RT)
inline double z(double p, double T)
{
    return 1.0 + p/pStar*gammar_pi(p,T);
}

// Specific enthalpy [J/kg]
inline double h0(double p, double T)
{
    return R*TStar*gamma0_tau(p,T);
}


inline double hr(double p, double T)
{
    return R*TStar*gammar_tau(p,T);
}

inline double h(double p, double T)
{
    return h0(p,T) + hr(p,T);
}


// Specific internal energy [J/kg]
inline double u0(double p, double T)
{
    double tau = TStar/T;
    double pi = p/pStar;
    return R*T*(tau*gamma0_tau(p,T) - pi*gamma0_pi(p,T));
}


inline double ur(double p, double T)
{
    double tau = TStar/T;
    double pi = p/pStar;
    return R*T*(tau*gammar_tau(p,T) - pi*gammar_pi(p,T));
}

inline double u(double p, double T)
{
    return u0(p,T) + ur(p,T);
}


// Specific entropy [J/kg/K]
inline double s0(double p, double T)
{
    double tau = TStar/T;
    return R*(tau*gamma0_tau(p,T) - gamma0(p,T));
}


inline double sr(double p, double T)
{
    double tau = TStar/T;
    return R*(tau*gammar_tau(p,T) - gammar(p,T));
}

inline double s(double p, double T)
{
    return s0(p,T) + sr(p,T);
}


// Specific heat capacity at const. pressure [J/kg/K]
inline double cp0(double p, double T)
{
    double tau2 = TStar*TStar/(T*T);
    return -R*tau2*gamma0_tautau(p,T);
}


inline double cpr(double p, double T)
{
    double tau2 = TStar*TStar/(T*T);
    return -R*tau2*gammar_tautau(p,T);
}

inline double cp(double p, double T)
{
    return cp0(p,T) + cpr(p,T);
}

inline double dcpdT(double p, double T)
{
    double tau = TStar/T;
    return R*(2*tau*(gamma0_tautau(p,T) + gammar_tautau(p,T)) +
    tau*tau*(gamma0_tautautau(p,T) + gammar_tautautau(p,T)))*tau/T;
}


// Specific heat capacity at const. volume [J/kg/K]
inline double cv0(double p, double T)
{
    double tau = TStar/T;
    return -R*tau*tau*gamma0_tautau(p,T);
}


inline double cvr(double p, double T)
{
    double tau = TStar/T;
    double pi  = p/pStar;
    return -R*(tau*tau*gammar_tautau(p,T) +
    pow(1 + pi*gammar_pi(p,T) - tau*pi*gammar_pitau(p,T), 2)/(1 - pi*pi*gammar_pipi(p,T)));
}

inline double cv(double p, double T)
{
    return cv0(p,T) + cvr(p,T);
}


inline double w2(double p, double T)
{
    double tau = TStar/T;
    double pi  = p/pStar;
    double grp = gammar_pi(p,T);
    double grpp = gammar_pipi(p,T);
    double grpt = gammar_pitau(p,T);
    double grtt = gammar_tautau(p,T);
    double g0tt = gamma0_tautau(p,T);

    return R*T*(1 + 2*pi*grp + pi*pi*grp*grp)/
        ((1 - pi*pi*grpp) + pow(1 + pi*grp - tau*pi*grpt, 2)/(tau*tau*(g0tt + grtt)));
}


}

#endif
