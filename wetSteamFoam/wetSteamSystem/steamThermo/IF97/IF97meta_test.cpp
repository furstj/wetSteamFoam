#include <iostream>

#include "IF97meta.H"

using namespace std;
using namespace IF97::meta;

int main() {

    double p = 1e6;
    double T = 450;

    cout << "v     = "   << v(p,T) << endl;
    cout << "rho   = " << 1/v(p,T) << endl;
    cout << "vr/v0 = "   << vr(p,T)/v0(p,T) << endl;
    cout << endl;
    
    cout << "h     = "   << h(p,T) << endl;
    cout << "hr/h0 = "   << hr(p,T)/h0(p,T) << endl;
    cout << endl;
    
    cout << "u     = "   << u(p,T) << endl;
    cout << "ur/u0 = "   << ur(p,T)/u0(p,T) << endl;
    cout << endl;
    
    cout << "s     = "   << s(p,T) << endl;
    cout << "sr/s0 = "   << sr(p,T)/s0(p,T) << endl;
    cout << endl;
    
    cout << "cp     = "   << cp(p,T) << endl;
    cout << "cpr/cp0= "   << cpr(p,T)/cp0(p,T) << endl;
    cout << endl;
    
    cout << "cv     = "   << cv(p,T) << endl;
    cout << "cvr/cv0= "   << cvr(p,T)/cv0(p,T) << endl;
    cout << endl;

    cout << "w      = "   << sqrt(w2(p,T)) << endl;
    cout << endl;

    return 0;
};
