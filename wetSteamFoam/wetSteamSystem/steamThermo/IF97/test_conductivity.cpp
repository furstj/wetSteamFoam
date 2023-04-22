#include <iostream>

#include "IF97.H"

using namespace std;
using namespace IF97;

int main() {

    double T[] = {100, 200, 300, 400};
    double p[] = {1e5, 1e5, 1e5, 1e5};

    for (int i=0; i<4; i++)
    {
        double rho = 1/reg2meta::v(p[i],T[i]+273.15);
        cout << T[i] << "\t" << p[i] << "\t" << lambda(rho,T[i]+273.15) << "\n";
    }
    return 0;
};
