#include <iostream>

#include "IF97.H"

using namespace std;
using namespace IF97;

int main() {

    double T[] = {298.15, 298.15, 373.15, 433.14, 433.15, 873.15, 873.15, 873.15,
                  1173.15, 1173.15, 1173.15};
    double rho [] = {998, 1200, 1000, 1, 1000, 1, 100, 600, 1, 100, 400};

    for (int i=0; i<11; i++)
    {
        cout << T[i] << "\t" << rho[i] << "\t" << mu(rho[i],T[i]) << "\n";
    }
    return 0;
};
