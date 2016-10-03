#include <iostream>
#include "mtr0.h"
using namespace std;

double mc, s, T, U;
double NG, NJ0, NJ1, NJ2;

void fill_matr() {
    NG = 1. / sqrt(s * T * U);
    NJ0 = 2. / sqrt(s * T * U);
    NJ1 = 1. / sqrt(s * T * U);
    NJ2 = 1. / Mcc;
};

int main(void) {
    cout<<"delta"<<endl;
    return 0;
}
