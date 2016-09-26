#include <iostream>
#include <vector>
using namespace std;

#include "LHAPDF/LHAPDF.h"
using namespace LHAPDF;
using namespace std;

int main(void) {
    cout << "hadronic" << endl;
    const string setname="CT10";
    const int imem=0;
    const PDF *pdf=mkPDF(setname);
    return 0;
}
