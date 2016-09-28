#include <iostream>
#include "TH2D.h"
#include "src/ramC/Random.h"
using namespace std;

double f(double x, double y) {
    return 1/exp(x) + y*y;
}

double xMin = -2, xMax = 2, yMin = 0, yMax = 1;
int nbx = 30, nby = 20;
TH2D *h;

void check() {
    Random r(0);
    h->Reset();
    int nEv = 1e7;
    for (int i = 0; i < nEv; ++i) {
        double x = r.rand(xMin, xMax), y = r.rand(yMin, yMax);
        h->Fill(x, y, f(x, y));
    }
    h->Scale(1. / nEv);

    //    double wx=h->GetXaxis()->GetBinWidth(1), wy=h->GetYaxis()->GetBinWidth(1);
    //    cout<<" wx="<<wx<<" wy="<<wy<<endl;
    double x0 = r.rand(xMin, xMax), y0 = r.rand(yMin, yMax);
    cout << h->Interpolate(x0, y0) * nbx * nby / f(x0, y0) << endl;
}

int main(void) {
    h = new TH2D("h", "h", nbx, xMin, xMax, nby, yMin, yMax);
    for (int i = 0; i < 10; ++i)
        check();
    return 0;
}
