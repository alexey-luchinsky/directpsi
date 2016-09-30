#include <iostream>
#include "TH2D.h"
#include "src/ramC/Random.h"
#include "utils.h"
#include "TH1D.h"
#include "TFile.h"
using namespace std;

double sMin=Mcc*Mcc, sMax=100;
double alpha=3;

double f(double s) {
    return 1/s;
}

int main(void) {
    TH1D *h=new TH1D("h","h",10,sMin,sMax); h->Sumw2();
    int nEv=1e6;
    Random r;
    for(int iEv=0; iEv<nEv; ++iEv) {
        double x=r.rand(0,1);
        double s=sMin+(sMax-sMin)*pow(x,alpha);
        double wt=alpha*(sMax-sMin)*pow(x,alpha-1);
        if(iEv<10)
            cout<<" x="<<x<<" s="<<s<<" wt="<<wt<<endl;
        h->Fill(s,f(s)*wt);
    };
    h->Scale(1./nEv);
    saveHST(h,"h.hst",true);
    
    
    return 0;
}