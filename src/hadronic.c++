#include <iostream>
#include <vector>
#include "ramC/Random.h"
#include "LHAPDF/LHAPDF.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TH2D.h"
using namespace std;
using namespace LHAPDF;

const int nMatr = 26;
TH2D * hMatr[nMatr]; // hMatr[nT, s], nT=-(k1-k3)^2/(s-Mcc)^2 in [0,1])
TFile *hist_file;

string f_to_string(double v) {
    char c[30];
    sprintf(c, "%f", v);
    return string(c);
}

string i_to_string(int v) {
    char c[30];
    sprintf(c, "%d", v);
    return string(c);
}


void load_integrals() {
    hist_file=new TFile("interpolate.root","READ");
    for(int iH=0; iH<nMatr; ++iH) {
        string name="h"+i_to_string(iH);
        hMatr[iH]=(TH2D*)hist_file->Get(name.c_str());
    }
}


int main(void) {
    cout << "hadronic" << endl;
    const string setname = "CT10";
    const int imem = 0;
    const PDF *pdf = mkPDF(setname);
    
    load_integrals();
    
    Random random_generator;
    TFile out_file("PP.root","RECREATE");
    TNtuple tup("tup","tup","hatS:pT2:nT:mtr2:pdf1:pdf2:wt");

    double Mcc = 3.1, Mcc2 = Mcc*Mcc, q2 = Mcc2;
    double S = 700.;
    double sMin=hMatr[0]->GetYaxis()->GetBinLowEdge(1);
    if(sMin<Mcc2) {
        cout<<" root file sMin="<<sMin<<" lower than Mcc2="<<Mcc2<<". Setting sMin=Mcc2"<<endl;
        sMin=Mcc2;
    };
    double sMax=hMatr[0]->GetYaxis()->GetBinUpEdge(hMatr[0]->GetYaxis()->GetNbins());
    if(sMax>S) {
        cout<<" sMax="<<sMax<<" higher than S="<<S<<" setting sMax=S"<<endl;
    };
    int nEv = 1e6;
    for (int iEv = 0; iEv < nEv; ++iEv) {
        double wt = 1;
        double s = random_generator.rand(sMin, sMax);
        wt *= (sMax - sMin)/S; // hat s
        // y
        double yMax = log(S / s) / 2;
        double y = random_generator.rand(-yMax, yMax);
        wt *= 2 * yMax; // y
        double x1 = sqrt(s / S) * exp(y), pdf1 = pdf->xfxQ2(0, x1, q2) / x1;
        double x2 = sqrt(s / S) * exp(-y), pdf2 = pdf->xfxQ2(0, x2, q2)/ x2;
        
        double nT=random_generator.rand(0,1), t=(Mcc2-s)*nT, u=Mcc2-s-t, 
                pT2=t*u/s;
        
        double mtr2=0, mtr;
        if(iEv<10) 
            cout<<" nT="<<nT<<" s="<<s<<" wt="<<wt<<endl;
        for(int iH=1; iH<nMatr; ++iH) {
            mtr=hMatr[iH]->Interpolate(nT,s);
            mtr2 += pow(mtr,2);
        };
        tup.Fill(s,pT2,nT,mtr2,pdf1,pdf2,wt);
    };
    tup.Write(); out_file.Save();
    return 0;
}
