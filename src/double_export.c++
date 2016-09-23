#include <iostream>
#include "TFile.h"
#include "TNtuple.h"
#include "TH2D.h"

using namespace std;

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

double wave_function(double q2, double delta) {
    const double PI=acos(-1.);
    return exp(-q2/delta/delta)*pow(delta,3)/sqrt(PI);
}

void saveHST2(TH2D *hist, TString name, bool print = false) {
    if (print) cout << " Saving " << name << endl;
    FILE *file = fopen(name.Data(), "w");
    for(int iX=1; iX<=hist->GetNbinsX(); ++iX) {
        for(int iY=1; iY<=hist->GetNbinsY(); ++iY) {
            double x_=hist->GetXaxis()->GetBinCenter(iX);
            double y_=hist->GetYaxis()->GetBinCenter(iY);
            double z_=hist->GetBinContent(iX,iY);
            fprintf(file,"%e %e %e\n",x_,y_,z_);
        }
    }
    fclose(file);
};


int main(void) {
    cout << "double_export" << endl;

    TFile in_file("matr.root", "READ");
    TNtuple *tup = (TNtuple*) in_file.Get("tup");
    double sMin = tup->GetMinimum("s"), sMax = tup->GetMaximum("s");
    int nEv = tup->GetEntries();
    cout << " nEv=" << nEv << " sMin=" << sMin << " sMax=" << sMax << endl;

    TFile out_file("out.root","RECREATE");
    out_file.cd();

    const int nMatr = 26;
    float matr[nMatr];
    TH2D * hMatr[nMatr];
    for (int iH = 0; iH < nMatr; ++iH) {
        string name = "matr_" + i_to_string(iH);
        tup->SetBranchAddress(name.c_str(), &matr[iH]);
        name = "h" + i_to_string(iH);
        hMatr[iH] = new TH2D(name.c_str(), name.c_str(), 50, -1, 1, 10, sMin, sMax);
        hMatr[iH]->Sumw2();
    }
    float weight, wf, cosPsi, q2, s;
    tup->SetBranchAddress("wt", &weight);
    tup->SetBranchAddress("cosPsi", &cosPsi);
    tup->SetBranchAddress("s", &s);

    for (int i = 0; i < nEv; ++i) {
        tup->GetEntry(i);
        double delta=0.4;
        double WF = wave_function(q2, delta);
        for (int iH = 0; iH < nMatr; ++iH)
            hMatr[iH]->Fill(cosPsi, s, matr[iH] * weight * WF);
    }

    for(int i=0; i<nMatr; ++i)
        hMatr[i]->Write();
//    saveHST2(hMatr[1],"h1.hst2");
//    saveHST2(hMatr[2],"h2.hst2");
    
    hMatr[1]->Scale(1./tup->GetEntries());
    saveHST2(hMatr[1],"h1.hst2");

//    for(int iH=2; iH<5; ++iH) {
////        hMatr[iH]->Scale(1./tup->GetEntries());
//        string name="h"+i_to_string(iH)+".hst2";
//        saveHST2(hMatr[iH],name.c_str());
//    }
    out_file.Save();

    return 0;
}
