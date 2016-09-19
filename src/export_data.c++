#include <iostream>
#include "TFile.h"
#include "TNtuple.h"
#include "TH1D.h"

using namespace std;

void saveHST(TH1D *hist, TString name, bool print = false) {
    if (print) cout << " Saving " << name << endl;
    FILE *file = fopen(name.Data(), "w");
    for (int i = 1; i <= hist->GetNbinsX(); ++i) {
        fprintf(file, "%e %e %e\n", hist->GetBinCenter(i), hist->GetBinContent(i) / hist->GetBinWidth(i), hist->GetBinError(i) / hist->GetBinWidth(i));
    };
    if (print) cout << "\t Histogram sum=" << hist->GetSum() << endl;
    fclose(file);
};

string f_to_string(double v) {
  char c[20];
  sprintf(c,"%f",v);
  return string(c);
}

string i_to_string(int v) {
  char c[20];
  sprintf(c,"%d",v);
  return string(c);
}


int main(int argc, char **argv) {
    cout<<"Hello, world!"<<endl;
    if(argc<2) {
        cout<<"data file required!"<<endl;
        return -1;
    }
    TFile file(argv[1],"READ");
    TNtuple *tup=(TNtuple*)file.Get("tup");

    const int nMatr = 26;
    float matr[nMatr];
    TH1D *hMatr[nMatr];
    for(int iH=0; iH<nMatr; ++iH) {
        string name="matr_"+i_to_string(iH);
        tup->SetBranchAddress(name.c_str(), &matr[iH]);
        name="h"+i_to_string(iH);
        hMatr[iH]=new TH1D(name.c_str(), name.c_str(), 20, -1, -1);
        hMatr[iH]->Sumw2();
    }


    
    float weight, wf, cosPsi, q2;
    tup->SetBranchAddress("wt",&weight);
    tup->SetBranchAddress("cosPsi",&cosPsi);
    for(int i=0; i<tup->GetEntries(); ++i) {
        tup->GetEntry(i);
        for(int iH=0; iH<nMatr; ++iH)
            hMatr[iH]->Fill(cosPsi,matr[iH]*weight);
    }
    for(int iH=0; iH<nMatr; ++iH) {
        hMatr[iH]->Scale(1./tup->GetEntries());
        string name="matr_"+i_to_string(iH)+".hst";
        saveHST(hMatr[iH],name.c_str());
    }
}
