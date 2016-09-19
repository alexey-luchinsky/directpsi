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
    cout<<tup->GetEntries()<<" entries "<<endl;
    TH1D *hist0=new TH1D("h0","h0",20,-1,1);  hist0->Sumw2();
    TH1D *hist10=new TH1D("h10","h10",20,-1,1);  hist10->Sumw2();
    
    float weight, wf, cosPsi, q2, matr0, matr10;
    tup->SetBranchAddress("wt",&weight);
    tup->SetBranchAddress("cosPsi",&cosPsi);
    tup->SetBranchAddress("matr_0",&matr0);
    tup->SetBranchAddress("matr_10",&matr10);
    for(int i=0; i<tup->GetEntries(); ++i) {
        tup->GetEntry(i);
//        if(i<10) 
//            cout<<"i="<<i<<" wt="<<weight<<" wf="<<wf<<" matr0="<<matr0<<endl;
        hist0->Fill(cosPsi,matr0*weight);
        hist10->Fill(cosPsi,matr10*weight);
    }
    hist0->Scale(1./tup->GetEntries());  saveHST(hist0,"matr0.hst",true);
    hist10->Scale(1./tup->GetEntries());  saveHST(hist10,"matr10.hst",true);
}
