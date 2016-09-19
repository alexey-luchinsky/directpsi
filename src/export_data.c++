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
    TH1D *hist=new TH1D("h","h",20,-1,1);
    hist->Sumw2();
    
    tup->Project("h","cosPsi","matr_1*wt*wf");
    hist->Scale(1./tup->GetEntries());
    saveHST(hist,"h1.hst",true);
    
    hist->Reset();
    tup->Project("h","cosPsi","matr_2*wt*wf");
    hist->Scale(1./tup->GetEntries());
    saveHST(hist,"h2.hst",true);
}
