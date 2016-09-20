#include <iostream>
#include "TFile.h"
#include "TNtuple.h"
#include "TH1D.h"
#include <tclap/CmdLine.h>


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

double wave_function(double q2, double delta) {
    const double PI=acos(-1.);
    return exp(-q2/delta/delta)*pow(delta,3)/sqrt(PI);
}

string in_name;
string out_name;
float delta;
bool init_args(int argc, char **argv) {
    cout<<"init_args"<<endl;
    try {
        TCLAP::CmdLine cmd("Extracts histograms from root file",' ',"0.9");
        TCLAP::ValueArg<string> in_cmd("i","in","input file name",false,"matr.root","string",cmd);
        TCLAP::ValueArg<string> out_cmd("o","out","prefix for out file names",false,"matr_","string",cmd);
        TCLAP::ValueArg<float> delta_cmd("d","delta","wave function parameter",false,0.4,"float",cmd);
        cmd.parse(argc, argv);
        in_name=in_cmd.getValue();
        out_name=out_cmd.getValue();
        delta=delta_cmd.getValue();
    } catch (TCLAP::ArgException e) {
        cout<<" ArgException "<<e.error()<<" in argument "<<e.argId()<<endl;
        return false;
    }
    return true;
}

int main(int argc, char **argv) {
    if(!init_args(argc, argv)) return -1;
    cout<<" in_name="<<in_name<<endl;
    cout<<" out_name="<<out_name<<endl;
    cout<<" delta="<<delta<<endl;
    TFile file(in_name.c_str(),"READ");
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
        double WF=wave_function(q2,delta);
        for(int iH=0; iH<nMatr; ++iH)
            hMatr[iH]->Fill(cosPsi,matr[iH]*weight*WF);
    }
    for(int iH=0; iH<nMatr; ++iH) {
        hMatr[iH]->Scale(1./tup->GetEntries());
        string name=out_name+i_to_string(iH)+".hst";
        saveHST(hMatr[iH],name.c_str());
    }
}
