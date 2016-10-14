void saveHST(TH1D *hist, TString name, bool print) {
    if (print) cout << " Saving " << name << endl;
    FILE *file = fopen(name.Data(), "w");
    for (int i = 1; i <= hist->GetNbinsX(); ++i) {
        fprintf(file, "%e %e %e\n", hist->GetBinCenter(i), hist->GetBinContent(i) / hist->GetBinWidth(i), hist->GetBinError(i) / hist->GetBinWidth(i));
    };
    if (print) cout << "\t Histogram sum=" << hist->GetSum() << endl;
    fclose(file);
};

void exportQ2() {
	TFile file("matr.root");
	TH1D *h1=new TH1D("h1","h1",20,0,2.2);
	tup->Project("h1","q2","matr_1*wt*(0.3<xs && xs<0.31 && -0.1<cosPsi && cosPsi<0.1)");
	saveHST(h1,"hQ2_1.hst",true);	
}
