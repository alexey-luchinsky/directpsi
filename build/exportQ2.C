void saveHST(TH1D *hist, TString name, bool print) {
    if (print) cout << " Saving " << name << endl;
    FILE *file = fopen(name.Data(), "w");
    for (int i = 1; i <= hist->GetNbinsX(); ++i) {
        fprintf(file, "%e %e %e %e\n", 
        	hist->GetBinCenter(i), 
        	hist->GetBinContent(i) / hist->GetBinWidth(i), 
        	hist->GetBinWidth(i)/2,
        	hist->GetBinError(i) / hist->GetBinWidth(i));
    };
    if (print) cout << "\t Histogram sum=" << hist->GetSum() << endl;
    fclose(file);
};

bool exportQ2() {
	TFile *file=new TFile("matr.root");
	if(!file->IsOpen()) return false;
	TH1D *h1=new TH1D("h1","h1",20,0,2.2); h1->Sumw2();
	tup->Project("h1","q2","matr_1*wt*(0.3<xs && xs<0.31 && -0.1<cosPsi && cosPsi<0.1)");
	saveHST(h1,"partonic/hQ2_1.hst",true);
	h1->Reset();
	tup->Project("h1","q2","matr_15*wt*(0.3<xs && xs<0.31 && -0.1<cosPsi && cosPsi<0.1)");
	saveHST(h1,"partonic/hQ2_15.hst",true);
	return true;
}

string i_to_string(int v) {
    char c[30];
    sprintf(c, "%d", v);
    return string(c);
}

const int nMatr = 26;
TH2D * hMatr[nMatr]; // hMatr[nT, xs], nT=-(k1-k3)^2/(s-Mcc)^2 in [0,1])
float sMin, sMax, alpha, delta;

bool load_integrals(char *in_fileName) {
    hist_file = new TFile(in_fileName, "READ");
    if (!hist_file->IsOpen()) return false;

    // read statistics
    TVectorD stats = *((TVectorD*) hist_file->Get("stats"));
    sMin = stats(0);
    sMax = stats(1);
    alpha = stats(2);
    delta = stats(3);

    for (int iH = 0; iH < nMatr; ++iH) {
        string name = "h" + i_to_string(iH);
        hMatr[iH] = (TH2D*) hist_file->Get(name.c_str());
    }
	cout<<" sMin="<<sMin<<" sMax="<<sMax<<" alpha="<<alpha<<" delta="<<delta<<endl;
    return true;
}

double Mcc=3.1;

bool exportEN(char *in_fileName, char *out_name) {
	if(!load_integrals(in_fileName)) {
		return false;
	}
	TH1D *hEn=new TH1D("hEn","hEn",20,3.2,5.5); hEn->Sumw2();
	int nEv=5e6;
	TRandom rand;
	for(int iEv=0; iEv<nEv; ++iEv) {
		double wt=1;
		double xs=rand.Rndm();
		double s=sMin+(sMax-sMin)*pow(xs,alpha); wt *= alpha * (sMax - sMin) * pow(xs, alpha - 1);
		double nT=rand.Rndm(); wt *= (s-Mcc*Mcc);
        // conversion to dsdt
        double PI = acos(-1.);
        wt *= 1. / (64 * PI * s)*4 / s;
        // symmetry, etc
        wt *= 1. / (2 * 2 * 8 * 2 * 8);
        // tranfer to nb
        wt *= 0.389*1e6;

        double mtr2 = 0, mtr;
        for (int iH = 1; iH < nMatr; ++iH) {
            double alphaQCD = 0.3;
            mtr = hMatr[iH]->Interpolate(nT, xs);
            mtr2 += pow(4 * PI*alphaQCD, 3) * pow(mtr, 2);
        };
	
		hEn->Fill(sqrt(s),wt*mtr2);
	};
	hEn->Scale(1./nEv);
	saveHST(hEn,out_name,true);
}

void all() {
	if(!exportQ2())
		cout<<"Cannot open file matr.root, try run calcInts.exe");
	if(!exportEN("inter001.root","partonic/hEn_001.hst"))
		cout<<"Cannot open file inter001.root. Try run ./double_export.exe -s 20 -d 0.01 -o inter001.root";
	if(!exportEN("inter04.root","partonic/hEn_04.hst"))
		cout<<"Cannot open file inter04.root. Try run ./double_export.exe -s 20 -d 0.4 -o inter04.root";
};