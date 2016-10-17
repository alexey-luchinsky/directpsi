void gen_bins(double min, double max, int nBins, double alpha, string name="bins.txt") {
	cout<<"Saving bins to "<<name<<endl;
	ofstream file(name.c_str());
	
	for(double x=0; x<=1.; x +=1./nBins) {
		file<<min+(max-min)*pow(x,alpha)<<endl;
		cout<<min+(max-min)*pow(x,alpha)<<endl;
	};
	
	file.close();
	gApplication->Terminate();
}