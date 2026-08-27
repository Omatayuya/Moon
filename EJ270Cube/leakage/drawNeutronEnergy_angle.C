int drawNeutronEnergy_angle(){
	gROOT->Reset();
	gROOT->SetBatch();
	gROOT->ForceStyle();

	/*-----------------------------------*/
	vector<TString> vStrDet{"EJ-270_1", "EJ-270_2"};
	TString strPhi("380MV");
	vector<TString> vStrHydrogen{"H_0ppm", "H_10ppm", "H_20ppm", "H_50ppm", "H_100ppm", "H_200ppm", "H_500ppm", "H_1000ppm", "H_2000ppm", "H_5000ppm", "H_10000ppm"};
	
	constexpr Int_t nTheta = 18;

	vector<TColor*> vColor{myBlue, myRed, myGreen, myCyan, myOrange, myPurple};
	/*-----------------------------------*/

	vector<TCanvas*> vCan;
	vector<TLegend*> vLeg;

	for(size_t i=0; i<vStrHydrogen.size(); ++i){
		stringstream ssIn;
		ssIn << "Neutron_" << strPhi << "_" << vStrHydrogen.at(i) << "_angle.root";
		ifstream ifs(ssIn.str());
		if(ifs.fail()){ cerr << "--- File " << ssIn.str() << " not found" << endl; return 1; }
		auto fin = new TFile(ssIn.str().c_str());

		vCan.push_back(new TCanvas());
		vLeg.push_back(new TLegend(0.4, 0.15, 0.8, 0.46));

		for(Int_t j=0; j<nTheta; ++j){
			stringstream ssName;
			ssName << "neutron_E*Current_" << j;
			auto h = static_cast<TH1D*>(fin->Get(ssName.str().c_str()));

			h->SetStats(0);
			h->SetTitle("Neutron energy ("+vStrHydrogen.at(i)+")");
			h->SetLineColor(49-j);
			h->SetMaximum(0.5e-1);
			h->SetMinimum(1e-6);
			h->Draw("hist same");

			stringstream ssLeg;
			ssLeg << 1.*j*90/18 << "-" << 1.*(j+1)*90/18 << " deg";
			vLeg.back()->AddEntry(h, ssLeg.str().c_str(), "l");
		}
		vLeg.back()->SetNColumns(2);
		vLeg.back()->Draw();
		gPad->SetLogx(); gPad->SetLogy();
	}


	//pdf file
	if(true){
		stringstream ssPdf;
		ssPdf << "NeutronEnergy_angle.pdf";
		TString fPdfOut = ssPdf.str();

		if(vCan.size()==1) vCan.at(0)->Print(fPdfOut);
		else{
			for(size_t i=0; i<vCan.size(); i++){
				if(i==0) vCan.at(i)->Print(fPdfOut+"(");
				else if(i==vCan.size()-1) vCan.at(i)->Print(fPdfOut+")");
				else vCan.at(i)->Print(fPdfOut);
			}
		}
	}

	return 0;
}

