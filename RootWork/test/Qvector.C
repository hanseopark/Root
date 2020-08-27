void Qvector(){

	TCanvas *c1 = new TCanvas("c1","co",1000,500);
	c1->Divide(2,1);

	TRandom3 *trandom = new TRandom3;
	TH1D *h1 = new TH1D("","",100,-pi/2,pi/2);
	TH1D *h2 = new TH1D("","",100,-5,5);
	TH2D *h3 = new TH2D("","",100,-pi/2,pi/2,100,-pi,pi);

	Double_t x= trandom->Rndm();
	Double_t y = trandom->Rndm();
	
	

	Int_t n = 1000;

	for(Int_t k = 0; k<n; k++){

	Double_t qx, qy = 0;
	for(Int_t i = 0; i<n; i++){
	Double_t r = 10*trandom->Rndm();
	Double_t phi = pi*(trandom->Rndm()*2.0-1.0);

	
	qx += r*cos(2*phi);
	qy += r*sin(2*phi);
	
	
	}

	Double_t Phi = atan2(qy,qx)/2;
	Double_t Phi2 = tan(Phi);
	h1 -> Fill(Phi);
	h2 -> Fill(Phi2);
	
	}

	c1->cd(1);
	h1-> Draw();

	c1->cd(2);
	h2-> Draw();
}
