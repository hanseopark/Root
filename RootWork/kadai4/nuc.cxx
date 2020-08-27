{

	Double_t R = 1.2*pow(197,1/3);

	c1 = new TCanvas("c1", "kadai4", 1500,1000);
	c1->Divide(3,2);
	
	h1 = new TH2D("h1","x,y",100,-2,2,100,-2,2);
	h2 = new TH2D("h2","x,z",100,-2,2,100,-2,2);
	h3 = new TH1D("h3","R", 100, 0,3);
	h4 = new TH1D("h4","theta", 100, 0, pi);
	h5 = new TH1D("h5", "phi", 100, 0, 2*pi);
	
	TRandom3 *trandom = new TRandom3;
	Double_t x,y,z = 0;
	Double_t r, theta, phi =0;

for(Int_t i = 0; i <2000; i++){

	x = 20* (trandom->Rndm()*2.0-1.0);
	y = 20* (trandom->Rndm()*2.0-1.0);
	z = 20* (trandom->Rndm()*2.0-1.0);
	
	r = sqrt(x*x+y*y+z*z);
	if(R>r){
		theta = acos(z/r);
		h1 -> Fill(x,y);
		h2 -> Fill(x,z);
		h3 -> Fill(r);
		h4 -> Fill(theta);
		if(x<0 && y<0){
		phi = atan(y/x)+pi;
		}
		if(x<0 && y>0){
		phi = atan(y/x)+pi;
		}
		if(x>0 && y<0){
		phi = atan(y/x)+2*pi;
		}
		if(x>0 && y>0){
		phi = atan(y/x);
		}
		h5 -> Fill(phi);
	//atan2(y,x);
	}
	else{
	i=i-1;
	}
}
	c1->cd(1);
	h1->Draw("");

	c1->cd(2);
	h2->Draw("");

	c1->cd(3);
	h3->Draw();

	c1->cd(4);
	h4->Draw();

	c1->cd(5);
	h5->Draw();
	
}
