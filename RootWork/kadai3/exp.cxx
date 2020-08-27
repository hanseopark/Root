{
	TRandom3 *trandom = new TRandom3;
	Double_t x,y,f,g = 0;
	Int_t slope = 0; //0~inf
	
	Double_t x1, y2 =0;
	cout << "slope: ?" <<endl;
	cin >>slope;

	h1 = new TH1D("h1", "Uniform",100,0,1000);
	h2 = new TH1D("h2", "Exp(Inverse function)",100,0,100);
	h3 = new TH1D("h3", "Exp(Rejection)",100,0,100);

	for(Int_t i=0; i<10000; i++){
	x = 1000* trandom->Rndm();
	y = trandom->Rndm();  //uniform random number 0~1
	f = exp(-x/slope); //exp function
	g = slope*y;

	x1 = 1000* trandom->Rndm();
	y2 = trandom->Rndm();
	
	h1 ->Fill(x);
	x = -slope*log(1-y);

	h2 ->Fill(x);

	if(y2 < exp(-x1/slope)){
	h3 ->Fill(x1);
	}

	}

	c1 = new TCanvas("c1", "exp function", 1200, 800);
	c1 ->Divide(2,1);
	

/*	c1->cd(1);
	h1 -> SetMinimum(0);
	h1 -> Draw();
*/
	c1->cd(1);
	h2 -> SetMinimum(0);
	h2 -> Draw();

	c1->cd(2);
	h3 -> SetMinimum(0);
	h3 -> Draw();

}
