
void circle(){
	TCanvas *c1 = new TCanvas("c1","entry_10000",1200,800);
	c1 -> Divide(2,1);

	TH1D *h2 = new TH1D("","",100,2.8,3.5);
	TH1D *h4 = new TH1D("","",100,2.8,3.5);

	h2 -> SetTitle("1000");
	h4 -> SetTitle("10000");

	TRandom3 *trandom = new TRandom3;
	Double_t sum1,sum2,sum3,sum4 = 0;
	Double_t area1,area2,area3,area4 = 0;

	for(Int_t k = 0; k<10000; k++){
		
		sum2 = 0;
		area2 = 0;
		sum4 =0;
		area4 = 0;
		
		for(Int_t i = 1; i<=1000; i++){
		Double_t x2= trandom->Rndm()*2.0-1.0;
		Double_t y2= trandom->Rndm()*2.0-1.0;
		Double_t D2 = pow(x2,2)+pow(y2,2);

		if(D2<1){
		
			sum2 = sum2 + 1;
			area2 = sum2/i*4;
		

		}
		}

		h2 -> Fill(area2);

	for(Int_t i = 1; i<=10000; i++){
		Double_t x4= trandom->Rndm()*2.0-1.0;
		Double_t y4= trandom->Rndm()*2.0-1.0;
		Double_t D4 = pow(x4,2)+pow(y4,2);

		if(D4<1){
			
			sum4 = sum4 + 1;
			area4 = sum4/i*4;

		}
		}
		h4 -> Fill(area4);
	}		
	c1 -> cd(1);
	h2 -> Draw();
	
	c1 -> cd(2);
	h4 -> Draw();
}
