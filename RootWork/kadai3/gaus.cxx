Double_t func(double x){
	Double_t Func = 0;
	Func = exp(-x*x/2);

	return Func;
}

Double_t integral(double start, double end){
	int num = 1000;
	double ans = 0;
	double xi[1000] = {0,};
	double dx = (end-start)/num;

	for(int i = 0; i<num; i++){
		
		xi[i] = start + i*dx;
		ans += func(xi[i])*dx;
	}

	return ans;
	
}

Double_t Func(double x){
	Double_t Func = 0;
	Func = integral(-10,x)/integral(-10,10);

	return Func;
}

void gaus()
{
	TRandom3 *trandom = new TRandom3;
	TCanvas *c1 = new TCanvas("c1","",1200,800);
	
	TH1D *h1 = new TH1D("h1","gaus(Inverse)", 200, -10 ,10);
	TH1D *h2 = new TH1D("h2","gaus(Rejection)", 200, -10 ,10);
	TH1D *h3 = new TH1D("h3","gaus(Inverse2)", 200, -10 ,10);
	TH1D *h4 = new TH1D("h4","gaus(test)", 200, -10 ,10);
	
	

	float xx[200],yy[200];

	for(Int_t i=0; i<200; i++){
		Double_t k= -10 + i*0.1;
		Double_t x =0;
		Double_t y = trandom->Rndm();
		xx[i] = k;
		yy[i]= Func(k);
		h3->Fill(k,func(k));
	}
	TGraph *g1 = new TGraph(200,yy,xx);

	for(int i =0;i<10000; i++){
	Double_t x2,y2=0;
	Double_t y = trandom->Rndm();
	Double_t x1 = g1 -> Eval(y); //star

	x2 = 10*(trandom->Rndm()*2.0-1.0);
	y2 = trandom->Rndm();
	
	h1 ->Fill(x1);
	Double_t x3 = h3->GetRandom();
	h4 -> Fill(x3);

	if(y2 < exp(-x2*x2/2)){
	h2 ->Fill(x2);
	}
	}
	c1->Divide(2,1);
	
	c1->cd(1);
	g1->Draw("AP");

	h1 ->SetLineColor(1);
	h2 ->SetLineColor(2);
	h3 ->SetLineColor(3);
	h4 ->SetLineColor(4);
	
	c1->cd(2);
	h1->Draw();
	h2 ->Draw("same");
	h3 ->Draw("same");
	h4 ->Draw("same");
}
