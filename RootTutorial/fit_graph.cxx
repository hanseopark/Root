{
	g = new TGraph();

	double x,y;

	for(int i=0; i<10; i++){
	x= 0.5*i;
	y= 4*x*x+exp(-x/3) + gRandom->Gaus(0,1);
	g->SetPoint(i, x, y);
	}

	f = new TF1("f","[0]*x^2 + [1]*exp(-x*[2])",0 ,5);
	f -> SetParNames("quar","amp","decay");
	f -> SetParameters(2, 1);
	
	g->Fit(f);


	g->SetMarkerStyle(22);

	g->Draw("AP");	

}
