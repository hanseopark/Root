{
	g = new TGraph();

	for(int i = 0; i<100; i++){
		g->SetPoint(i, gRandom->Gaus(5,1), gRandom->Gaus(0,2));
	}

	g->SetMarkerStyle(22);

	g->Draw("APL");
}
