{
	g = new TGraph("data.txt");
	g->SetMarkerStyle(21);
	g->SetLineColor(kRed);
	g->SetTitle("Data Plot with TGRAPH;gw;ms");
	g->Draw("ALP");
}
