{
	c1 = new TCanvas("c1", "My root plots", 1000, 500);
	c1 ->Divide(2,1);

	c2 = new TCanvas("c2", "The other plots", 500, 500);
	c2 ->Divide(2,1);

	mc(0);

	c1->cd(1);
 	f1 = new TF1("f1", "sin(x)", 0, 5);
	f1 ->SetTitle("My Function");
	f1 ->SetLineColor(kBlue);
	f1 ->Draw();

	c1 ->cd(2);
	f2 = new TF1("f2", "cos(x)", 0, 5);
	f2 ->SetTitle("f2 function");
	f2 ->SetLineColor(kRed);
	f2 ->Draw();

	c2 ->cd(1);
	f3 = new TF2("f3", "x*x+y*y", -10, 10, -10, 10);
	f3 ->SetTitle("f3 function");
	f3 ->Draw("surf");


}

