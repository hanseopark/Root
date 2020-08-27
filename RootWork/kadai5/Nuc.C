#include<math.h>

double rho(double r){
	double A_Pb = 207.2;
	double R = 1.18*pow(A_Pb,1/3.)-0.48;
	double a = 0.54;
	double rho;
	rho = 1/(1+exp((r-R)/a));
	return rho;
}

void Nuc(){

	TCanvas *c1 = new TCanvas("c1", "nuc", 1500,1000);
	TCanvas *c2 = new TCanvas("c2","c2", 1500,500);
	c1->Divide(3,2);
	c2->Divide(3,1);
	
	TH3D *h1 = new TH3D("h1", "", 100 ,-10,10, 100 ,-10,10, 100 ,-10,10);
	TH1D *h2 = new TH1D("h2", "r",100 , 0, 10);
	TH1D *h3 = new TH1D("h3", "theta", 100, 0 , 4);
	TH1D *h4 = new TH1D("h4", "phi", 100, 0, 8);
	TH1D *h5 = new TH1D("h5", "wounded_nucleon_b=0",100, 0, 1000);
	double A_Pb=207.2;
	double R = 1.18*pow(A_Pb,1/3.)-0.48; //fm
	double a = 0.54; //fm

	double x2,y2,z2 = 0;
	double r1,theta, phi =0;

	double b_max = sqrt(4.2/pi); //mb
	int count1 = 0;
	int count2 = 0;
	int count3 = 0;

	//first
	for(int i =1; i< A_Pb*10; i++){
 	 x2 = R* (gRandom->Rndm()*2.0-1.0);
         y2 = R* (gRandom->Rndm()*2.0-1.0);
         z2 = R* (gRandom->Rndm()*2.0-1.0);
 
         r1 = sqrt(x2*x2+y2*y2+z2*z2);
 	 theta = acos(z2/r1);
	 if(x2<0 && y2<0){
	phi = atan(y2/x2)+pi;
	}
	if(x2<0 && y2>0){
	phi = atan(y2/x2)+pi;
	}
	if(x2>0 && y2<0){
	phi = atan(y2/x2)+2*pi;
	}
	if(x2>0 && y2>0){
	phi = atan(y2/x2);
	}

	 	if(rho(r1)>gRandom->Rndm()){
		h1 ->Fill(r1*sin(theta)*cos(phi),r1*sin(theta)*sin(phi),r1*cos(theta));
		h2 ->Fill(r1);
		h3 ->Fill(theta);
		h4 ->Fill(phi);
		}
		else{
		i= i-1;
		}

	
	}

	
	//second
	Double_t X[1000];
	Double_t Y1[1000];
	Double_t Y2[1000];
	Double_t Y3[1000];
	for(Int_t k= 1; k < 1001; k++){
	for(int i =1; i< A_Pb; i++){
 	 x2 = R* (gRandom->Rndm()*2.0-1.0);
         y2 = R* (gRandom->Rndm()*2.0-1.0);
         z2 = R* (gRandom->Rndm()*2.0-1.0);
 
         r1 = sqrt(x2*x2+y2*y2+z2*z2);
 	 theta = acos(z2/r1);
	 if(x2<0 && y2<0){
	phi = atan(y2/x2)+pi;
	}
	if(x2<0 && y2>0){
	phi = atan(y2/x2)+pi;
	}
	if(x2>0 && y2<0){
	phi = atan(y2/x2)+2*pi;
	}
	if(x2>0 && y2>0){
	phi = atan(y2/x2);
	}

	
	 	if(rho(r1)>gRandom->Rndm()){
		double D_0 = sqrt(pow(r1*sin(theta)*cos(phi)-0,2)+pow(r1*sin(theta)*sin(phi),2));
		double D_3 = sqrt(pow(r1*sin(theta)*cos(phi)-3,2)+pow(r1*sin(theta)*sin(phi),2));
		double D_6 = sqrt(pow(r1*sin(theta)*cos(phi)-6,2)+pow(r1*sin(theta)*sin(phi),2));
		if(D_0<b_max){
		count1 = count1+1;
		}
		
		if(D_3<b_max){
		count2 = count2+1;
		}
		
		if(D_6<b_max){
		count3 = count3+1;
		}

		}
		else{
		i= i-1;
		}

	
	}
	X[k] = k;
	Y1[k] = count1;
	Y2[k] = count2;
	Y3[k] = count3;
	count1 = 0;
	count2 = 0;
	count3 = 0;
	}
	
	TGraph *gimp0 = new TGraph(1000,X,Y1);
	gimp0 -> SetTitle("Impact parameter = 0");
	TGraph *gimp3 = new TGraph(1000,X,Y2);
	gimp3 -> SetTitle("Impact parameter = 3");
	TGraph *gimp6 = new TGraph(1000,X,Y3);
	gimp6 -> SetTitle("Impact parameter = 6");

	TF1* f1 = new TF1("f1", "rho(x)",0,10);
	
	c1 ->cd(1);
	f1 ->Draw();

	c1 ->cd(2);
	h1 ->Draw();
	
	c1 ->cd(3);
	h2 ->Draw();

	c1 ->cd(4);
	h3 ->Draw();

	c1 ->cd(5);
	h4 ->Draw();

	c2 ->cd(1);
	gimp0 ->Draw();
	
	c2 ->cd(2);
	gimp3 ->Draw();
	
	c2 ->cd(3);
	gimp6 ->Draw();
}

