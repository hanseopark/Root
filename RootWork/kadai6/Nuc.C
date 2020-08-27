#include<math.h>

double rho(double r){
	double A_Pb = 207.2;
	double R = 1.18*pow(A_Pb,1/3)-0.48;
	double a = 0.45;
	double rho;
	rho = 1/(1+exp((r-R)/a));
	return rho;
}

void Nuc(){

	TCanvas *c1 = new TCanvas("c1", "nuc", 1500,1000);
	c1->Divide(3,2);
	
	TH3D *h1 = new TH3D("h1", "", 100 ,-5,5, 100 ,-5,5, 100 ,-5,5);
	TH1D *h2 = new TH1D("h2", "r",100 , 0, 3);
	TH1D *h3 = new TH1D("h3", "theta", 100, 0 , 10);
	TH1D *h4 = new TH1D("h4", "phi", 100, 0, 10);
	double A_Pb=207.2;
	double R = 1.18*pow(A_Pb,1/3)-0.48; //fm
	double a = 0.45; ///fm

	double x2,y2,z2 = 0;
	double r1,theta, phi =0;

	double b_max = sqrt(42/pi); //mb
	int count1 = 0;
	int count2 = 0;
	int count3 = 0;

	for(int i =0; i< A_Pb; i++){
 	 x2 = 5* (gRandom->Rndm()*2.0-1.0);
         y2 = 5* (gRandom->Rndm()*2.0-1.0);
         z2 = 5* (gRandom->Rndm()*2.0-1.0);
 
         r1 = sqrt(x2*x2+y2*y2+z2*z2);
         theta = pi*gRandom->Rndm();
         phi = 2*pi*gRandom->Rndm();
 		
	 	if(rho(r1)>gRandom->Rndm()){	
		h1 ->Fill(x2,y2,z2);
		h2 ->Fill(r1);
		h3 ->Fill(theta);
		h4 ->Fill(phi);
		
		if(r1<b_max){
		count1 = count1+1;
		}
		
		if(r1<b_max-3){
		count2 = count2+1;
		}
		
		if(r1<b_max-6){
		count3 = count3+1;
		}

		}
		else{
		i= i-1;
		}

	
	}
	
	cout << "Count1: "<< count1 <<endl;
	cout << "Count2: "<< count2 <<endl;
	cout << "count3: "<< count3 <<endl;
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
}
