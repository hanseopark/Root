#include <TCanvas.h>

Double_t Rho(Double_t r){
	IntT A_Au = 197;
	Double_t R = 1.18*pow(A_Au,1/3.)-0.48;
	Double_t a = 0.54;
	Double_t rho;
	rho = 1 / (1+ exp((r-R)/a));
	return rho;
}

void Fita2(){

	TCanvas *c1 = new TCanvas("c1","c1",1000,1000);
	c1 -> Divide(2,2);

	TRandom3 *trandom = new TRandom3;
	Int_t A_Au = 197;
	Int_t countA = 0;
	Int_t countB = 0;

	for(Int_t n=0; n<1000; n++){
	Phi = pi*(trandom->Rndm()*2.0-1.0);
	imp = 14* trandom ->Rndm();
		
	while(countA < A_Au){

		Double_t r = 14*trandom->Rndm();
		Double_t phi = pi*(trandom->Rndm()*2.0-1.0);
		Doulbe_t Acent[2] = {r, Phi};
		Double_t Temp[2] = {r,phi}

	}

	
	}


}
