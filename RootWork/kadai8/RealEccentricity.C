#include <TCanvas.h>
#include <TGraph.h>
#include <TH2D.h>
#include <TGraph2D.h>

Double_t Rho(Double_t r){
	Double_t A_Au = 197.1;
	Double_t R = 1.18*pow(A_Au,1/3.)-0.48;
	Double_t a = 0.54;
	Double_t rho;
	rho = 1 / (1 + exp((r - R)/a));
	return rho;
}

void RealEccentricity(){
	
	TRandom3 *trandom = new TRandom3;
	Int_t A_Au = 197.;
	Double_t imp = 0;
	Double_t phi = 0;
	Double_t posA[3][A_Au];
	Double_t posB[3][A_Au];
	Double_t posA_r[3][197];
	Double_t posB_r[3][197];
	Double_t posA_p[3][197]={0,};
	Double_t posB_p[3][197]={0,};
	Double_t posA_pr[3][197]={0,};
	Double_t posB_pr[3][197]={0,};

	Int_t countA = 0;
	Int_t countB = 0;
	Int_t count_b = 0;
	Int_t count_p = 0;
	Int_t count_pr = 0;
	Int_t count_pA = 0;
	Int_t count_pB = 0;
	Int_t count_pAr = 0;
	Int_t count_pBr = 0;
	Int_t count_s = 0;

	Double_t sum_Ax =0;
	Double_t sum_Bx =0;
	Double_t sum_Ay =0;
	Double_t sum_By =0;
	Double_t avex_2 = 0;
	Double_t avey_2 = 0;
	Double_t ecc = 0;
	
	TCanvas *c1 = new TCanvas("c1", "c1", 1000, 800);
	c1 -> Divide(2, 1);
	TCanvas *c2 = new TCanvas("c2", "c2", 1200, 600);
	c2 -> Divide(2, 1);
	TCanvas *c3 = new TCanvas("c3", "c3", 1200, 600);
	c3 -> Divide(1, 1);
	
	TH2D *hposA = new TH2D("hposA", "hposA;x;y", 100, -20, 20, 100, -20, 20);
	hposA -> SetMarkerStyle(kFullCircle);
	hposA -> SetMarkerColor(kBlue);
	hposA -> SetTitle("Imp = 0,x-y surface");

	TH2D *hposB = new TH2D("hposB", "hposB;x;y", 100, -20, 20, 100, -20, 20);
	hposB -> SetMarkerStyle(kFullSquare);
	hposB -> SetMarkerColor(kRed);
	
	TH2D *hposAB = new TH2D("hposAB", "hposAB;x;y", 100, -20, 20, 100, -20, 20);
	hposAB -> SetMarkerStyle(kFullSquare);
	hposAB -> SetMarkerColor(kOrange);
	hposAB -> SetTitle("Imp = 0,x-y surface");

	TH2D *hposA2 = new TH2D("hposA2", "hposA;x;y", 100, -20, 20, 100, -20, 20);
	hposA2 -> SetMarkerStyle(kFullCircle);
	hposA2 -> SetMarkerColor(kBlue);
	hposA2 -> SetTitle("Imp = R, x-y surface");

	TH2D *hposB2 = new TH2D("hposB2", "hposB;x;y", 100, -20, 20, 100, -20, 20);
	hposB2 -> SetMarkerStyle(kFullSquare);
	hposB2 -> SetMarkerColor(kRed);

	TH2D *hposAB2 = new TH2D("hposAB2", "hposAB2;x;y", 100, -20, 20, 100, -20, 20);
	hposAB2 -> SetMarkerStyle(kFullSquare);
	hposAB2 -> SetMarkerColor(kOrange);
	hposAB2 -> SetTitle("Imp = R, x-y surface");

	TH1D *hecc = new TH1D("","",100,-20,20);
	TH2D *hres = new TH2D("","",100, -2,2,100,0,2*pi);
	hres -> SetTitle("atan(a),phi");
	Double_t X[1000],X3[1000], Y[1000], Y2[1000], Y3[1000], Y4[1000];
	
	for(Int_t n=1; n<130; n++){
		phi = 2*pi*trandom->Rndm();
		countA = 0;
		countB = 0;
		Double_t participant_A[197] = {0,};
		Double_t participant_B[197] = {0,};
		Double_t Distance_A[197] = {0,};
		Double_t Distance_B[197] = {0.};
		Double_t Distance[197][197] = {0,};

		while(countA < A_Au){// Generate Au A
			Double_t Acent[3] = {(Double_t)0+(Double_t)imp/2*cos(phi), imp/2*sin(phi), 0}; //Position of center of Au A
			Double_t Temp[3] = {trandom->Rndm()*14-7, trandom->Rndm()*14-7, trandom->Rndm()*14-7};
			for(Int_t j=0; j<3; j++){
				Distance_A[countA] += pow(Acent[j] - Temp[j], 2);
			}
			Distance_A[countA] = sqrt(Distance_A[countA]);
			Double_t r_A = 0;
			for(Int_t j=0; j<3; j++){
				r_A += pow(Temp[j],2);
			}
			r_A = sqrt(r_A);
			if(Rho(r_A) > (trandom->Rndm())){
				for(Int_t k=0; k<3; k++){
					posA[k][countA] = Acent[k]+Temp[k];
				}
					countA++;
			}
		}
		while(countB < A_Au){// Generate Au B
			Double_t Bcent[3] = {(Double_t)0-(Double_t)imp/2*cos(phi), -imp/2*sin(phi), 0}; //Position of center of Au B
			Double_t Temp[3] = {trandom->Rndm()*14-7, trandom->Rndm()*14-7, trandom->Rndm()*14-7};
			for(Int_t j=0; j<3; j++){
				Distance_B[countB] += pow(Bcent[j] -Temp[j], 2);
			}
			Distance_B[countB] = sqrt(Distance_B[countB]);
			Double_t r_B = 0;
			for(Int_t j=0; j<3; j++){
				r_B += pow(Temp[j],2);
			}
			r_B = sqrt(r_B);
			if(Rho(r_B) > (trandom->Rndm())){
				for(Int_t k=0; k<3; k++){
					posB[k][countB] = Bcent[k]+Temp[k];
				}
					countB++;
			}
		}
		for(countA =0; countA < A_Au; countA++){
		for(countB =0; countB < A_Au; countB++){
			for(Int_t k =0; k <3; k++){
			Distance[countA][countB] += pow(posA[k][countA]-posB[k][countB],2);
			}
			Distance[countA][countB] = sqrt(Distance[countA][countB]);
			Double_t b_max = sqrt(4.2/pi);
			if (Distance[countA][countB] < b_max){
			count_b ++;
			participant_A[countA]=1;
			participant_B[countB]=1;
			}
		}
		}

		Double_t posAB[3][394] = {0,};
		Double_t posAB2[3][394] = {0,};
		for(countA =0; countA<A_Au; countA++){
			count_pA += participant_A[countA];
			for(Int_t k =0; k<3; k++){
			posA_p[k][countA] = participant_A[countA] * posA[k][countA];	
			}	
			sum_Ax += pow(posA_p[0][countA],2);
			sum_Ay += pow(posA_p[1][countA],2);
		
			if (participant_A[countA] == 1) {
				if(imp == 0) {
						posAB[0][countA]=posA[0][countA];
						posAB[1][countA]=posA[1][countA];
						hposA -> Fill(posA[0][countA], posA[1][countA]);
						hposAB -> Fill(posAB[0][countA], posAB[1][countA]);
				}
				if(imp == 6.4){ 
					hposA2 -> Fill(posA[0][countA], posA[1][countA]);
					posAB2[0][countA]=posA[0][countA];
					posAB2[1][countA]=posA[1][countA];
					hposAB2 -> Fill(posAB2[0][countA], posAB2[1][countA]);
				}
			}
		}
		for(countB =0; countB<A_Au; countB++){
			count_pB += participant_B[countB];
			for(Int_t k =0; k<3; k++){
			posB_p[k][countB] = participant_B[countB] * posB[k][countB];	
			}
			sum_Bx += pow(posB_p[0][countB],2);
			sum_By += pow(posB_p[1][countB],2);
			
			if (participant_B[countB] == 1){ 
				if(imp == 0) {
						posAB[0][countB+197]=posB[0][countB];
						posAB[1][countB+197]=posB[1][countB];
						hposB -> Fill(posB[0][countB], posB[1][countB]);
						hposAB -> Fill(posAB[0][countB+197], posAB[1][countB+197]);
				}
				if(imp == 6.4) {
					hposB2 -> Fill(posB[0][countB], posB[1][countB]);
					posAB2[0][countB+197]=posB[0][countB];
					posAB2[1][countB+197]=posB[1][countB];
					hposAB2 -> Fill(posAB2[0][countB+197], posAB2[1][countB+197]);
				}
			}
		}
		Double_t SX,SY,SXY,SXP = 0;
		Double_t a = 0;
		for(Int_t l = 0; l<394; l++){
			SX +=posAB[0][l];
			SY +=posAB[1][l];
			SXY +=posAB[0][l]*posAB[1][l];
			SXP += pow(posAB[0][l],2);
		}
		a = (394*SXY-SX*SY)/(394*SXP-pow(SX,2));

		hres -> Fill(atan(a),phi);
		
		count_p = count_pA + count_pB;
		count_s = 2*A_Au-count_p;
	
		avex_2 = (sum_Ax+sum_Bx)/count_p;
		avey_2 = (sum_Ay+sum_By)/count_p;

		ecc = abs((avex_2-avey_2)/(avex_2+avey_2));
		
		X[n] = n*0.1-0.1;
		Y[n] = count_b;
		Y2[n] = count_p;
		Y3[n] = count_s;
		Y4[n] = ecc;
		imp = 0.1*n;
		count_b =0;
		count_p =0;
		count_pA = 0;
		count_pB = 0;
		count_s =0;
		ecc = 0;
	}

	TF1 *f1 = new TF1("f1","[0]*x+[1]",-10,10);
	TF1 *f2 = new TF1("f2","[0]*x+[1]",-10,10);
	c1 -> cd(1);
	hposA -> Draw("");
	hposB -> Draw("SAME");
	
	c1 -> cd(2);
	hposA2 -> Draw("");
	hposB2 -> Draw("SAME");

	c2 -> cd(1);
	hposAB ->Draw("");
	hposAB -> Fit(f1);
	
	c2 -> cd(2);
	hposAB2 ->Draw("");
	hposAB2 ->Fit(f2);
	
	c3 -> cd(1);
	hres -> Draw("colz");
	return;
}
