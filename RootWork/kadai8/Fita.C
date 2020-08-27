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
// root fitting confirm 


void Fita(){
	
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
	
	TCanvas *c1 = new TCanvas("c1", "c1", 1000, 1000);
	c1 -> Divide(2, 2);
	TCanvas *c2 = new TCanvas("c2", "c2", 1000, 1000);
	c2 -> Divide(2, 2);
	TCanvas *c3 = new TCanvas("c3", "c3", 1000, 1000);
	c3 -> Divide(1, 1);
	
	TCanvas *c4 = new TCanvas("c4","c4", 1000, 1000);
	c4 -> Divide(2, 2);

	TH2D *hposAB = new TH2D("hposAB", "hposAB;x;y", 100, -20, 20, 100, -20, 20);
	hposAB -> SetMarkerStyle(kFullSquare);
	hposAB -> SetMarkerColor(kOrange);
	hposAB -> SetTitle("Imp = 0,x-y surface");

	TH2D *hposAB2 = new TH2D("hposAB2", "hposAB2;x;y", 100, -20, 20, 100, -20, 20);
	hposAB2 -> SetMarkerStyle(kFullSquare);
	hposAB2 -> SetMarkerColor(kOrange);
	hposAB2 -> SetTitle("Imp = R, x-y surface");

	TH1D *hres1 = new TH1D("hres1","",100,-2,2);
	hres1 -> SetTitle("atan(a1)");

	TH1D *hres2 = new TH1D("hres2","",100,-2,2);
	hres2 -> SetTitle("atan(Fitting(a1))");

	TH1D *hres3 = new TH1D("hres1","",100,-2,2);
	hres3 -> SetTitle("atan(1/a1)");
	
	TH1D *hres4 = new TH1D("hres1","",100,-2,2);
	hres4 -> SetTitle("atan(Fitting(a2))");
	
	TH2D *hres11 = new TH2D("","",100, -2, 2,100,-pi,pi);
	hres11 -> SetTitle("atan(a1),phi");
	
	TH2D *hres22 = new TH2D("","",100, -2, 2,100,-pi,pi);
	hres22 -> SetTitle("atan(fitting(a1)),phi");
	
	TH2D *hres33 = new TH2D("","",100, -2, 2,100,-pi,pi);
	hres33 -> SetTitle("atan(1/a1),phi");

	TH2D *hres44 = new TH2D("","",100,-2,2,100,-pi,pi);
	hres44 -> SetTitle("atan(Fitting(a2)),phi");

	TH2D *hres444 = new TH2D("","",100,-2,2,100,-pi,pi);
	hres444 -> SetTitle("atan(1/Fitting(a2))rev,phi");
	
	TH2D *hres4444 = new TH2D("","",100,-2,2,100,-pi,pi);
	hres4444 -> SetTitle("atan(1/Fitting(a2)),phi");
		
	TF1 *f1 = new TF1("f1","[0]*x+[1]",-10,10);
	TF1 *f2 = new TF1("f2","[0]*x+[1]",-10,10);
	
	for(Int_t n=1; n<700; n++){
	
		hposAB->Reset();
		hposAB2->Reset();
		
		phi = pi*(trandom->Rndm()*2.0-1.0);
		imp = 10;//14*trandom->Rndm();
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
		Double_t posABR[3][394][700] = {0,};
		
		for(countA =0; countA<A_Au; countA++){
			count_pA += participant_A[countA];
			for(Int_t k =0; k<3; k++){
			posA_p[k][countA] = participant_A[countA] * posA[k][countA];	
			}	
			sum_Ax += pow(posA_p[0][countA],2);
			sum_Ay += pow(posA_p[1][countA],2);
		
			if (participant_A[countA] == 1) {
						posAB[0][countA]=posA[0][countA];
						posAB[1][countA]=posA[1][countA];
						hposAB -> Fill(posAB[0][countA], posAB[1][countA]);
						hposAB2 -> Fill(posAB[1][countA], posAB[0][countA]);
			
				posABR[0][countA][n] = posA[0][countA];
				posABR[1][countA][n] = posA[1][countA];
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
						posAB[0][countB+197]=posB[0][countB];
						posAB[1][countB+197]=posB[1][countB];
						hposAB -> Fill(posAB[0][countB+197], posAB[1][countB+197]);
						hposAB2 -> Fill(posAB[1][countB+197], posAB[0][countB+197]);
				
				posABR[0][countB+197][n] = posB[0][countB];
				posABR[1][countB+197][n] = posB[1][countB];
			}
		}
		Double_t SX[700],SY[700],SXY[700],SXP[700] = {0,};
		count_p = count_pA + count_pB;
		Double_t a = 0;
		for(Int_t l = 0; l<394; l++){
			SX[n] +=posABR[0][l][n];
			SY[n] +=posABR[1][l][n];
			SXY[n] +=posABR[0][l][n]*posABR[1][l][n];
			SXP[n] += pow(posABR[0][l][n],2);
		}
		a = (count_p*SXY[n]-SX[n]*SY[n])/(count_p*SXP[n]-pow(SX[n],2));
		hres1->Fill(atan(a));

		hres11->Fill(atan(a),phi);
		
		if(atan(1/a)>0){
		
		hres3->Fill(atan(1/a)-pi/2);
		hres33->Fill(atan(1/a)-pi/2,phi);
		
		}

		if(atan(1/a<0)){
		
		hres3->Fill(atan(1/a)+pi/2);
		hres33->Fill(atan(1/a)+pi/2,phi);

		}
		count_s = 2*A_Au-count_p;
	
		
		count_b =0;
		count_p =0;
		count_pA = 0;
		count_pB = 0;
		count_s =0;
		ecc = 0;
		
		hposAB->Draw("");
		hposAB->Fit("f1","q");
		
		hres2->Fill(atan(f1->GetParameter(0)));
		hres22->Fill(atan(f1->GetParameter(0)),phi);

		hposAB2->Draw("");
		hposAB2->Fit("f2","q");
		hres4->Fill(atan(f2->GetParameter(0)));
		hres44->Fill(atan(f2->GetParameter(0)),phi);
	
		if(atan(1/(f2->GetParameter(0)))>0){
		
			hres444->Fill(atan(1/(f2->GetParameter(0)))-pi/2,phi);
		}	
		if(atan(1/(f2->GetParameter(0)))<0){
		
			hres444->Fill(atan(1/(f2->GetParameter(0)))+pi/2,phi);
		}

			hres4444->Fill(atan(1/(f2->GetParameter(0))),phi);
	
	}

	c1 -> cd(1);
	hres1 -> Draw("colz");
	
	c1 -> cd(2);
	hres2 -> Draw("colz");
	c1 -> cd(3);
	hres3 -> Draw("colz");
	c1 -> cd(4);
	hres4 -> Draw("colz");
	c2 -> cd(1);
	hres11 -> Draw("colz");
	c2 -> cd(2);
	hres22 -> Draw("colz");
	c2 -> cd(3);
	hres33 -> Draw("colz");
	c2 -> cd(4);
	hres44 -> Draw("colz");
	
	hres22 -> SetMarkerStyle(kFullSquare);
	hres22 -> SetMarkerColor(kRed);

	hres4444 -> SetMarkerStyle(kFullSquare);
	hres4444 -> SetMarkerColor(kBlue);


	c3 -> cd(1);
	hres22->SetTitle("atan(Fitting(a1),phi).rev");
	hres22 -> Draw("");
	hres4444-> Draw("same");

/*	c3-> cd(2);
	hres444->Draw("colz");

	c3 ->cd(3);
	hres22 ->Draw("colz");

	c3 ->cd(4);
	hres4444->Draw("colz");
*/	
	return;
}
