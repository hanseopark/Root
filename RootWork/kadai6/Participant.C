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

void Participant(){
	
	TRandom3 *trandom = new TRandom3;
	Int_t A_Au = 197.;
	Double_t imp = 0;
	Double_t posA[3][A_Au];
	Double_t posB[3][A_Au];
	Double_t posA_p[3][197]={0,};
	Double_t posB_p[3][197]={0,};

	Int_t countA = 0;
	Int_t countB = 0;
	Int_t count_b = 0;
	Int_t count_p = 0;
	Int_t count_pA = 0;
	Int_t count_pB = 0;
	Int_t count_s = 0;

	Double_t sum_Ax =0;
	Double_t sum_Bx =0;
	Double_t sum_Ay =0;
	Double_t sum_By =0;
	Double_t avex_2 = 0;
	Double_t avey_2 = 0;
	Double_t ecc = 0;

	TCanvas *c1 = new TCanvas("c1", "c1", 1500, 1000);
	c1 -> Divide(1, 1);
	TCanvas *c2 = new TCanvas("c2", "c2", 1200, 600);
	c2 -> Divide(2, 1);
	TCanvas *c3 = new TCanvas("c3", "c3", 1000, 800);
	c3 -> Divide(1, 1);

	TH3D *hposA = new TH3D("hposA", "hposA;x;y", 100, -20, 20, 100, -20, 20, 100, -20, 20);
	hposA -> SetMarkerStyle(kFullCircle);
	hposA -> SetMarkerColor(kBlue);
	hposA -> SetTitle("Imp = 0");

	TH3D *hposB = new TH3D("hposB", "hposB;x;y", 100, -20, 20, 100, -20, 20, 100, -20, 20);
	hposB -> SetMarkerStyle(kFullSquare);
	hposB -> SetMarkerColor(kRed);

	TH3D *hposA2 = new TH3D("hposA2", "hposA;x;y", 100, -20, 20, 100, -20, 20, 100, -20, 20);
	hposA2 -> SetMarkerStyle(kFullCircle);
	hposA2 -> SetMarkerColor(kBlue);
	hposA2 -> SetTitle("Imp = R");

	TH3D *hposB2 = new TH3D("hposB2", "hposB;x;y", 100, -20, 20, 100, -20, 20, 100, -20, 20);
	hposB2 -> SetMarkerStyle(kFullSquare);
	hposB2 -> SetMarkerColor(kRed);

	TH3D *hposA3 = new TH3D("hposA3", "hposA;x;y", 100, -20, 20, 100, -20, 20, 100, -20, 20);
	hposA3 -> SetMarkerStyle(kFullCircle);
	hposA3 -> SetMarkerColor(kBlue);
	hposA3 -> SetTitle("Imp = 2R");

	TH3D *hposB3 = new TH3D("hposB3", "hposB;x;y", 100, -20, 20, 100, -20, 20, 100, -20, 20);
	hposB3 -> SetMarkerStyle(kFullSquare);
	hposB3 -> SetMarkerColor(kRed);

	Double_t X[1000],X3[1000], Y[1000], Y2[1000], Y3[1000], Y4[1000];
	for(Int_t n=1; n<130; n++){
		countA = 0;
		countB = 0;
		Double_t participant_A[197] = {0,};
		Double_t participant_B[197] = {0,};
		Double_t Distance_A[197] = {0,};
		Double_t Distance_B[197] = {0.};
		Double_t Distance[197][197] = {0,};
		sum_Ax=sum_Bx=sum_Ay=sum_By=0;

		while(countA < A_Au){// Generate Au A
			Double_t Acent[3] = {(Double_t)0+(Double_t)imp/2, 0, 0}; //Position of center of Au A
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
			Double_t Bcent[3] = {(Double_t)0-(Double_t)imp/2, 0, 0}; //Position of center of Au B
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
			Distance[countA][countB] = pow(posA[0][countA]-posB[0][countB],2)+pow(posA[2][countA]-posB[2][countB],2);

			Distance[countA][countB] = sqrt(Distance[countA][countB]);
			Double_t b_max = sqrt(4.2/pi);
			if (Distance[countA][countB] < b_max){
			count_b ++;
			participant_A[countA]=1;
			participant_B[countB]=1;
			}
		}
		}
		
		for(countA =0; countA<A_Au; countA++){
			count_pA += participant_A[countA];
			for(Int_t k =0; k<3; k++){
			posA_p[k][countA] = participant_A[countA] * posA[k][countA];	
			}	
			sum_Ax += pow(posA_p[0][countA],2);
			sum_Ay += pow(posA_p[1][countA],2);
		
			if (participant_A[countA] == 1) {
				if(imp == 0) hposA -> Fill(posA[0][countA], posA[1][countA],posA[2][countA]);
				if(imp == 6.4) hposA2 -> Fill(posA[0][countA], posA[1][countA],posA[2][countA]);
				if(imp == 12.8) hposA3 -> Fill(posA[0][countA], posA[1][countA],posA[2][countA]);
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
				if(imp == 0) hposB -> Fill(posB[0][countB], posB[1][countB],posB[2][countB]);
				if(imp == 6.4) hposB2 -> Fill(posB[0][countB], posB[1][countB],posB[2][countB]);
				if(imp == 12.8) hposB3 -> Fill(posB[0][countB], posB[1][countB],posB[2][countB]);
			}
		}
			
		count_p = count_pA + count_pB;
		count_s = 2*A_Au-count_p;
	
		avex_2 = (sum_Ax+sum_Bx)/count_p;
		avey_2 = (sum_Ay+sum_By)/count_p;

		ecc = (avex_2-avey_2)/(avex_2+avey_2);
		
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

	TGraph *gbinary = new TGraph(130,X,Y);
	gbinary -> SetMarkerStyle(6);
	gbinary -> SetLineColor(3);
	gbinary -> SetMarkerSize(3);
	gbinary -> SetMarkerColor(kRed);
	gbinary -> SetTitle("Binary Collision,Participant,Spectator");
	gbinary -> GetXaxis() -> SetTitle("Impact Parameter");
	gbinary -> GetYaxis() -> SetTitle("collision");

	TGraph *gparticipant = new TGraph(130,X,Y2);
	gparticipant -> SetMarkerStyle(6);
	gparticipant -> SetLineColor(6);
	gparticipant -> SetMarkerSize(3);
	gparticipant -> SetMarkerColor(kBlue);
	gparticipant -> SetTitle("Participant");

	TGraph *gspectator = new TGraph(130,X,Y3);
	gspectator -> SetMarkerStyle(6);
	gspectator -> SetLineColor(9);
	gspectator -> SetMarkerSize(3);
	gspectator -> SetMarkerColor(kOrange);

	TGraph *geccentricity = new TGraph(130,X,Y4);
	geccentricity -> SetMarkerStyle(6);
	geccentricity -> SetLineColor(9);
	geccentricity -> SetMarkerSize(3);
	geccentricity -> SetMarkerColor(kOrange);
	geccentricity -> SetTitle("Eccentricity");
	geccentricity -> GetXaxis() -> SetTitle("Impact Parameter");
	geccentricity -> GetYaxis() -> SetTitle("Eccentricity");

	c1 -> cd(1);
	gbinary -> Draw("ACP");
	gparticipant ->Draw("same");
	gspectator -> Draw("same");

	c2 -> cd(1);
	hposA ->Project3D("yx") -> Draw("");
	hposB ->Project3D("yx") -> Draw("same");
	
	c2 -> cd(2);
	hposA2 -> Project3D("yx") -> Draw("");
	hposB2 -> Project3D("yx") -> Draw("SAME");
/*
	c2 -> cd(3);
	hposA3 -> Draw("");
	hposB3 -> Draw("SAME");
*/
	c3 -> cd(1);
	geccentricity ->Draw("ACP");

	return;
}
