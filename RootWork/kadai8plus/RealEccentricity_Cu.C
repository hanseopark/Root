#include <TCanvas.h>
#include <TGraph.h>
#include <TH2D.h>
#include <TGraph2D.h>

Double_t Rho(Double_t r){
	Double_t A_Cu = 63.5;
	Double_t R = 1.18*pow(A_Cu,1/3.)-0.48;
	Double_t a = 0.54;
	Double_t rho;
	rho = 1 / (1 + exp((r - R)/a));
	return rho;
}

void RealEccentricity_Cu(){
	
	TRandom3 *trandom = new TRandom3;
	Int_t A_Cu = 64.;
	Double_t imp = 0;
	Double_t phi = 0;
	Double_t posA[3][A_Cu];
	Double_t posB[3][A_Cu];
	Double_t posA_r[3][64];
	Double_t posB_r[3][64];
	Double_t posA_p[3][64]={0,};
	Double_t posB_p[3][64]={0,};
	Double_t posA_pr[3][64]={0,};
	Double_t posB_pr[3][64]={0,};

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
	
	TH2D *hposA = new TH2D("hposA", "hposA;x;y", 100, -20, 20, 100, -20, 20);
	hposA -> SetMarkerStyle(kFullCircle);
	hposA -> SetMarkerColor(kBlue);
	hposA -> SetTitle("Imp = 0,x-y surface");

	TH2D *hposB = new TH2D("hposB", "hposB;x;y", 100, -20, 20, 100, -20, 20);
	hposB -> SetMarkerStyle(kFullSquare);
	hposB -> SetMarkerColor(kRed);

	TH2D *hposA2 = new TH2D("hposA2", "hposA;x;y", 100, -20, 20, 100, -20, 20);
	hposA2 -> SetMarkerStyle(kFullCircle);
	hposA2 -> SetMarkerColor(kBlue);
	hposA2 -> SetTitle("Imp = R, x-y surface");

	TH2D *hposB2 = new TH2D("hposB2", "hposB;x;y", 100, -20, 20, 100, -20, 20);
	hposB2 -> SetMarkerStyle(kFullSquare);
	hposB2 -> SetMarkerColor(kRed);

	TH1D *hecc = new TH1D("","",100,-20,20);
	TH1D *hres = new TH1D("","",100, 0, 100);
	hres -> SetTitle("imp=4.2,phi=random");
	Double_t X[1000],X3[1000], Y[1000], Y2[1000], Y3[1000], Y4[1000];
	for(Int_t n=1; n<130; n++){
		phi = 2*pi*trandom->Rndm();
		countA = 0;
		countB = 0;
		Double_t participant_A[64] = {0,};
		Double_t participant_B[64] = {0,};
		Double_t Distance_A[64] = {0,};
		Double_t Distance_B[64] = {0.};
		Double_t Distance[64][64] = {0,};

		while(countA < A_Cu){// Generate Au A
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
		while(countB < A_Cu){// Generate Au B
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
		for(countA =0; countA < A_Cu; countA++){
		for(countB =0; countB < A_Cu; countB++){
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
		for(countA =0; countA<A_Cu; countA++){
			count_pA += participant_A[countA];
			for(Int_t k =0; k<3; k++){
			posA_p[k][countA] = participant_A[countA] * posA[k][countA];	
			}	
			sum_Ax += pow(posA_p[0][countA],2);
			sum_Ay += pow(posA_p[1][countA],2);
		
			if (participant_A[countA] == 1) {
				if(imp == 0) hposA -> Fill(posA[0][countA], posA[1][countA]);
				if(imp == 4.2) hposA2 -> Fill(posA[0][countA], posA[1][countA]);
			}
		}
		for(countB =0; countB<A_Cu; countB++){
			count_pB += participant_B[countB];
			for(Int_t k =0; k<3; k++){
			posB_p[k][countB] = participant_B[countB] * posB[k][countB];	
			}
			sum_Bx += pow(posB_p[0][countB],2);
			sum_By += pow(posB_p[1][countB],2);
			
			if (participant_B[countB] == 1){ 
				if(imp == 0) hposB -> Fill(posB[0][countB], posB[1][countB]);
				if(imp == 4.2) hposB2 -> Fill(posB[0][countB], posB[1][countB]);
			}
		}
			
		count_p = count_pA + count_pB;
		count_s = 2*A_Cu-count_p;
	
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

	Double_t resolution_phiave = 0;
	//imp = 4.2
	for(Int_t t=0; t<100; t++){
		imp = 4.2;
		phi = 2*pi*trandom->Rndm();
		countA = 0;
		countB = 0;
		Double_t participant_A[64] = {0,};
		Double_t participant_B[64] = {0,};
		Double_t participant_Ar[64] = {0,};
		Double_t participant_Br[64] = {0,};
		Double_t Distance_A[64] = {0,};
		Double_t Distance_B[64] = {0.};
		Double_t Distance_Ar[64] = {0,};
		Double_t Distance_Br[64] = {0,};
		Double_t Distance[64][64] = {0,};
		Double_t Distance_r[64][64] = {0,};

		while(countA < A_Cu){// Generate Au A (x-z surface)
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
		countA = 0;
		while(countA < A_Cu){// Generate Au A (rot phi surface)
			Double_t Acent[3] = {(Double_t)0+(Double_t)imp/2*cos(phi), imp/2*sin(phi), 0}; //Position of center of Au A
			Double_t Temp[3] = {trandom->Rndm()*14-7, trandom->Rndm()*14-7, trandom->Rndm()*14-7};
			for(Int_t j=0; j<3; j++){
				Distance_Ar[countA] += pow(Acent[j] - Temp[j], 2);
			}
			Distance_Ar[countA] = sqrt(Distance_A[countA]);
			Double_t r_A = 0;
			for(Int_t j=0; j<3; j++){
				r_A += pow(Temp[j],2);
			}
			r_A = sqrt(r_A);
			if(Rho(r_A) > (trandom->Rndm())){
				for(Int_t k=0; k<3; k++){
					posA_r[k][countA] = Acent[k]+Temp[k];
				}
					countA++;
			}
		}
		while(countB < A_Cu){// Generate Au B (x-z surface)
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
		countB=0;
		while(countB < A_Cu){// Generate Au B (rot phi surface)
			Double_t Bcent[3] = {(Double_t)0-(Double_t)imp/2*cos(phi), -imp/2*sin(phi), 0}; //Position of center of Au B
			Double_t Temp[3] = {trandom->Rndm()*14-7, trandom->Rndm()*14-7, trandom->Rndm()*14-7};
			for(Int_t j=0; j<3; j++){
				Distance_Br[countB] += pow(Bcent[j] -Temp[j], 2);
			}
			Distance_Br[countB] = sqrt(Distance_B[countB]);
			Double_t r_B = 0;
			for(Int_t j=0; j<3; j++){
				r_B += pow(Temp[j],2);
			}
			r_B = sqrt(r_B);
			if(Rho(r_B) > (trandom->Rndm())){
				for(Int_t k=0; k<3; k++){
					posB_r[k][countB] = Bcent[k]+Temp[k];
				}
				countB++;
			}

		}
		for(countA =0; countA < A_Cu; countA++){
		for(countB =0; countB < A_Cu; countB++){
			Distance[countA][countB] = pow(posA[0][countA]-posB[0][countB],2)+pow(posA[2][countA]-posB[2][countB],2);

			Distance[countA][countB] = sqrt(Distance[countA][countB]);
			Double_t b_max = sqrt(4.2/pi);
			if (Distance[countA][countB] < b_max){
			participant_A[countA]=1;
			participant_B[countB]=1;
			}
			for(Int_t k=0; k<3; k++){
			Distance_r[countA][countB] += pow(posA_r[k][countA]-posB_r[k][countB],2);
			}
			Distance_r[countA][countB] = sqrt(Distance_r[countA][countB]);
			if (Distance_r[countA][countB] < b_max){
			participant_Ar[countA]=1;
			participant_Br[countB]=1;
			}	
		}
		}
	
		Double_t sum_pxr = 0;
		Double_t sum_pyr = 0;
		Double_t sum_px = 0;
		Double_t sum_py = 0;
		Double_t resolution_x = 0;
		Double_t resolution_y = 0;
		Double_t resolution_phi = 0;

		for(countA =0; countA<A_Cu; countA++){
			count_pA += participant_A[countA];
			count_pAr += participant_Ar[countA];
			for(Int_t k =0; k<3; k++){
			posA_p[k][countA] = participant_A[countA] * posA[k][countA];	
			}
			for(Int_t k =0; k<3; k++){
			posA_pr[k][countA] = participant_Ar[countA]*posA_r[k][countA];
			}
			sum_Ax += pow(posA_p[0][countA],2);
			sum_Ay += pow(posA_p[1][countA],2);
		
			if (participant_A[countA] == 1) {
			}
		}
		for(countB =0; countB<A_Cu; countB++){
			count_pB += participant_B[countB];
			count_pBr += participant_Br[countB];
			for(Int_t k =0; k<3; k++){
			posB_p[k][countB] = participant_B[countB] * posB[k][countB];	
			}
			for(Int_t k =0; k<3; k++){
			posB_pr[k][countB] = participant_Br[countB] * posB_r[k][countB];	
			}
			sum_Bx += pow(posB_p[0][countB],2);
			sum_By += pow(posB_p[1][countB],2);
	
			if (participant_B[countB] == 1){ 
			}
		}

		for(Int_t k = 0; k<A_Cu; k++){
			sum_px += posA_p[0][k]+posB_p[0][k];
			sum_py += posA_p[1][k]+posB_p[1][k];
			sum_pxr += posA_pr[0][k]+posB_pr[0][k];
			sum_pyr += posA_pr[1][k]+posB_pr[1][k];
		}

		count_p = count_pA + count_pB;
		count_pr = count_pAr + count_pBr;
		count_s = 2*A_Cu-count_p;
	
		resolution_x = abs(sum_pxr/count_pr-sum_px/count_p);
		resolution_y = abs(sum_pyr/count_pr-sum_py/count_p);
		resolution_phi = atan(resolution_x/resolution_y);
		resolution_phiave += resolution_phi;

		avex_2 = (sum_Ax+sum_Bx)/count_p;
		avey_2 = (sum_Ay+sum_By)/count_p;

		ecc = abs((avex_2-avey_2)/(avex_2+avey_2));
	
		hecc ->Fill(imp, ecc);
		hres ->Fill(t,resolution_phi);	
		count_b =0;
		count_p =0;
		count_pr = 0;
		count_pA = 0;
		count_pB = 0;
		count_pAr =0;
		count_pBr =0;
		count_s =0;
		ecc = 0;
	}

	resolution_phiave = resolution_phiave/100;
	cout << "phiの分解能は: "<< resolution_phiave <<endl;


/*	TGraph *gbinary = new TGraph(130,X,Y);
	gbinary -> SetMarkerStyle(6);
	gbinary -> SetLineColor(3);
	gbinary -> SetMarkerSize(3);
	gbinary -> SetMarkerColor(kRed);
	gbinary -> SetTitle("Binary");
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
*/
//	c1 -> cd(1);
//	gbinary -> Draw("ACP");
//	gparticipant ->Draw("same");
//	gspectator -> Draw("same");

	c1 -> cd(1);
	hposA -> Draw("");
	hposB -> Draw("SAME");

	c1 -> cd(2);
	hposA2 -> Draw("");
	hposB2 -> Draw("SAME");

	c2 -> cd(1);
	hecc ->Draw("");

	c2 -> cd(2);
	hres ->Draw("p");
	
	return;
}
