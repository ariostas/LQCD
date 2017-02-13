#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <iostream>
#include <fstream>
#include <TString.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TComplex.h>
#include <TRandom3.h>

#include "GaugeField.h"

#endif

using namespace TMath;
using namespace std;

// Define common objects
vector<vector<TComplex> > id, nu, nu2;
TRandom3 *randGen = new TRandom3();
TF1 *func = new TF1("x_0 dist", "TMath::Sqrt(1.-x*x)*TMath::Exp(2./3.*[0] * x)", -.99999, .99999);

// Declare functions
void initialize(GaugeField<vector<vector<TComplex> > >*, TString);
void singleUpdate(GaugeField<vector<vector<TComplex> > >*,  UInt_t, UInt_t, UInt_t, UInt_t, UInt_t);
vector<vector<TComplex> > randomMatrix(vector<vector<TComplex> >);
void update(GaugeField<vector<vector<TComplex> > >*, Int_t);
void measure(GaugeField<vector<vector<TComplex> > >*, Int_t);
TGraph* readData(TString);
void graph(TGraph*, TGraph*, const TString, TCanvas*, const TString, const TString, const TString);

vector<vector<TComplex> > multiply(vector<vector<TComplex> >, vector<vector<TComplex> >, Int_t);
vector<vector<TComplex> > multiply(TComplex, vector<vector<TComplex> >, Int_t);
vector<vector<TComplex> > add(vector<vector<TComplex> >, vector<vector<TComplex> >);
vector<vector<TComplex> > dagger(vector<vector<TComplex> >);
TComplex determinant(vector<vector<TComplex> >, Int_t);
TComplex trace(vector<vector<TComplex> >);

// Define parameters
Int_t Ns = 4, Nt = 8;
Int_t Ntherm = 0, Ncorr = 1, Nconfig = 50;
Double_t beta = 5.8;

// Define storage variables
Double_t entryArray[1000];
Double_t plaqArray[1000];

// Main function
void lqcd(){

	cout << "\n\nStarting process...\n\n";

	vector<TComplex> row1; row1.push_back(1); row1.push_back(0); row1.push_back(0);
	vector<TComplex> row2; row2.push_back(0); row2.push_back(1); row2.push_back(0);
	vector<TComplex> row3; row3.push_back(0); row3.push_back(0); row3.push_back(1);
	id.push_back(row1); id.push_back(row2); id.push_back(row3);
	row1[0] = 0; row2[1] = 0; row3[2] = 0;
	nu.push_back(row1); nu.push_back(row2); nu.push_back(row3);
	vector<TComplex> row; row.push_back(0); row.push_back(0);
	nu2.push_back(row); nu2.push_back(row);

	GaugeField<vector<vector<TComplex> > > field(Ns, Ns, Ns, Nt);

	initialize(&field, "cold");
	update(&field, Ntherm);
	for(Int_t nconfig = 0; nconfig < Nconfig; nconfig++){
		measure(&field, nconfig);
		if(nconfig < Nconfig-1) update(&field, Ncorr);
	}

	TGraph *chromaData = readData("chromaData");

	TGraph *myData = new TGraph(Nconfig, entryArray, plaqArray);

	TCanvas *can = new TCanvas("Canvas", "Canvas", 1600, 900);
	graph(chromaData, myData, "Chroma vs script comparison", can, "nth measurement", "Average plaquette (|Re(Tr(U_{#mu#nu}))|)", "comparison");

}

void initialize(GaugeField<vector<vector<TComplex> > > *f, TString mode){

	if(mode == "cold"){
		for(Int_t nx = 0; nx < Ns; nx++){
			for(Int_t ny = 0; ny < Ns; ny++){
				for(Int_t nz = 0; nz < Ns; nz++){
					for(Int_t nt = 0; nt < Nt; nt++){
						for(Int_t dir = 0; dir < 4; dir++){
							(*f)(nx,ny,nz,nt,dir) = id;
						}
					}
				}
			}
		}
	}

	else if(mode == "hot"){

		cout << "Not yet implemented" << endl;

	}

}

void singleUpdate(GaugeField<vector<vector<TComplex> > > *f, UInt_t x, UInt_t y, UInt_t z, UInt_t t, UInt_t d){

	// Sum of staples
	vector<vector<TComplex> > A = nu;

	for(UInt_t dir = 0; dir < 4; dir++){
		if(d == dir) continue;

		vector<vector<TComplex> > tempA1, tempA2;

		UInt_t tempx = x, tempy = y, tempz = z, tempt = t;

		if(d == 0) tempx++;
		else if(d == 1) tempy++;
		else if(d == 2) tempz++;
		else if(d == 3) tempt++;

		tempA1 = (*f)(tempx, tempy, tempz, tempt, dir);

		if(dir == 0) tempx++;
		else if(dir == 1) tempy++;
		else if(dir == 2) tempz++;
		else if(dir == 3) tempt++;
		if(d == 0) tempx--;
		else if(d == 1) tempy--;
		else if(d == 2) tempz--;
		else if(d == 3) tempt--;

		tempA1 = multiply(tempA1, dagger((*f)(tempx, tempy, tempz, tempt, d)), 3);

		if(dir == 0) tempx--;
		else if(dir == 1) tempy--;
		else if(dir == 2) tempz--;
		else if(dir == 3) tempt--;

		tempA1 = multiply(tempA1, dagger((*f)(tempx, tempy, tempz, tempt, dir)), 3);

		if(d == 0) tempx++;
		else if(d == 1) tempy++;
		else if(d == 2) tempz++;
		else if(d == 3) tempt++;
		if(dir == 0) tempx--;
		else if(dir == 1) tempy--;
		else if(dir == 2) tempz--;
		else if(dir == 3) tempt--;

		tempA2 = dagger((*f)(tempx, tempy, tempz, tempt, dir));

		if(d == 0) tempx--;
		else if(d == 1) tempy--;
		else if(d == 2) tempz--;
		else if(d == 3) tempt--;

		tempA2 = multiply(tempA2, dagger((*f)(tempx, tempy, tempz, tempt, d)), 3);
		tempA2 = multiply(tempA2, (*f)(tempx, tempy, tempz, tempt, dir), 3);

		A = add(A,tempA1);
		A = add(A,tempA2);

	}

	vector<vector<TComplex> > tempA;
	vector<TComplex> row1; row1.push_back(TComplex(0)); row1.push_back(TComplex(0));
	tempA.push_back(row1); tempA.push_back(row1);

	// Construct R
	vector<vector<TComplex> > R = id, tempR;

	tempA[0][0] = A[0][0]; tempA[0][1] = A[0][1];
	tempA[1][0] = A[1][0]; tempA[1][1] = A[1][1];

	tempR = randomMatrix(tempA);

	R[0][0] = tempR[0][0]; R[0][1] = tempR[0][1];
	R[1][0] = tempR[1][0]; R[1][1] = tempR[1][1];

	// Construct S
	vector<vector<TComplex> > S = id, tempS;

	vector<vector<TComplex> > AR = multiply(A,R,3);

	tempA[0][0] = AR[0][0]; tempA[0][1] = AR[0][2];
	tempA[1][0] = AR[2][0]; tempA[1][1] = AR[2][2];

	tempS = randomMatrix(tempA);

	S[0][0] = tempS[0][0]; S[0][2] = tempS[0][1];
	S[2][0] = tempS[1][0]; S[2][2] = tempS[1][1];

	// Construct T
	vector<vector<TComplex> > T = id, tempT;

	vector<vector<TComplex> > ARS = multiply(A,S,3);

	tempA[0][0] = ARS[1][1]; tempA[0][1] = ARS[1][2];
	tempA[1][0] = ARS[2][1]; tempA[1][1] = ARS[2][2];

	tempT = randomMatrix(tempA);

	T[1][1] = tempT[0][0]; T[1][2] = tempT[0][1];
	T[2][1] = tempT[1][0]; T[2][2] = tempT[1][1];

	// Construct new link

	(*f)(x,y,z,t,d) = multiply(T,multiply(S,multiply(R,(*f)(x,y,z,t,d),3),3),3);

}

void update(GaugeField<vector<vector<TComplex> > > *f, Int_t NSweeps){

	for(Int_t n = 0; n < NSweeps; n++){

		for(Int_t dir = 0; dir < 4; dir++){
			for(Int_t nx = 0; nx < Ns; nx++){
				for(Int_t ny = 0; ny < Ns; ny++){
					for(Int_t nz = 0; nz < Ns; nz++){
						for(Int_t nt = 0; nt < Nt; nt++){
							singleUpdate(f, nx, ny, nz, nt, dir);
						}
					}
				}
			}
		}

	}

}

vector<vector<TComplex> > randomMatrix(vector<vector<TComplex> > A){

	if(A.size() != 2) cout << "Input must be block of original matrix" << endl;

	TComplex a = TComplex::Sqrt(determinant(A, 2));

	vector<vector<TComplex> > V = multiply(1./a, A, 2);

	// func->SetParameter(0, a.Rho()*beta);
	// Double_t x0 = func->GetRandom(-.99999,.99999), x1, x2, x3;

	Double_t x0, x1, x2, x3;

	while(true){
		Double_t r = randGen->Uniform(1), r1 = randGen->Uniform(1), r2 = randGen->Uniform(1), r3 = randGen->Uniform(1);

		Double_t lambda2 = -1./(3./2.*a.Rho()*beta)*(Log(r1)+Cos(2*Pi()*r2)*Cos(2*Pi()*r2)*Log(r3));

		if(r*r > 1 - lambda2) continue;

		x0 = 1-2*lambda2;
		break;
	}

	randGen->Sphere(x1, x2, x3, Sqrt(1.-x0*x0));

	vector<TComplex> row1, row2;
	row1.push_back(TComplex(x0,x3)); row1.push_back(TComplex(x2,x1));
	row2.push_back(TComplex(-x2,x1)); row2.push_back(TComplex(x0,-x3));

	vector<vector<TComplex> > result;
	result.push_back(row1); result.push_back(row2);

	return result;

}

void measure(GaugeField<vector<vector<TComplex> > > *f, Int_t nEntry){

	// Measure the average value of the trace of plaquettes
	Double_t avPlaq = 0, nPlaq = 0;
	Double_t avSpPlaq = 0, nSpPlaq = 0;
	Double_t avTmPlaq = 0, nTmPlaq = 0;
	for(Int_t dir1 = 0; dir1 < 4; dir1++){
		for(Int_t dir2 = dir1+1; dir2 < 4; dir2++){
			for(Int_t nx = 0; nx < Ns; nx++){
				for(Int_t ny = 0; ny < Ns; ny++){
					for(Int_t nz = 0; nz < Ns; nz++){
						for(Int_t nt = 0; nt < Nt; nt++){

							Int_t tempx = nx, tempy = ny, tempz = nz, tempt = nt;

							vector<vector<TComplex> > plaq = (*f)(tempx, tempy, tempz, tempt, dir1);

							if(dir1 == 0) tempx++;
							else if(dir1 == 1) tempy++;
							else if(dir1 == 2) tempz++;
							else if(dir1 == 3) tempt++;

							plaq = multiply(plaq, (*f)(tempx, tempy, tempz, tempt, dir2), 3);

							if(dir2 == 0) tempx++;
							else if(dir2 == 1) tempy++;
							else if(dir2 == 2) tempz++;
							else if(dir2 == 3) tempt++;
							if(dir1 == 0) tempx--;
							else if(dir1 == 1) tempy--;
							else if(dir1 == 2) tempz--;
							else if(dir1 == 3) tempt--;

							plaq = multiply(plaq, dagger((*f)(tempx, tempy, tempz, tempt, dir1)), 3);

							if(dir2 == 0) tempx--;
							else if(dir2 == 1) tempy--;
							else if(dir2 == 2) tempz--;
							else if(dir2 == 3) tempt--;

							plaq = multiply(plaq, dagger((*f)(tempx, tempy, tempz, tempt, dir2)), 3);

							Double_t tr = 1 - 1./3.*trace(plaq).Rho();

							// cout << determinant(plaq, 3) << endl;

							avPlaq += tr; nPlaq++;
							if(dir1 != 3 && dir2 != 3){avSpPlaq += tr; nSpPlaq++;}
							if(dir1 == 3 || dir2 == 3){avTmPlaq += tr; nTmPlaq++;}

						}
					}
				}
			}
		}
	}

	avPlaq /= nPlaq;
	avSpPlaq /= nSpPlaq;
	avTmPlaq /= nTmPlaq;

	entryArray[nEntry] = nEntry;
	plaqArray[nEntry] = avPlaq;

	cout << "Average: " << avPlaq << endl;
	cout << "Average spatial: " << avSpPlaq << endl;
	cout << "Average temporal: " << avTmPlaq << endl << endl;


}

vector<vector<TComplex> > multiply(vector<vector<TComplex> > mat1, vector<vector<TComplex> > mat2, Int_t size){

	vector<vector<TComplex> > result = id;

	if(size == 2){
		result = nu2;
		for(Int_t n = 0; n < 2; n++){
			for(Int_t m = 0; m < 2; m++){
				result[n][m] = mat1[n][0]*mat2[0][m] + mat1[n][1]*mat2[1][m];
			}
		}
	}

	else{
		for(Int_t n = 0; n < 3; n++){
			for(Int_t m = 0; m < 3; m++){
				result[n][m] = mat1[n][0]*mat2[0][m] + mat1[n][1]*mat2[1][m] + mat1[n][2]*mat2[2][m];
			}
		}
	}

	return result;
}

vector<vector<TComplex> > multiply(TComplex s, vector<vector<TComplex> > mat, Int_t size){

	vector<vector<TComplex> > result = id;

	for(UInt_t n = 0; n < size; n++){
		for(UInt_t m = 0; m < size; m++){
			result[n][m] = s*mat[n][m];
		}
	}

	return result;
}

vector<vector<TComplex> > add(vector<vector<TComplex> > mat1, vector<vector<TComplex> > mat2){

	vector<vector<TComplex> > result = id;

	for(Int_t n = 0; n < 3; n++){
		for(Int_t m = 0; m < 3; m++){
			result[n][m] = mat1[n][m] + mat2[n][m];
		}
	}

	return result;
}

vector<vector<TComplex> > dagger(vector<vector<TComplex> > mat){

	vector<vector<TComplex> > result = id;

	for(Int_t n = 0; n < 3; n++){
		for(Int_t m = 0; m < 3; m++){
			result[n][m] = TComplex::Conjugate(mat[m][n]);
		}
	}

	return result;
}

TComplex determinant(vector<vector<TComplex> > mat, Int_t size){
	if(size == 2) return mat[0][0]*mat[1][1] - mat[0][1]*mat[1][0];
	else return mat[0][0]*(mat[1][1]*mat[2][2] - mat[1][2]*mat[2][1]) - mat[0][1]*(mat[1][0]*mat[2][2] - mat[1][2]*mat[2][0]) + mat[0][2]*(mat[1][0]*mat[2][1] - mat[1][1]*mat[2][0]);
}

TComplex trace(vector<vector<TComplex> > mat){
	return mat[0][0] + mat[1][1] + mat[2][2];
}

TGraph* readData(TString dataset){

    TString fileName = dataset + ".txt";

    Int_t nEntries = 0;
    Double_t plaq = 0;
    Double_t plaqArr[1000], entryArr[1000];

    ifstream ifs(fileName); if(!ifs.is_open()){cout << "Error. File " << fileName << " not found. Exiting...\n"; return NULL;}

    while(ifs >> plaq){

        plaqArr[nEntries] = plaq;
        entryArr[nEntries] = nEntries;
        nEntries++;

    }

    ifs.close();

    TGraph *gr = new TGraph(nEntries, entryArr, plaqArr);

    return gr;

}

void graph(TGraph *graph1, TGraph *graph2, const TString graphName, TCanvas *can, const TString xTitle, const TString yTitle, const TString name){

    //gStyle->SetOptStat(2210);
    gStyle->SetOptFit(1111);
    gStyle->SetOptStat(kFALSE);
    //gStyle->SetOptFit(1100);

    if(!graph1 || !graph2){
        cout << "Error: Graph \"" << graphName << "\" not defined" << endl;
        return;
    }

    graph1->SetLineWidth(2);
    graph1->SetMarkerStyle(20);
    graph1->SetMarkerSize(1.5);
    graph1->SetMarkerColor(kBlue);
    graph1->SetLineColor(kBlack);
    graph2->SetLineWidth(2);
    graph2->SetMarkerStyle(20);
    graph2->SetMarkerSize(1.5);
    graph2->SetMarkerColor(kRed);
    graph2->SetLineColor(kBlack);


    gStyle->SetLegendBorderSize(0);
    TLegend *leg = new TLegend(0.655,0.675,0.885,0.875);
    leg->SetTextSize(0.055);
    leg->AddEntry(graph1, "Chroma data","p");
    leg->AddEntry(graph2, "My script's data","p");

    TMultiGraph *mg = new TMultiGraph();

    mg->Add(graph1);
    mg->Add(graph2);

    // graph1->Draw("ap");
    // graph2->Draw("same ap");
    mg->Draw("ap");
    leg->Draw("same");

    // add axis labels
    mg->GetXaxis()->SetTitle(xTitle);
    mg->GetXaxis()->CenterTitle();
    mg->GetXaxis()->SetTitleSize(0.055);
    mg->GetXaxis()->SetTitleOffset(0.82);
    //graph->GetXaxis()->SetLabelOffset(0.010);
    mg->GetXaxis()->SetLabelSize(0.05);
    mg->GetYaxis()->SetTitle(yTitle);
    mg->GetYaxis()->CenterTitle();
    mg->GetYaxis()->SetTitleSize(0.055);
    mg->GetYaxis()->SetTitleOffset(0.9);
    mg->GetYaxis()->SetLabelSize(0.05);
    gStyle->SetTitleSize(0.08, "t");
    mg->SetTitle(graphName);
    // can->Update();

    can->SaveAs(name + ".png");

    //can->Clear();
}
