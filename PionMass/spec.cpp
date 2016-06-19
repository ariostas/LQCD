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
#include <TF2.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <TString.h>
#include "TMath.h"
#include <TPaveStats.h>
#include <TText.h>
#include <TLatex.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TGraphErrors.h>
#include <TMatrixD.h>
#include "TRandom3.h"
#include "Math/WrappedTF1.h"
#include "Math/BrentRootFinder.h"

#endif

using namespace TMath;
using namespace std;

// Declare functions
vector<vector<Double_t> > readData(TString);
vector<Double_t> average(vector<vector<Double_t> >);
vector<Double_t> averageBoot(vector<vector<Double_t> >);
vector<Double_t> errorStd(vector<vector<Double_t> >);
vector<Double_t> errorJack(vector<vector<Double_t> >);
vector<Double_t> errorBoot(vector<vector<Double_t> >);
vector<vector<Double_t> > effMass(vector<vector<Double_t> >);
TMatrixD weightMatrix(vector<vector<Double_t> >);
Double_t findMass(vector<vector<Double_t> >, Int_t, Int_t, Double_t, Double_t);
Double_t findMassBoot(vector<vector<Double_t> >, Int_t, Int_t, Double_t&, Double_t, Double_t);
void graph(vector<TGraphErrors*>, vector<TString>, const TString, const TString, const TString, vector<vector<Double_t> > *mass = 0, TString flags = "");

// Declare common objects
TRandom3 *randGen = new TRandom3();
TF1 massFunc("Mass Function", "cosh(x*(double([1])-double([2])/2.0))/cosh(x*(double([1])+1.0-double([2])/2.0))-double([0])", 0, 5);
ROOT::Math::WrappedTF1 wmf(massFunc);
ROOT::Math::BrentRootFinder brf;
Double_t intArr[100], halfIntArr[100], zeroArr[100];

// Main function
void spec(){

    TH1::StatOverflows(kTRUE);
    
    cout << "\n\nStarting process...\n\n";

    brf.SetFunction(wmf, 0, 5);

    for(Int_t x = 0; x < 100; x++){
        intArr[x] = x; halfIntArr[x] = x+0.5; zeroArr[x] = 0;
    }

    vector<vector<Double_t> > correlators020 = readData("pion.D2000.DG2p8284_1.DG2p8284_1.SS");
    vector<vector<Double_t> > correlators010 = readData("pion.D1000.DG2p8284_1.DG2p8284_1.SS");
    vector<vector<Double_t> > correlators005 = readData("pion.D500.DG2p8284_1.DG2p8284_1.SS");
    vector<vector<Double_t> > correlators002 = readData("pion.D200.DG2p8284_1.DG2p8284_1.SS");

    // vector<vector<Double_t> > correlators020 = readData("pion.D2000.P_1.P_1.PP");
    // vector<vector<Double_t> > correlators010 = readData("pion.D1000.P_1.P_1.PP");
    // vector<vector<Double_t> > correlators005 = readData("pion.D500.P_1.P_1.PP");
    // vector<vector<Double_t> > correlators002 = readData("pion.D200.P_1.P_1.PP");

    vector<Double_t> avCorr020 = averageBoot(correlators020), errCorr020 = errorBoot(correlators020);
    vector<Double_t> avCorr010 = averageBoot(correlators010), errCorr010 = errorBoot(correlators010);
    vector<Double_t> avCorr005 = averageBoot(correlators005), errCorr005 = errorBoot(correlators005);
    vector<Double_t> avCorr002 = averageBoot(correlators002), errCorr002 = errorBoot(correlators002);

    vector<TGraphErrors*> corrGraphs;
    corrGraphs.push_back(new TGraphErrors(avCorr020.size(), intArr, &avCorr020[0], zeroArr, &errCorr020[0]));
    corrGraphs.push_back(new TGraphErrors(avCorr010.size(), intArr, &avCorr010[0], zeroArr, &errCorr010[0]));
    corrGraphs.push_back(new TGraphErrors(avCorr005.size(), intArr, &avCorr005[0], zeroArr, &errCorr005[0]));
    corrGraphs.push_back(new TGraphErrors(avCorr002.size(), intArr, &avCorr002[0], zeroArr, &errCorr002[0]));

    vector<TString> massLabels = {"m = 0.20", "m = 0.10", "m = 0.05", "m = 0.02"};

    graph(corrGraphs, massLabels, "#Deltat", "C(#Deltat)", "corrGraph", 0, "log");

    vector<vector<Double_t> > masses020 = effMass(correlators020);
    vector<vector<Double_t> > masses010 = effMass(correlators010);
    vector<vector<Double_t> > masses005 = effMass(correlators005);
    vector<vector<Double_t> > masses002 = effMass(correlators002);

    vector<Double_t> avMass020 = averageBoot(masses020), errMass020 = errorBoot(masses020);
    vector<Double_t> avMass010 = averageBoot(masses010), errMass010 = errorBoot(masses010);
    vector<Double_t> avMass005 = averageBoot(masses005), errMass005 = errorBoot(masses005);
    vector<Double_t> avMass002 = averageBoot(masses002), errMass002 = errorBoot(masses002);

    vector<TGraphErrors*> massGraphs;
    massGraphs.push_back(new TGraphErrors(avMass020.size(), halfIntArr, &avMass020[0], zeroArr, &errMass020[0]));
    massGraphs.push_back(new TGraphErrors(avMass010.size(), halfIntArr, &avMass010[0], zeroArr, &errMass010[0]));
    massGraphs.push_back(new TGraphErrors(avMass005.size(), halfIntArr, &avMass005[0], zeroArr, &errMass005[0]));
    massGraphs.push_back(new TGraphErrors(avMass002.size(), halfIntArr, &avMass002[0], zeroArr, &errMass002[0]));

    vector<vector<Double_t> > masses(4);
    Double_t m1, m2, m3, m4, e1, e2, e3, e4;
    m1 = findMassBoot(correlators020, 4, 12, e1, 1.99081, 1.34377e-03);
    m2 = findMassBoot(correlators010, 4, 12, e2, 1.89002, 1.51954e-03);
    m3 = findMassBoot(correlators005, 4, 12, e3, 1.83778, 1.61974e-03);
    m4 = findMassBoot(correlators002, 4, 12, e4, 1.80585, 1.68432e-03);

    masses.at(0).push_back(m1); masses.at(0).push_back(e1); masses.at(0).push_back(4); masses.at(0).push_back(12);
    masses.at(1).push_back(m2); masses.at(1).push_back(e2); masses.at(1).push_back(4); masses.at(1).push_back(12);
    masses.at(2).push_back(m3); masses.at(2).push_back(e3); masses.at(2).push_back(4); masses.at(2).push_back(12);
    masses.at(3).push_back(m4); masses.at(3).push_back(e4); masses.at(3).push_back(4); masses.at(3).push_back(12);

    graph(massGraphs, massLabels, "#Deltat", "m_{eff}(#Deltat)", "massGraph", &masses);


    // cout << m1 << " +- " << e1 << endl;
    // cout << m2 << " +- " << e2 << endl;
    // cout << m3 << " +- " << e3 << endl;
    // cout << m4 << " +- " << e4 << endl;
    // findMass(correlators010, 4, 14, mass, err);
    // findMass(correlators005, 4, 14, mass, err);
    // findMass(correlators002, 4, 14, mass, err);

    // TH1D *hist = new TH1D("corr", "corr", 100, 0, 1);

    // for(Int_t x = 0; x < correlators020.size(); x++){
    //     hist->Fill(Log(correlators020.at(x).at(10)));
    // }

    // TCanvas *can = new TCanvas("sds", "sds", 1300, 900);

    // hist->Draw();

    TMatrixD w1 = weightMatrix(correlators020);
    for(Int_t x = 0; x < 32; x++){
        for(Int_t y = 0; y < 32; y++){
            cout << w1(x,y) << "  ";
        }
        cout << endl;
        
    }

}

// Read correlator data from output of strip_hadspec
vector<vector<Double_t> > readData(TString filename){

    cout << "Reading " << filename << "..." << endl;

    ifstream ifs(filename); if(!ifs.is_open()){cout << "Error. File " << filename << " not found. Exiting...\n"; assert(0);}

    Int_t NConfigs = 0, Nt = 0, IsComplex = 0, Ns = 0, temp = 0;
    ifs >> NConfigs >> Nt >> IsComplex >> Ns >> temp;

    vector<vector<Double_t> > corr(NConfigs);

    if(NConfigs == 0){cout << "Error reading the file. Exiting..." << endl; return corr;}

    Int_t deltaT = 0;
    Double_t re = 0, im = 0;
    
    Int_t configCounter = -1;
    while(ifs >> deltaT >> re >> im){

        if(deltaT == 0) configCounter++;

        corr.at(configCounter).push_back(re);

    }

    ifs.close();

    // Check integrity of the data
    for(UInt_t x = 0; x < corr.size(); x++){
        if(corr.at(x).size() != UInt_t(Nt)) cout << "Error: Number of correlators different to temporal extent\nCheck the code\n";
    }

    return corr;

}

// Compute average over all configurations
vector<Double_t> average(vector<vector<Double_t> > corr){
    vector<Double_t> av;
    for(UInt_t x = 0; x < corr.at(0).size(); x++){
        Double_t sum = 0, N = corr.size();
        for(UInt_t y = 0; y < corr.size(); y++){
            sum += corr.at(y).at(x);
        }
        av.push_back(sum/N);
    }
    return av;
}

// Compute average over all configurations (with bootstrapping)
vector<Double_t> averageBoot(vector<vector<Double_t> > corr){
    vector<Double_t> av(corr.at(0).size());
    Double_t N = corr.size();
    for(Int_t x = 0; x < 5*N; x++){
        vector<vector<Double_t> > temp;
        for(Int_t y = 0; y < N; y++){
            temp.push_back(corr.at(randGen->Integer(N)));
        }
        vector<Double_t> tempAv = average(temp);
        for(UInt_t y = 0; y < av.size(); y++){
            av.at(y) += tempAv.at(y);
        }
        temp.clear();
    }
    for(UInt_t y = 0; y < av.size(); y++){
        av.at(y) = 1./(5.*N)*av.at(y);
    }    
    return av;
}

// Compute error (standard deviation)
vector<Double_t> errorStd(vector<vector<Double_t> > corr){
    vector<Double_t> av = average(corr), err;
    for(UInt_t x = 0; x < corr.at(0).size(); x++){
        Double_t sum = 0, N = corr.size();
        for(UInt_t y = 0; y < corr.size(); y++){
            sum += (corr.at(y).at(x)-av.at(x))*(corr.at(y).at(x)-av.at(x));
        }
        err.push_back(Sqrt(1./N/(N-1.)*sum));
    }
    return err;
}

// Compute error (Jacknife)
vector<Double_t> errorJack(vector<vector<Double_t> > corr){
    vector<Double_t> av = average(corr), err(corr.at(0).size());
    Double_t N = corr.size();
    for(UInt_t x = 0; x < corr.size(); x++){
        vector<vector<Double_t> > temp;
        for(UInt_t y = 0; y < corr.size(); y++){
            if(x == y) continue;
            temp.push_back(corr.at(x));
        }
        vector<Double_t> tempAv = average(temp);
        for(UInt_t y = 0; y < err.size(); y++){
            err.at(y) += (tempAv.at(y)-av.at(y))*(tempAv.at(y)-av.at(y));
        }
        temp.clear();
    }
    for(UInt_t y = 0; y < err.size(); y++){
        err.at(y) = Sqrt((N-1.)/N*err.at(y));
    }    
    return err;
}

// Compute error (Bootstrapping)
vector<Double_t> errorBoot(vector<vector<Double_t> > corr){
    vector<Double_t> av = averageBoot(corr), err(corr.at(0).size());
    Double_t N = corr.size();
    for(Int_t x = 0; x < 5*N; x++){
        vector<vector<Double_t> > temp;
        for(Int_t y = 0; y < N; y++){
            temp.push_back(corr.at(randGen->Integer(N)));
        }
        vector<Double_t> tempAv = average(temp);
        for(UInt_t y = 0; y < err.size(); y++){
            err.at(y) += (tempAv.at(y)-av.at(y))*(tempAv.at(y)-av.at(y));
        }
        temp.clear();
    }
    for(UInt_t y = 0; y < err.size(); y++){
        err.at(y) = Sqrt(1./(5.*N)*err.at(y));
    }    
    return err;
}

vector<vector<Double_t> > effMass(vector<vector<Double_t> > corr){
    vector<vector<Double_t> > masses;
    vector<Double_t> temp;
    for(UInt_t x = 0; x < corr.size();x++){
        masses.push_back(temp);
    }
    massFunc.FixParameter(2, corr.at(0).size());
    for(UInt_t x = 0; x < corr.size();x++){
        for(UInt_t y = 0; y < corr.at(0).size()-1; y++){
            massFunc.FixParameter(1, y);
            massFunc.FixParameter(0, corr.at(x).at(y)/corr.at(x).at(y+1));
            // masses.at(x).push_back(massFunc->GetX(0., 0., 2., 1e-15, 100000));
            // masses.at(x).push_back(Log(corr.at(x).at(y)/corr.at(x).at(y+1)));
            brf.Solve();
            masses.at(x).push_back(brf.Root());
        }
    }
    return masses;
}

Double_t findMass(vector<vector<Double_t> > corr, Int_t nmin, Int_t nmax, Double_t massGuess, Double_t AGuess){

    vector<Double_t> av = averageBoot(corr);
    Double_t N = av.size();
    TMatrixD w = weightMatrix(corr);

    TString func;

    for(Int_t x = nmin; x <= nmax; x++){
        for(Int_t y = nmin; y <= nmax; y++){
            func.Append(TString::Format("+(%1.6g-2.0*y*exp(-%1.1f*x/2.0)*cosh((%1.1f/2.0-%1.1f)*x))", av.at(x), N, N, Double_t(x)));
            func.Append(TString::Format("*%1.6g", w(x,y)));
            func.Append(TString::Format("*(%1.6g-2.0*y*exp(-%1.1f*x/2.0)*cosh((%1.1f/2.0-%1.1f)*x))", av.at(y), N, N, Double_t(y)));
        }
    }

    TF2 *fitFunc = new TF2("Fit function", func, massGuess-0.5, massGuess+0.5, AGuess-0.005, AGuess+0.005);

    Double_t mass, A;

    fitFunc->GetMinimumXY(mass, A);

    delete fitFunc;

    // cout << mass << "  " << A << endl;

    return mass;
    
}

Double_t findMassBoot(vector<vector<Double_t> > corr, Int_t nmin, Int_t nmax, Double_t &err, Double_t massGuess, Double_t AGuess){

    Double_t N = corr.size(), NB = 2*N;
    vector<Double_t> tempMasses;
    for(Int_t x = 0; x < NB; x++){
        if(x == 0) cout << TString::Format("Performing bootstrapping on masses... %3.0f\%\n", Double_t(x)/Double_t(NB)*100.);
        else cout << TString::Format("\e[APerforming bootstrapping on masses... %3.0f\%\n", Double_t(x)/Double_t(NB)*100.);
        if(x == NB-1) cout << TString::Format("\e[APerforming bootstrapping on masses... \033[1;32mdone\033[0m\n");
        vector<vector<Double_t> > temp;
        for(Int_t y = 0; y < N; y++){
            temp.push_back(corr.at(randGen->Integer(N)));
        }
        tempMasses.push_back(findMass(temp, nmin, nmax, massGuess, AGuess));
        temp.clear();
    }
    Double_t mass = 0;
    for(UInt_t x = 0; x < NB; x++){
        mass += tempMasses.at(x);
    }
    mass /= NB;
    err = 0;
    for(UInt_t x = 0; x < NB; x++){
        err += (tempMasses.at(x)-mass)*(tempMasses.at(x)-mass);
    }
    err = Sqrt(1./NB*err);

    return mass;

}

TMatrixD weightMatrix(vector<vector<Double_t> > corr){

    Int_t n = corr.at(0).size();
    Double_t N = corr.size();
    vector<Double_t> av = average(corr);
    TMatrixD covMatrix(n, n);

    for(Int_t x = 0; x < n; x++){
        for(Int_t y = 0; y < n; y++){
            for(Int_t z = 0; z < N; z++){
                covMatrix(x,y) += (corr.at(z).at(x)-av.at(x))*(corr.at(z).at(y)-av.at(y))*1e25;
            }
            covMatrix(x,y) = 1.0/((N-1.0)*N)*covMatrix(x,y);
        }
    }

    TMatrixD inv = covMatrix.Invert();
    inv *= 1e-25;

    return inv;

}

void graph(vector<TGraphErrors*> graphs, vector<TString> massLabels, const TString xTitle, const TString yTitle, const TString name, vector<vector<Double_t> > *mass, TString flags){

    //gStyle->SetOptStat(2210);
    gStyle->SetOptFit(1111);
    gStyle->SetOptStat(kFALSE);
    //gStyle->SetOptFit(1100);

    if(graphs.size() != massLabels.size()){cout << "Error: Size of vectors does not match." << endl; return;}

    TCanvas *can = new TCanvas(name, name, 1300, 900);

    gPad->SetLogy( (flags.Contains("log") ? 1 : 0) );

    vector<TF1*> massFit;
    if(mass != 0){
        cout << mass << endl;
        for(UInt_t x = 0; x < mass->size(); x++){
            massFit.push_back(new TF1(TString::Format("mass name%i",x), TString::Format("%1.6f", mass->at(x).at(0)), mass->at(x).at(2), mass->at(x).at(3)));
            massFit.push_back(new TF1(TString::Format("mass+ name%i",x), TString::Format("%1.6f", mass->at(x).at(0)+mass->at(x).at(1)), mass->at(x).at(2), mass->at(x).at(3)));
            massFit.push_back(new TF1(TString::Format("mass- name%i",x), TString::Format("%1.6f", mass->at(x).at(0)-mass->at(x).at(1)), mass->at(x).at(2), mass->at(x).at(3)));
        }
        
    }

    vector<Int_t> colors, markers;
    colors.push_back(kRed); colors.push_back(kBlack); colors.push_back(kBlue); colors.push_back(kGreen+2); colors.push_back(kRed); colors.push_back(kBlack); 
    markers.push_back(23); markers.push_back(22); markers.push_back(33); markers.push_back(21); markers.push_back(20);

    for(UInt_t x = 0; x < graphs.size(); x++){
        graphs.at(x)->SetLineWidth(2);
        graphs.at(x)->SetMarkerStyle(markers.at(x%5));
        graphs.at(x)->SetMarkerSize(1.5);
        graphs.at(x)->SetMarkerColor(colors.at(x%6));
        graphs.at(x)->SetLineColor(colors.at(x%6));

        if(mass != 0){
            massFit.at(3*x)->SetLineColor(colors.at(x%6));
            massFit.at(3*x+1)->SetLineColor(colors.at(x%6));
            massFit.at(3*x+2)->SetLineColor(colors.at(x%6));
            massFit.at(3*x+1)->SetLineStyle(7);
            massFit.at(3*x+2)->SetLineStyle(7);
        }
    }

    TMultiGraph *mg = new TMultiGraph();
    mg->SetTitle("");

    gStyle->SetLegendBorderSize(0);
    TLegend *leg = new TLegend(0.655,0.675,0.885,0.875);
    leg->SetTextSize(0.055);
    for(UInt_t x = 0; x < graphs.size(); x++){
        leg->AddEntry(graphs.at(graphs.size()-x-1), massLabels.at(graphs.size()-x-1),"lep");
        mg->Add(graphs.at(x));
    }
    
    mg->Draw("acp");
    leg->Draw("same");
    if(mass != 0){
        for(UInt_t x = 0; x < massFit.size(); x++){
            massFit.at(x)->Draw("same");
        }
    }

    // add axis labels
    mg->GetXaxis()->SetTitle(xTitle);
    mg->GetXaxis()->CenterTitle();
    mg->GetXaxis()->SetTitleSize(0.055);
    mg->GetXaxis()->SetTitleOffset(0.85);
    //graph->GetXaxis()->SetLabelOffset(0.010);
    mg->GetXaxis()->SetLabelSize(0.05);
    mg->GetYaxis()->SetTitle(yTitle);
    mg->GetYaxis()->CenterTitle();
    mg->GetYaxis()->SetTitleSize(0.055);
    mg->GetYaxis()->SetTitleOffset(0.9);
    mg->GetYaxis()->SetLabelSize(0.05);
    gStyle->SetTitleSize(0.08, "t");
    mg->SetTitle("");
    // can->Update();

    // TF1 *f1 = new TF1(name, "2.0*[1]*TMath::Exp(-16.*[0])*TMath::CosH((16.-x)*[0])", 4, 10);
    // f1->SetParameter(0, 2);
    // f1->SetParameter(1, 1e-3);
    // f1->SetLineColor(kBlack);
    // graphs.at(2)->Fit(f1, "ME", "", 4, 10);
    // graphs.at(3)->Fit(f1, "ME", "", 4, 10);
    // cout << f1->GetParameter(0) << endl;

    can->SaveAs(name + ".png");

    //can->Clear();
}