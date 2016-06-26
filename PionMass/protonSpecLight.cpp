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
vector<Double_t> effMassBoot(vector<vector<Double_t> >, vector<Double_t>&);
TMatrixD weightMatrix(vector<vector<Double_t> >);
Double_t findMass(vector<vector<Double_t> >, Int_t, Int_t, Double_t, Double_t, Double_t *chi2 = 0);
Double_t findMassBoot(vector<vector<Double_t> >, Int_t, Int_t, Double_t&, Double_t, Double_t, Double_t *chi2 = 0, Double_t *chi2err = 0);
void scanFitRange(vector<vector<Double_t> >);
void graph(vector<TGraphErrors*>, vector<TString>, const TString, const TString, const TString, vector<vector<Double_t> > *mass = 0, TString flags = "");

// Declare common objects
TRandom3 *randGen = new TRandom3();
Double_t intArr[100], halfIntArr[100], zeroArr[100];

// Main function
void protonSpecLight(){

    TH1::StatOverflows(kTRUE);
    
    cout << "\n\nStarting process...\n\n";

    for(Int_t x = 0; x < 100; x++){
        intArr[x] = x; halfIntArr[x] = x+0.5; zeroArr[x] = 0;
    }

    vector<vector<Double_t> > correlators092 = readData("proton.D-9199.DG2_1.DG2_1.SS");

    vector<Double_t> avCorr092 = averageBoot(correlators092), errCorr092 = errorBoot(correlators092);

    vector<TGraphErrors*> corrGraphs;
    corrGraphs.push_back(new TGraphErrors(avCorr092.size(), intArr, &avCorr092[0], zeroArr, &errCorr092[0]));

    vector<TString> massLabels = {"m = -0.92"};

    graph(corrGraphs, massLabels, "#Deltat", "C(#Deltat)", "protonCorrGraph", 0, "log");

    vector<Double_t> errMass092, avMass092 = effMassBoot(correlators092, errMass092);

    vector<TGraphErrors*> massGraphs;
    massGraphs.push_back(new TGraphErrors(avMass092.size(), halfIntArr, &avMass092[0], zeroArr, &errMass092[0]));
    vector<vector<Double_t> > masses(1);
    Double_t m1, e1;
    m1 = findMassBoot(correlators092, 5, 14, e1, 7.33767e-01, 5.11457e-05);

    masses.at(0).push_back(m1); masses.at(0).push_back(e1); masses.at(0).push_back(5); masses.at(0).push_back(14);

    graph(massGraphs, massLabels, "#Deltat", "m_{eff}(#Deltat)", "protonMassGraph", &masses);

    // scanFitRange(correlators092);

}

// Read correlator data from output of strip_hadspec
vector<vector<Double_t> > readData(TString filename){

    cout << "Reading " << filename << "..." << endl;

    const TString fullFilename = "./VariationalLight/" + filename;

    ifstream ifs(fullFilename); if(!ifs.is_open()){cout << "Error. File " << fullFilename << " not found. Exiting...\n"; assert(0);}

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

vector<Double_t> effMassBoot(vector<vector<Double_t> > corr, vector<Double_t> &err){

    Double_t N = corr.size(), NB = 5.*N;
    vector<vector<Double_t> > allMasses;
    for(Int_t x = 0; x < NB; x++){
        vector<vector<Double_t> > tempcorr;
        for(Int_t y = 0; y < N; y++){
            tempcorr.push_back(corr.at(randGen->Integer(N)));
        }
        vector<Double_t> tempmass, tempAv = average(tempcorr);
        for(UInt_t y = 0; y < tempAv.size()-1; y++){
            Double_t m = Log(tempAv.at(y)/tempAv.at(y+1));
            if(y > Double_t(tempAv.size())/2.) m *= -1.;
            tempmass.push_back(m);

            // massFunc.FixParameter(1, y);
            // massFunc.FixParameter(0, tempAv.at(y)/tempAv.at(y+1));
            // brf.Solve();
            // tempmass.push_back(brf.Root());
        }
        tempcorr.clear();
        tempAv.clear();
        allMasses.push_back(tempmass);
        tempmass.clear();
    }

    vector<Double_t> masses(corr.at(0).size()-1);
    for(UInt_t x = 0; x < NB; x++){
        for(UInt_t y = 0; y < allMasses.at(0).size(); y++){
            masses.at(y) += allMasses.at(x).at(y);
        }
    }
    for(UInt_t x = 0; x < masses.size(); x++){
        masses.at(x) /= NB;
    }

    vector<Double_t> error(corr.at(0).size()-1);
    for(UInt_t x = 0; x < NB; x++){
        for(UInt_t y = 0; y < allMasses.at(0).size(); y++){
            error.at(y) += (allMasses.at(x).at(y)-masses.at(y))*(allMasses.at(x).at(y)-masses.at(y));
        }
    }
    for(UInt_t x = 0; x < error.size(); x++){
        error.at(x) = Sqrt(1./NB*error.at(x));
    }

    err = error;
    return masses;
}

Double_t findMass(vector<vector<Double_t> > corr, Int_t nmin, Int_t nmax, Double_t massGuess, Double_t AGuess, Double_t *chi2){

    vector<Double_t> av = averageBoot(corr);
    Double_t N = av.size();
    TMatrixD w = weightMatrix(corr);

    TString func;

    for(Int_t x = nmin; x <= nmax; x++){
        for(Int_t y = nmin; y <= nmax; y++){
            func.Append(TString::Format("+(%1.6g-y*TMath::Exp(-%1.1f*x))", av.at(x), float(x)));
            func.Append(TString::Format("*%1.6g", w(x,y)));
            func.Append(TString::Format("*(%1.6g-y*TMath::Exp(-%1.1f*x))", av.at(y), float(y)));
            // func.Append(TString::Format("+(%1.6g-2.0*y*exp(-%1.1f*x/2.0)*cosh((%1.1f/2.0-%1.1f)*x))", av.at(x), N, N, Double_t(x)));
            // func.Append(TString::Format("*%1.6g", w(x,y)));
            // func.Append(TString::Format("*(%1.6g-2.0*y*exp(-%1.1f*x/2.0)*cosh((%1.1f/2.0-%1.1f)*x))", av.at(y), N, N, Double_t(y)));
        }
    }

    // cout << func << endl;

    TF2 *fitFunc = new TF2("Fit function", func, massGuess-0.5, massGuess+0.5, AGuess-0.0005, AGuess+0.0005);

    Double_t mass, A;

    fitFunc->GetMinimumXY(mass, A);

    if(chi2 != 0) (*chi2) = fitFunc->Eval(mass, A)/Double_t(nmax - nmin - 1);

    // cout << "chi2: " << fitFunc->Eval(mass, A) << endl;

    delete fitFunc;

    // massFunc.FixParameter(1, nmin);
    // massFunc.FixParameter(0, av.at(nmin)/av.at(nmin+1));
    // brf.Solve();

    // cout << mass << "  " << brf.Root() << endl;

    return mass;
    
}

Double_t findMassBoot(vector<vector<Double_t> > corr, Int_t nmin, Int_t nmax, Double_t &err, Double_t massGuess, Double_t AGuess, Double_t *chi2, Double_t *chi2err){

    Double_t N = corr.size(), NB = Nint(N/4);
    vector<Double_t> tempMasses, tempchi2(NB);
    for(Int_t x = 0; x < NB; x++){
        if(x == 0) cout << TString::Format("Performing bootstrap on masses... %3.0f\%\n", Double_t(x)/Double_t(NB)*100.);
        else cout << TString::Format("\e[APerforming bootstrap on masses... %3.0f\%\n", Double_t(x)/Double_t(NB)*100.);
        if(x == NB-1) cout << TString::Format("\e[APerforming bootstrap on masses... \033[1;32mdone\033[0m\n");
        vector<vector<Double_t> > temp;
        for(Int_t y = 0; y < N; y++){
            temp.push_back(corr.at(randGen->Integer(N)));
        }
        tempMasses.push_back(findMass(temp, nmin, nmax, massGuess, AGuess, &(tempchi2.at(x))));
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

    Double_t c2 = 0, c2e = 0;
    for(UInt_t x = 0; x < NB; x++){
        c2 += tempchi2.at(x);
    }
    c2 /= NB;

    for(UInt_t x = 0; x < NB; x++){
        c2e += (tempchi2.at(x)-c2)*(tempchi2.at(x)-c2);
    }
    c2e = Sqrt(1./NB*c2e);

    if(chi2 != 0) (*chi2) = c2;
    if(chi2err != 0) (*chi2err) = c2e;

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
                covMatrix(x,y) += (corr.at(z).at(x)-av.at(x))*(corr.at(z).at(y)-av.at(y));
            }
            covMatrix(x,y) = 1.0/((N-1.0)*N)*covMatrix(x,y)*1e30;
        }
    }

    TMatrixD inv = covMatrix.Invert();

    inv *= 1e30;
    return inv;

}

void scanFitRange(vector<vector<Double_t> > corr){

    Double_t arr[100];
    for(Int_t x = 1; x < 101; x++){
        arr[x-1] = x;
    }

    Int_t N = 11;
    vector<Double_t> m(N), me(N), c2(N), c2e(N);

    for(Int_t x = 0; x < N; x++){
        m.at(x) = findMassBoot(corr, 12-x, 14, me.at(x), 7.33767e-01, 5.11457e-05, &(c2.at(x)), &(c2e.at(x)));
    }

    TGraphErrors *massScan = new TGraphErrors(N, arr, &m[0], zeroArr, &me[0]);
    TGraphErrors *chiScan = new TGraphErrors(N, arr, &c2[0], zeroArr, &c2e[0]);

    for(UInt_t x = 0; x < m.size(); x++){
        cout << m.at(x) << "  " << c2.at(x) << endl;
    }

    TCanvas *can = new TCanvas("scan", "scan", 1300, 900);

    can->Divide(1,2);
    can->cd(1);
    massScan->SetLineWidth(2);
    massScan->SetMarkerStyle(21);
    massScan->SetTitle("");
    massScan->GetXaxis()->SetTitle("Fit range (#Deltat away from Nt/2)");
    massScan->GetXaxis()->CenterTitle();
    massScan->GetXaxis()->SetLabelSize(0.04);
    massScan->GetXaxis()->SetTitleSize(0.05);
    massScan->GetYaxis()->SetTitle("Mass");
    massScan->GetYaxis()->CenterTitle();
    massScan->GetYaxis()->SetLabelSize(0.04);
    massScan->GetYaxis()->SetTitleSize(0.05);
    chiScan->SetLineWidth(2);
    chiScan->SetMarkerStyle(21);
    chiScan->SetTitle("");
    chiScan->GetXaxis()->SetTitle("Fit range (#Deltat away from Nt/2)");
    chiScan->GetXaxis()->CenterTitle();
    chiScan->GetXaxis()->SetLabelSize(0.04);
    chiScan->GetXaxis()->SetTitleSize(0.05);
    chiScan->GetYaxis()->SetTitle("#chi^{2}/d.o.f");
    chiScan->GetYaxis()->CenterTitle();
    chiScan->GetYaxis()->SetLabelSize(0.04);
    chiScan->GetYaxis()->SetTitleSize(0.05);

    massScan->Draw("ap");
    can->cd(2);
    chiScan->Draw("ap");

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
        // cout << mass << endl;
        for(UInt_t x = 0; x < mass->size(); x++){
            massFit.push_back(new TF1(TString::Format("mass_%s%i", name.Data(), x), TString::Format("%1.6f", mass->at(x).at(0)), mass->at(x).at(2), mass->at(x).at(3)));
            massFit.push_back(new TF1(TString::Format("mass1_%s%i", name.Data(), x), TString::Format("%1.6f", mass->at(x).at(0)+mass->at(x).at(1)), mass->at(x).at(2), mass->at(x).at(3)));
            massFit.push_back(new TF1(TString::Format("mass2_%s%i", name.Data(), x), TString::Format("%1.6f", mass->at(x).at(0)-mass->at(x).at(1)), mass->at(x).at(2), mass->at(x).at(3)));
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
    // f1->SetParameter(0, 0.7);
    // f1->SetParameter(1, 1e-1);
    // f1->SetLineColor(kBlack);
    // graphs.at(0)->Fit(f1, "ME", "", 5, 12);
    // cout << f1->GetParameter(0) << endl;

    can->SaveAs(name + ".png");
    can->SaveAs("./Plots/" + name + ".cpp");

    //can->Clear();
}
