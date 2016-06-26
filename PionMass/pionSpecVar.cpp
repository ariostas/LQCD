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
#include <TVectorD.h>
#include "TRandom3.h"
#include "Math/WrappedTF1.h"
#include "Math/BrentRootFinder.h"

#endif

using namespace TMath;
using namespace std;

// Declare functions
vector<vector<vector<Double_t> > > readData(TString);
vector<vector<Double_t> > average(vector<vector<vector<Double_t> > >);
vector<vector<Double_t> > averageBoot(vector<vector<vector<Double_t> > >);
vector<vector<Double_t> > errorBoot(vector<vector<vector<Double_t> > >);
vector<TMatrixD> constructMatrices(vector<vector<Double_t> >);
vector<vector<Double_t> > findEigenvalues(vector<TMatrixD>);
vector<vector<Double_t> > findEvalBoot(vector<vector<vector<Double_t> > >, vector<vector<Double_t> >&);
vector<vector<Double_t> > findEffMassesBoot(vector<vector<vector<Double_t> > >, vector<vector<Double_t> >&);
vector<TGraphErrors*> makeGraphs(vector<vector<Double_t> >, vector<vector<Double_t> >, TString pos = "int");
vector<vector<Double_t> > effMass(vector<vector<Double_t> >);
// TMatrixD weightMatrices(vector<vector<Double_t> >);
void graph(vector<TGraphErrors*>, vector<TString>, const TString, const TString, const TString, vector<vector<Double_t> > *mass = 0, TString flags = "");

// Declare common objects
TRandom3 *randGen = new TRandom3();
TF1 massFunc("Mass Function", "cosh(x*(double([1])-double([2])/2.0))/cosh(x*(double([1])+1.0-double([2])/2.0))-double([0])", 0, 5);
ROOT::Math::WrappedTF1 wmf(massFunc);
ROOT::Math::BrentRootFinder brf;
Double_t intArr[100], halfIntArr[100], zeroArr[100];
vector<TString> corrTypes = {"P_1.P_1.PP", "DG2_1.P_1.SP", "DG4_1.P_1.SP", "P_1.DG2_1.PS", "DG2_1.DG2_1.SS", "DG4_1.DG2_1.SS",
                             "P_1.DG4_1.PS", "DG2_1.DG4_1.SS", "DG4_1.DG4_1.SS"};

// Main function
void pionSpecVar(){

    TH1::StatOverflows(kTRUE);
    
    cout << "\n\nStarting process...\n\n";

    brf.SetFunction(wmf, 0, 5);

    for(Int_t x = 0; x < 100; x++){
        intArr[x] = x; halfIntArr[x] = x+0.5; zeroArr[x] = 0;
    }

    vector<vector<vector<Double_t> > > correlators020 = readData("2000");
    vector<vector<vector<Double_t> > > correlators010 = readData("1000");
    vector<vector<vector<Double_t> > > correlators005 = readData("500");
    vector<vector<vector<Double_t> > > correlators002 = readData("200");

    vector<vector<Double_t> > evalErrors020, evals020 = findEvalBoot(correlators020, evalErrors020);
    vector<vector<Double_t> > evalErrors010, evals010 = findEvalBoot(correlators010, evalErrors010);
    vector<vector<Double_t> > evalErrors005, evals005 = findEvalBoot(correlators005, evalErrors005);
    vector<vector<Double_t> > evalErrors002, evals002 = findEvalBoot(correlators002, evalErrors002);

    vector<vector<Double_t> > massErrors020, masses020 = findEffMassesBoot(correlators020, massErrors020);
    vector<vector<Double_t> > massErrors010, masses010 = findEffMassesBoot(correlators010, massErrors010);
    vector<vector<Double_t> > massErrors005, masses005 = findEffMassesBoot(correlators005, massErrors005);
    vector<vector<Double_t> > massErrors002, masses002 = findEffMassesBoot(correlators002, massErrors002);

    vector<TGraphErrors*> graphs020 = makeGraphs(evals020, evalErrors020);
    vector<TGraphErrors*> graphs010 = makeGraphs(evals010, evalErrors010);
    vector<TGraphErrors*> graphs005 = makeGraphs(evals005, evalErrors005);
    vector<TGraphErrors*> graphs002 = makeGraphs(evals002, evalErrors002);

    vector<TGraphErrors*> massgraphs020 = makeGraphs(masses020, massErrors020);
    vector<TGraphErrors*> massgraphs010 = makeGraphs(masses010, massErrors010);
    vector<TGraphErrors*> massgraphs005 = makeGraphs(masses005, massErrors005);
    vector<TGraphErrors*> massgraphs002 = makeGraphs(masses002, massErrors002);

    vector<TString> names = {"Ground state", "1st excited state", "Remaining states"};

    graph(graphs020, names, "#Deltat", "Eigenvalues(#Deltat)", "pionVarCorrs020", 0, "log");
    graph(graphs010, names, "#Deltat", "Eigenvalues(#Deltat)", "pionVarCorrs010", 0, "log");
    graph(graphs005, names, "#Deltat", "Eigenvalues(#Deltat)", "pionVarCorrs005", 0, "log");
    graph(graphs002, names, "#Deltat", "Eigenvalues(#Deltat)", "pionVarCorrs002", 0, "log");

    graph(massgraphs020, names, "#Deltat", "M(#Deltat)", "pionVarMass020");
    graph(massgraphs010, names, "#Deltat", "M(#Deltat)", "pionVarMass010");
    graph(massgraphs005, names, "#Deltat", "M(#Deltat)", "pionVarMass005");
    graph(massgraphs002, names, "#Deltat", "M(#Deltat)", "pionVarMass002");

}

// Read correlator data from output of strip_hadspec
vector<vector<vector<Double_t> > > readData(TString mass){

    cout << "Reading correlators with quark mass " << mass << "..." << endl;

    vector<vector<vector<Double_t> > > corrs;

    for(UInt_t x = 0; x < corrTypes.size(); x++){

        TString filename = "./VariationalData/pion.D" + mass + "." + corrTypes.at(x);

        ifstream ifs(filename); if(!ifs.is_open()){cout << "Error. File " << filename << " not found. Exiting...\n"; assert(0);}

        Int_t NConfigs = 0, Nt = 0, IsComplex = 0, Ns = 0, temp = 0;
        ifs >> NConfigs >> Nt >> IsComplex >> Ns >> temp;

        vector<vector<Double_t> > corr(NConfigs);

        if(NConfigs == 0){cout << "Error reading the file. Exiting..." << endl; assert(0);}

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

        corrs.push_back(corr);
        corr.clear();

    }

    return corrs;

}

// Compute average over all configurations
vector<vector<Double_t> > average(vector<vector<vector<Double_t> > > corr){
    vector<vector<Double_t> > av(corr.size());
    for(UInt_t x = 0; x < corr.size(); x++){
        for(UInt_t y = 0; y < corr.at(0).at(0).size(); y++){
            Double_t sum = 0, N = corr.at(x).size();
            for(UInt_t z = 0; z < corr.at(x).size(); z++){
                sum += corr.at(x).at(z).at(y);
            }
            av.at(x).push_back(sum/N);
        }
    }
    return av;
}

// Compute average over all configurations (with bootstrapping)
vector<vector<Double_t> > averageBoot(vector<vector<vector<Double_t> > > corr){
    vector<vector<Double_t> > av(corr.size(),vector<Double_t>(corr.at(0).at(0).size()));
    Double_t N = corr.at(0).size(), NB = 5*N;
    for(Int_t y = 0; y < NB; y++){
        vector<vector<vector<Double_t> > > tempcorr(corr.size());
        for(Int_t z = 0; z < N; z++){
            Int_t num = randGen->Integer(N);
            for(UInt_t w = 0; w < corr.size(); w++){
                tempcorr.at(w).push_back(corr.at(w).at(num));
            }
        }
        vector<vector<Double_t> > tempAv = average(tempcorr);
        for(UInt_t z = 0; z < av.at(0).size(); z++){
            for(UInt_t w = 0; w < corr.size(); w++){
                av.at(w).at(z) += tempAv.at(w).at(z);
            }
        }
        tempcorr.clear();
        tempAv.clear();
    }
    for(UInt_t z = 0; z < av.at(0).size(); z++){
        for(UInt_t w = 0; w < corr.size(); w++){
            av.at(w).at(z) /= NB;
        }
    }
    return av;
}

// Compute error (Bootstrapping)
vector<vector<Double_t> > errorBoot(vector<vector<vector<Double_t> > > corr){
    vector<vector<Double_t> > av = averageBoot(corr), err(corr.size(),vector<Double_t>(corr.at(0).at(0).size()));
    Double_t N = corr.at(0).size(), NB = 5*N;
    for(Int_t y = 0; y < NB; y++){
        vector<vector<vector<Double_t> > > tempcorr;
        for(Int_t z = 0; z < N; z++){
            Int_t num = randGen->Integer(N);
            for(UInt_t w = 0; w < corr.size(); w++){
                tempcorr.at(w).push_back(corr.at(w).at(num));
            }
        }
        vector<vector<Double_t> > tempAv = average(tempcorr);
        for(UInt_t z = 0; z < err.at(0).size(); z++){
            for(UInt_t w = 0; w < corr.size(); w++){
                err.at(w).at(z) += (tempAv.at(w).at(z)-av.at(w).at(z))*(tempAv.at(w).at(z)-av.at(w).at(z));
            }
        }
        tempcorr.clear();
        tempAv.clear();
    }
    for(UInt_t z = 0; z < err.at(0).size(); z++){
        for(UInt_t w = 0; w < corr.size(); w++){
            err.at(w).at(z) = Sqrt(1./NB*err.at(w).at(z));
        }
    }
    return err;
}

vector<TMatrixD> constructMatrices(vector<vector<Double_t> > avcorr){

    UInt_t matSize = Sqrt(avcorr.size());
    vector<TMatrixD> matrices(avcorr.at(0).size(),TMatrixD(matSize,matSize));

    for(UInt_t x = 0; x < avcorr.at(0).size(); x++){
        for(UInt_t y = 0; y < avcorr.size(); y++){
            matrices.at(x)(y%matSize,y/matSize) = avcorr.at(y).at(x);
        }
    }

    return matrices;

}

vector<vector<Double_t> > findEigenvalues(vector<TMatrixD> mat){

    vector<vector<Double_t> > eval(mat.at(0).GetNrows());
    TMatrixD init = mat.at(0);
    init.Invert();

    for(UInt_t x = 0; x < mat.size(); x++){
        TMatrixD tempmat = init * mat.at(x);
        TVectorD tempeval;
        tempmat.EigenVectors(tempeval);
        for(UInt_t y = 0; y < tempeval.GetNrows(); y++){
            eval.at(y).push_back(tempeval(y));
        }
    }
    return eval;

}

vector<vector<Double_t> > findEvalBoot(vector<vector<vector<Double_t> > > corr, vector<vector<Double_t> > &err){

    Double_t N = corr.at(0).size(), NB = 5*N;
    vector<vector<vector<Double_t> > > evals;
    for(Int_t x = 0; x < NB; x++){
        vector<vector<vector<Double_t> > > tempcorr(corr.size());
        for(Int_t z = 0; z < N; z++){
            Int_t num = randGen->Integer(N);
            for(UInt_t w = 0; w < corr.size(); w++){
                tempcorr.at(w).push_back(corr.at(w).at(num));
            }
        }
        vector<vector<Double_t> > tempAv = average(tempcorr);
        tempcorr.clear();
        vector<TMatrixD> tempmat = constructMatrices(tempAv);
        tempAv.clear();
        evals.push_back(findEigenvalues(tempmat));
        tempmat.clear();
    }

    vector<vector<Double_t> > av(evals.at(0).size(),vector<Double_t>(evals.at(0).at(0).size()));
    for(UInt_t x = 0; x < evals.at(0).size(); x++){
        for(UInt_t y = 0; y < evals.at(0).at(0).size(); y++){
            for(UInt_t z = 0; z < NB; z++){
                av.at(x).at(y) += evals.at(z).at(x).at(y);
            }
            av.at(x).at(y) /= NB;
        }
    }

    vector<vector<Double_t> > error(evals.at(0).size(),vector<Double_t>(evals.at(0).at(0).size()));
    for(UInt_t x = 0; x < evals.at(0).size(); x++){
        for(UInt_t y = 0; y < evals.at(0).at(0).size(); y++){
            for(UInt_t z = 0; z < NB; z++){
                error.at(x).at(y) += (evals.at(z).at(x).at(y)-av.at(x).at(y))*(evals.at(z).at(x).at(y)-av.at(x).at(y));
            }
            error.at(x).at(y) = Sqrt(1./NB*error.at(x).at(y));
        }
    }
    err = error;
    return av;

}

vector<vector<Double_t> > findEffMassesBoot(vector<vector<vector<Double_t> > > corr, vector<vector<Double_t> > &err){

    Double_t N = corr.at(0).size(), NB = 5*N;
    vector<vector<vector<Double_t> > > allMasses;
    for(Int_t x = 0; x < NB; x++){
        vector<vector<vector<Double_t> > > tempcorr(corr.size());
        for(Int_t z = 0; z < N; z++){
            Int_t num = randGen->Integer(N);
            for(UInt_t w = 0; w < corr.size(); w++){
                tempcorr.at(w).push_back(corr.at(w).at(num));
            }
        }
        vector<vector<Double_t> > tempAv = average(tempcorr);
        tempcorr.clear();
        vector<TMatrixD> tempmat = constructMatrices(tempAv);
        tempAv.clear();
        vector<vector<Double_t> > tempeval = findEigenvalues(tempmat);
        tempmat.clear();

        vector<vector<Double_t> > tempMass(tempeval.size());
        for(UInt_t y = 0; y < tempMass.size(); y++){
            for(UInt_t z = 0; z < tempeval.at(0).size()-1; z++){
                massFunc.FixParameter(2, tempeval.at(0).size());
                massFunc.FixParameter(1, z);
                massFunc.FixParameter(0, tempeval.at(y).at(z)/tempeval.at(y).at(z+1));
                brf.Solve();
                tempMass.at(y).push_back(brf.Root());
            }
        }
        tempeval.clear();
        allMasses.push_back(tempMass);
    }

    vector<vector<Double_t> > av(allMasses.at(0).size(),vector<Double_t>(allMasses.at(0).at(0).size()));
    for(UInt_t x = 0; x < allMasses.at(0).size(); x++){
        for(UInt_t y = 0; y < allMasses.at(0).at(0).size(); y++){
            for(UInt_t z = 0; z < NB; z++){
                av.at(x).at(y) += allMasses.at(z).at(x).at(y);
            }
            av.at(x).at(y) /= NB;
        }
    }

    vector<vector<Double_t> > error(allMasses.at(0).size(),vector<Double_t>(allMasses.at(0).at(0).size()));
    for(UInt_t x = 0; x < allMasses.at(0).size(); x++){
        for(UInt_t y = 0; y < allMasses.at(0).at(0).size(); y++){
            for(UInt_t z = 0; z < NB; z++){
                error.at(x).at(y) += (allMasses.at(z).at(x).at(y)-av.at(x).at(y))*(allMasses.at(z).at(x).at(y)-av.at(x).at(y));
            }
            error.at(x).at(y) = Sqrt(1./NB*error.at(x).at(y));
        }
    }
    err = error;
    return av;

}

vector<TGraphErrors*> makeGraphs(vector<vector<Double_t> > evals, vector<vector<Double_t> > errors, TString pos){

    vector<TGraphErrors*> graphs;

    for(UInt_t x = 0; x < evals.size(); x++){
        if(pos == "half") graphs.push_back(new TGraphErrors(evals.at(0).size(), halfIntArr, &(evals.at(x))[0], zeroArr, &(errors.at(x))[0]));
        else graphs.push_back(new TGraphErrors(evals.at(0).size(), intArr, &(evals.at(x))[0], zeroArr, &(errors.at(x))[0]));
    }

    return graphs;

}

// TMatrixD weightMatrix(vector<vector<Double_t> > corr){

//     Int_t n = corr.at(0).size();
//     Double_t N = corr.size();
//     vector<Double_t> av = average(corr);
//     TMatrixD covMatrix(n, n);

//     for(Int_t x = 0; x < n; x++){
//         for(Int_t y = 0; y < n; y++){
//             for(Int_t z = 0; z < N; z++){
//                 covMatrix(x,y) += (corr.at(z).at(x)-av.at(x))*(corr.at(z).at(y)-av.at(y))*1e25;
//             }
//             covMatrix(x,y) = 1.0/((N-1.0)*N)*covMatrix(x,y);
//         }
//     }

//     TMatrixD inv = covMatrix.Invert();
//     inv *= 1e-25;

//     return inv;

// }

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
    
    if(name.Contains("Corr")) mg->SetMinimum(1e-18);
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
    can->SaveAs("./Plots/" + name + ".cpp");

    //can->Clear();
}
