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
#include <cmath>
#include <fstream>
#include <TString.h>
#include "TMath.h"
#include <TPaveStats.h>
#include <TText.h>
#include <TLatex.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TGraphErrors.h>

#endif

using namespace TMath;
using namespace std;

// Declare functions
TGraph* readData(TString);
void graph(TGraph *graph, const TString graphName, TCanvas *can, const TString xTitle, const TString yTitle, const TString name);
void graph(TGraph *graph1, TGraph *graph2, const TString graphName, TCanvas *can, const TString xTitle, const TString yTitle, const TString name);

TGraphErrors *graphLarge, *graphSmall, *allPeaks;

/*
 * MAIN FUNCTION
 */

void compare(){
    
    cout << "\n\nStarting process...\n\n";

    TGraph *data = readData("chromaData");

    TCanvas *can = new TCanvas("Canvas", "Canvas", 1600, 900);
    graph(data, "Chroma data", can, "nth measurement", "Average plaquette (|Re(Tr(U_{#mu#nu}))|)", "magField");

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

void graph(TGraph *graph, const TString graphName, TCanvas *can, const TString xTitle, const TString yTitle, const TString name){

    gStyle->SetFitFormat("4.3f");

    //gStyle->SetOptStat(2210);
    gStyle->SetOptFit(1111);
    gStyle->SetOptStat(kFALSE);
    //gStyle->SetOptFit(1100);

    if(!graph){
        cout << "Error: Graph \"" << graphName << "\" not defined" << endl;
        return;
    }

    graph->SetLineWidth(2);
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(1);
    graph->SetMarkerColor(kRed);
    graph->SetLineColor(kBlack);
    // TF1 *f = new TF1(graphName, "[0]+[1]*TMath::Power(TMath::BesselJ1([2]*x-[3])*TMath::BesselJ1([2]*x-[3])/(([2]*x-[3])*([2]*x-[3])),0.5)");
    // TF1 *f = new TF1(graphName, "[0]+[1]/TMath::Log([4]*TMath::BesselJ1([2]*x-[3])*TMath::BesselJ1([2]*x-[3])/(([2]*x-[3])*([2]*x-[3])) + [5])");
    // f->SetParameter(0, 250);
    // f->SetParameter(1, 100);
    // f->SetParameter(2, 0.395);
    // f->SetParameter(3, 1);
    // f->SetParameter(4, 10000);
    // f->SetParameter(5, 100);
    // f->SetParName(0, "Offset");
    // f->SetParName(1, "Peak intensity");
    // f->SetParName(2, "Aperture");
    // f->SetParName(3, "Center value");
    // f->SetParName(4, "Amp");
    // f->SetParLimits(0, 100, 300);
    // f->SetParLimits(1, 10, 100);
    // f->SetParLimits(2, 0.1, 1);
    // f->SetParLimits(3, -0.2, 2);
    // f->SetParLimits(4, 50, 10000);
    // f->SetParLimits(5, 1, 10000);
    // f->SetLineColor(kBlue);
    // f->SetLineWidth(4);
    //graph->Fit(f, "ME");

    TF1 *f1 = graph->GetFunction("fitRight");
    // TF1 *f2 = graph->GetFunction("linFit2");

    gStyle->SetLegendBorderSize(0);
    TLegend *leg = new TLegend(0.655,0.675,0.885,0.875);
    leg->SetTextSize(0.055);
    leg->AddEntry(graph, "Data","p");
    // leg->AddEntry(f1, "Fit","l");
    // leg->AddEntry(f2, "Anti-Stokes linear fit","l");
    
    graph->Draw("ap");
    leg->Draw("same");

    // add axis labels
    graph->GetXaxis()->SetTitle(xTitle);
    graph->GetXaxis()->CenterTitle();
    graph->GetXaxis()->SetTitleSize(0.055);
    graph->GetXaxis()->SetTitleOffset(0.82);
    //graph->GetXaxis()->SetLabelOffset(0.010);
    graph->GetXaxis()->SetLabelSize(0.05);
    graph->GetYaxis()->SetTitle(yTitle);
    graph->GetYaxis()->CenterTitle();
    graph->GetYaxis()->SetTitleSize(0.055);
    graph->GetYaxis()->SetTitleOffset(0.9);
    graph->GetYaxis()->SetLabelSize(0.05);
    gStyle->SetTitleSize(0.08, "t");
    graph->SetTitle(graphName);
    can->Update();

    can->SaveAs(name + ".png");

    //can->Clear();
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
    graph1->SetMarkerColor(kRed);
    graph1->SetLineColor(kBlack);
    graph2->SetLineWidth(2);
    graph2->SetMarkerStyle(20);
    graph2->SetMarkerSize(1.5);
    graph2->SetMarkerColor(kRed);
    graph2->SetLineColor(kBlack);

    TF1 *f1 = graph1->GetFunction("fitSmall");
    TF1 *f2 = graph2->GetFunction("fitLarge");

    f1->SetLineColor(kBlue);
    f1->SetLineWidth(4);
    f2->SetLineColor(kGreen+1);
    f2->SetLineWidth(4);
    f2->SetLineStyle(9);

    gStyle->SetLegendBorderSize(0);
    TLegend *leg = new TLegend(0.655,0.675,0.885,0.875);
    leg->SetTextSize(0.055);
    leg->AddEntry(graph1, "Data","lep");
    leg->AddEntry(f1, "Small Stokes fit","l");
    leg->AddEntry(f2, "Large Stokes fit","l");

    TMultiGraph *mg = new TMultiGraph();
    mg->SetTitle("Nitrogen");

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