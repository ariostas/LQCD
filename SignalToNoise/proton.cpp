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
vector<vector<vector<Double_t> > > readData(vector<TString>);
vector<vector<Double_t> > average(vector<vector<vector<Double_t> > >);
vector<vector<Double_t> > averageBoot(vector<vector<vector<Double_t> > >);
vector<Double_t> average(vector<vector<Double_t> >);
vector<Double_t> averageBoot(vector<vector<Double_t> >);
vector<vector<Double_t> > errorBoot(vector<vector<vector<Double_t> > >);
vector<TMatrixD> constructMatrices(vector<vector<Double_t> >);
vector<vector<Double_t> > findEigenvalues(vector<TMatrixD>&, vector<vector<TMatrixD> >&);
vector<vector<Double_t> > findEvalBoot(vector<vector<vector<Double_t> > >, vector<vector<Double_t> >&);
vector<vector<Double_t> > findEffMassesBoot(vector<vector<vector<Double_t> > >, vector<vector<Double_t> >&);
vector<TGraphErrors*> makeGraphs(vector<vector<Double_t> >, vector<vector<Double_t> >, TString pos = "int");
vector<vector<Double_t> > effMass(vector<vector<Double_t> >);
vector<Double_t> effMass(vector<Double_t>);
TMatrixD findSink(vector<vector<vector<Double_t> > >, TMatrixD, vector<TMatrixD>, Int_t);
TMatrixD weightMatrix(vector<vector<Double_t> >);
Double_t findMass(vector<vector<Double_t> >, Int_t, Int_t, Double_t, Double_t, Double_t *chi2 = 0);
Double_t findMassBoot(vector<vector<Double_t> >, Int_t, Int_t, Double_t&, Double_t, Double_t, Double_t *chi2 = 0, Double_t *chi2err = 0);
TMatrixD transpose(TMatrixD);
vector<Double_t> findCorrFromVec(TMatrixD, TMatrixD, vector<TMatrixD>);
vector<vector<Double_t> > findGndCorr(vector<vector<vector<Double_t> > >);
vector<vector<Double_t> > findGndMasses(vector<vector<vector<Double_t> > >, vector<vector<Double_t> >&);
void findOptimalSN(vector<vector<vector<Double_t> > >, TMatrixD&, TMatrixD&, vector<TMatrixD>, Int_t);
void scanFitRange(vector<vector<Double_t> >);
void graph(vector<TGraphErrors*>, vector<TString>, const TString, const TString, const TString, vector<vector<Double_t> > *mass = 0, TString flags = "");

// Declare common objects
TRandom3 *randGen = new TRandom3();
Double_t intArr[500], halfIntArr[500], zeroArr[500];
// vector<TString> corrTypes = {"P_1.P_1.PP", "DG2_1.P_1.SP", "DG4_1.P_1.SP", "P_1.DG2_1.PS", "DG2_1.DG2_1.SS", "DG4_1.DG2_1.SS",
//                              "P_1.DG4_1.PS", "DG2_1.DG4_1.SS", "DG4_1.DG4_1.SS"};

// Main function
void proton(){

    TH1::StatOverflows(kTRUE);
    
    cout << "\n\nStarting process...\n\n";

    for(Int_t x = 0; x < 500; x++){
        intArr[x] = x; halfIntArr[x] = x+0.5; zeroArr[x] = 0;
    }

    vector<TString> corrTypes;
    for(Int_t x = 1; x <= 3; x++){
         for(Int_t y = 1; y <= 3; y++){
            // corrTypes.push_back(TString::Format("proton_SrcDG%i_SnkDG%i_Interp4.dat", x, y));
            corrTypes.push_back(TString::Format("rho_x.D-2449.DG%i_1.DG%i_1.SS", x, y));
         }
    }

    vector<vector<vector<Double_t> > > correlators092 = readData(corrTypes);

    // vector<vector<Double_t> > gndcorr092 = findGndCorr(correlators092);
    // vector<vector<Double_t> > evalErrors092, evals092 = findEvalBoot(correlators092, evalErrors092);

    vector<vector<Double_t> > massErrors092, masses092 = findGndMasses(correlators092, massErrors092);

    // cout << masses092.at(2).size() << "  " << masses092.at(2).at(0) << endl;

    // vector<TGraphErrors*> graphs092 = makeGraphs(evals092, evalErrors092);

    vector<TGraphErrors*> massgraphs092 = makeGraphs(masses092, massErrors092);

    // vector<vector<Double_t> > masses(1);
    // Double_t m1, e1;
    // m1 = findMassBoot(gndcorr092, 4, 14, e1, 7.33767e-01, 5.02458e-01);
    // masses.at(0).push_back(m1); masses.at(0).push_back(e1); masses.at(0).push_back(4); masses.at(0).push_back(14);

    vector<TString> names = {"Smeared source and sink", "Variational Method", "Source optimized + S/N sink", "S/N optimized source and sink"};

    // graph(graphs092, names, "#Deltat", "Eigenvalues(#Deltat)", "protonVarCorrs092", 0, "log");

    graph(massgraphs092, names, "#Deltat", "M(#Deltat)", "protonVarMass092");

    // scanFitRange(gndcorr092);

}

// Read correlator data from output of strip_hadspec
vector<vector<vector<Double_t> > > readData(vector<TString> files){

    vector<vector<vector<Double_t> > > corrs;

    for(UInt_t x = 0; x < files.size(); x++){
        cout << "Reading file " << files.at(x) << "..." << endl;

        TString filename = "/home/arios/Documents/LQCDConfigs/cl3_16_48_b6p1_m0p2450/hadspec/" + files.at(x);

        ifstream ifs(filename); if(!ifs.is_open()){cout << "Error. File " << filename << " not found. Exiting...\n"; assert(0);}

        Int_t NConfigs = 0, Nt = 0, IsComplex = 0, Ns = 0, temp = 0;
        Double_t temp2 = 0;
        ifs >> NConfigs >> Nt >> IsComplex >> Ns >> temp >> temp2;

        vector<vector<Double_t> > corr(NConfigs);

        if(NConfigs == 0){cout << "Error reading the file. Exiting..." << endl; assert(0);}

        Int_t deltaT = 0;
        Double_t re = 0, im = 0;//, unkn = 0;
    
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

vector<vector<Double_t> > findEigenvalues(vector<TMatrixD> &mat, vector<vector<TMatrixD> > &evec){

    vector<vector<Double_t> > eval(mat.at(0).GetNrows());
    TMatrixD init = mat.at(0);
    TMatrixD init2 = mat.at(3);
    init.Invert();
    evec.clear();
    evec = vector<vector<TMatrixD> >(mat.size(),vector<TMatrixD>(0));

    for(UInt_t x = 0; x < mat.size(); x++){
        TMatrixD tempmat = init * mat.at(x);
        TMatrixD tempmat2 = init2 * mat.at(x);
        TVectorD tempeval, tempeval2;
        tempmat.EigenVectors(tempeval);
        const TMatrixD mvec = tempmat.EigenVectors(tempeval2);
        for(UInt_t y = 0; y < tempeval.GetNrows(); y++){
            eval.at(y).push_back(tempeval(y));
        }
        
        for(UInt_t y = 0; y < mvec.GetNcols(); y++){
        	TMatrixD vec(mvec.GetNrows(),1);
        	for(UInt_t z = 0; z < mvec.GetNrows(); z++){
        		vec(z,0) = mvec(z,y);
        	}
        	evec.at(y).push_back(vec);
        }

    }
    return eval;
}

vector<vector<Double_t> > findEvalBoot(vector<vector<vector<Double_t> > > corr, vector<vector<Double_t> > &err){

    Double_t N = corr.at(0).size(), NB = N;
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
        vector<vector<TMatrixD> > eigenvectors;
        evals.push_back(findEigenvalues(tempmat, eigenvectors));
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

// Find the ground state mass using different methods
// (Smeared operators, variational method, variational + S/N, S/N)
vector<vector<Double_t> > findGndMasses(vector<vector<vector<Double_t> > > corr, vector<vector<Double_t> > &err){

    Double_t N = corr.at(0).size(), NB = 100;
    vector<vector<vector<Double_t> > > masses(NB);
    for(Int_t x = 0; x < NB; x++){
        cout << x << endl;
        vector<vector<vector<Double_t> > > tempcorr(corr.size());
        for(Int_t z = 0; z < N; z++){
            Int_t num = randGen->Integer(N);
            for(UInt_t w = 0; w < corr.size(); w++){
                tempcorr.at(w).push_back(corr.at(w).at(num));
            }
        }
        vector<vector<Double_t> > tempAv = average(tempcorr);

        // First find the mass simply using smeared operators
        masses.at(x).push_back(effMass(tempAv.at(4)));

        // Then find mass using the variational method
        vector<TMatrixD> tempmat = constructMatrices(tempAv);
        vector<vector<TMatrixD> > eigenvectors;
        masses.at(x).push_back(effMass(findEigenvalues(tempmat, eigenvectors).at(0)));

        // Now find the mass using a S/N optimized sinks
        Int_t ts = 20;
        TMatrixD source = eigenvectors.at(0).at(ts);
        source *= 1./Sqrt((transpose(source) * source)(0,0));
        TMatrixD sink = findSink(tempcorr, source, tempmat, ts);
        vector<Double_t> tempcorr2 = findCorrFromVec(source, sink, tempmat);
        // for(UInt_t a = 0; a < tempcorr2.size(); a++){
        // 	cout << tempcorr2.at(a) << endl;
        // }
        // cout << endl;
        // cout << (transpose(source) * source)(0,0) << "  " << (transpose(sink) * sink)(0,0) << "  " << (transpose(sink) * source)(0,0) << endl;
        masses.at(x).push_back(effMass(tempcorr2));

        // Finally, find the mass using S/N optimized sources and sinks
        // TMatrixD start = source;
        TMatrixD start(5,1); start(1,0) = 1;
        TMatrixD sourceSN = start, sinkSN = start;
        findOptimalSN(tempcorr, sourceSN, sinkSN, tempmat, ts);
        vector<Double_t> tempcorr3 = findCorrFromVec(sourceSN, sinkSN, tempmat);
        masses.at(x).push_back(effMass(tempcorr3));

        cout << (transpose(sourceSN) * start)(0,0) << "  " << (transpose(sinkSN) * start)(0,0) << "  " << (transpose(sinkSN) * sourceSN)(0,0) << endl;

        // Delete everything
        eigenvectors.clear();
        tempcorr.clear();
        tempcorr2.clear();
        tempcorr3.clear();
        tempAv.clear();
        tempmat.clear();
    }

    vector<vector<Double_t> > av(masses.at(0).size(),vector<Double_t>(masses.at(0).at(0).size()));
    for(UInt_t x = 0; x < masses.at(0).size(); x++){
        for(UInt_t y = 0; y < masses.at(0).at(0).size(); y++){
            Double_t bad = 0;
            for(UInt_t z = 0; z < NB; z++){
                if(masses.at(z).at(x).at(y) != masses.at(z).at(x).at(y)) bad++;
                else av.at(x).at(y) += masses.at(z).at(x).at(y);
            }
            av.at(x).at(y) /= (NB-bad);
        }
    }

    vector<vector<Double_t> > error(masses.at(0).size(),vector<Double_t>(masses.at(0).at(0).size()));
    for(UInt_t x = 0; x < masses.at(0).size(); x++){
        for(UInt_t y = 0; y < masses.at(0).at(0).size(); y++){
            Double_t bad = 0;
            for(UInt_t z = 0; z < NB; z++){
                if(masses.at(z).at(x).at(y) != masses.at(z).at(x).at(y)) bad++;
                else error.at(x).at(y) += (masses.at(z).at(x).at(y)-av.at(x).at(y))*(masses.at(z).at(x).at(y)-av.at(x).at(y));
            }
            error.at(x).at(y) = Sqrt(1./(NB-bad)*error.at(x).at(y));
        }
    }
    err = error;
    return av;

}

vector<Double_t> effMass(vector<Double_t> corr){
    vector<Double_t> masses;
    for(UInt_t x = 0; x < corr.size();x++){
        Double_t m = Log(corr.at(x)/corr.at((x+1)%corr.size()));
        if(x > Double_t(corr.size())/2.) m *= -1.;
        masses.push_back(m);
    }
    return masses;
}

TMatrixD findSink(vector<vector<vector<Double_t> > > corr, TMatrixD evec, vector<TMatrixD> C, Int_t ts){

	Int_t matSize = evec.GetNrows();
	TMatrixD sink(matSize, 1);

	TMatrixD sigma2(matSize, matSize);
	for(UInt_t x = 0; x < corr.at(0).size(); x++){

		TMatrixD curlyC(matSize, matSize);
		for(UInt_t y = 0; y < corr.size(); y++){
			curlyC(y%matSize,y/matSize) = corr.at(y).at(x).at(ts);
		}

		sigma2 = sigma2 + (curlyC * evec * transpose(evec) * transpose(curlyC));

	}

	sigma2 *= 1./Double_t(corr.at(0).size());
	sigma2 *= 1e35;
	sigma2.Invert();
	sigma2 *= 1e35;

	Double_t A = (transpose(evec) * transpose(C.at(ts)) * sigma2 * sigma2 * C.at(ts) * evec)(0,0);
	A = 1./Sqrt(Abs(A));

	sink = A * sigma2 * C.at(ts) * evec;

	return sink;

}

void findOptimalSN(vector<vector<vector<Double_t> > > corr, TMatrixD &source, TMatrixD &sink, vector<TMatrixD> C, Int_t ts){

	Int_t matSize = source.GetNrows();

	for(Int_t iteration = 0; iteration < 50; iteration++){
		TMatrixD sigma2source(matSize, matSize), sigma2sink(matSize, matSize);
		for(UInt_t x = 0; x < corr.at(0).size(); x++){

			TMatrixD curlyC(matSize, matSize);
			for(UInt_t y = 0; y < corr.size(); y++){
				curlyC(y%matSize,y/matSize) = corr.at(y).at(x).at(ts);
			}

			sigma2source = sigma2source + (curlyC * source * transpose(source) * transpose(curlyC));
			sigma2sink = sigma2sink + (transpose(curlyC) * sink * transpose(sink) * curlyC);

		}

		sigma2source *= 1./Double_t(corr.at(0).size());
		sigma2sink *= 1./Double_t(corr.at(0).size());
		sigma2source *= 1e35;
		sigma2sink *= 1e35;
		sigma2source.Invert();
		sigma2sink.Invert();
		sigma2source *= 1e35;
		sigma2sink *= 1e35;

		Double_t Asources = (transpose(source) * transpose(C.at(ts)) * sigma2source * sigma2source * C.at(ts) * source)(0,0);
		Double_t Asinks = (transpose(sink) * C.at(ts) * sigma2sink * sigma2sink * transpose(C.at(ts)) * sink)(0,0);
		Asources = 1./Sqrt(Abs(Asources));
		Asinks = 1./Sqrt(Abs(Asinks));

		TMatrixD oldsource = source, oldsink = sink;

		source = Asinks * sigma2sink * transpose(C.at(ts)) * oldsink;
		sink = Asources * sigma2source * C.at(ts) * oldsource;

        source *= 1./Sqrt((transpose(source) * source)(0,0));
        sink *= 1./Sqrt((transpose(sink) * sink)(0,0));

		cout << iteration << "  " << (transpose(oldsource) * source)(0,0) << "  " << (transpose(oldsink) * sink)(0,0) << endl;
	}
}

TMatrixD transpose(TMatrixD mat){
	TMatrixD tr(mat.GetNcols(), mat.GetNrows());
	for(Int_t x = 0; x < mat.GetNcols(); x++){
		for(Int_t y = 0; y < mat.GetNrows(); y++){
			tr(x,y) = mat(y,x);
		}
	}
	return tr;
}

vector<Double_t> findCorrFromVec(TMatrixD source, TMatrixD sink, vector<TMatrixD> C){

	vector<Double_t> corr;

	for(UInt_t x = 0; x < C.size(); x++){
		corr.push_back((transpose(sink) * C.at(x) * source)(0,0));
	}
	return corr;

}

vector<vector<Double_t> > findGndCorr(vector<vector<vector<Double_t> > > corr){

    Double_t N = corr.at(0).size();
    vector<vector<Double_t> > gndcorr;
    for(Int_t z = 0; z < N; z++){
        vector<vector<Double_t> > tempcorr(corr.size());
        for(UInt_t w = 0; w < corr.size(); w++){
            tempcorr.at(w) = corr.at(w).at(z);
        }

        vector<TMatrixD> tempmat = constructMatrices(tempcorr);
        tempcorr.clear();
        vector<vector<TMatrixD> > eigenvectors;
        gndcorr.push_back(findEigenvalues(tempmat, eigenvectors).at(0));
        tempmat.clear();
    }

    return gndcorr;
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
        vector<vector<TMatrixD> > eigenvectors;
        vector<vector<Double_t> > tempeval = findEigenvalues(tempmat, eigenvectors);
        tempmat.clear();

        vector<vector<Double_t> > tempMass(tempeval.size());
        for(UInt_t y = 0; y < tempMass.size(); y++){
            for(UInt_t z = 0; z < tempeval.at(0).size()-1; z++){
                Double_t m = Log(tempeval.at(y).at(z)/tempeval.at(y).at(z+1));
                if(z > Double_t(tempeval.at(y).size())/2.) m *= -1.;
                tempMass.at(y).push_back(m);

                // massFunc.FixParameter(2, tempeval.at(0).size());
                // massFunc.FixParameter(1, z);
                // massFunc.FixParameter(0, tempeval.at(y).at(z)/tempeval.at(y).at(z+1));
                // brf.Solve();
                // tempMass.at(y).push_back(brf.Root());
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
        }
    }

    TF2 *fitFunc = new TF2("Fit function", func, massGuess-0.5, massGuess+0.5, AGuess-0.5, AGuess+0.5);

    Double_t mass, A;

    fitFunc->GetMinimumXY(mass, A);

    if(chi2 != 0) (*chi2) = fitFunc->Eval(mass, A)/Double_t(nmax - nmin - 1);

    // cout << "chi2: " << fitFunc->Eval(mass, A) << endl;

    delete fitFunc;

    // massFunc.FixParameter(1, nmin);
    // massFunc.FixParameter(0, av.at(nmin)/av.at(nmin+1));
    // brf.Solve();

    // cout << mass << endl;

    return mass;
    
}

Double_t findMassBoot(vector<vector<Double_t> > corr, Int_t nmin, Int_t nmax, Double_t &err, Double_t massGuess, Double_t AGuess, Double_t *chi2, Double_t *chi2err){

    Double_t N = corr.size(), NB = Nint(N/20);
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
            covMatrix(x,y) = 1.0/((N-1.0)*N)*covMatrix(x,y)*1e20;
        }
    }

    TMatrixD inv = covMatrix.Invert();
    inv *= 1e20;

    return inv;

}

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

void scanFitRange(vector<vector<Double_t> > corr){

    Double_t arr[100];
    for(Int_t x = 1; x < 101; x++){
        arr[x-1] = x;
    }

    Int_t N = 12;
    vector<Double_t> m(N), me(N), c2(N), c2e(N);

    for(Int_t x = 0; x < N; x++){
        m.at(x) = findMassBoot(corr, 12-x, 14, me.at(x), 7.33767e-01, 5.11457e-05, &(c2.at(x)), &(c2e.at(x)));
    }

    TGraphErrors *massScan = new TGraphErrors(N, arr, &m[0], zeroArr, &me[0]);
    TGraphErrors *chiScan = new TGraphErrors(N, arr, &c2[0], zeroArr, &c2e[0]);

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

    if(graphs.size() > massLabels.size()){cout << "Error: Size of vectors does not match." << endl; return;}

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
            if(x < mass->size()){
                massFit.at(3*x)->SetLineColor(colors.at(x%6));
                massFit.at(3*x+1)->SetLineColor(colors.at(x%6));
                massFit.at(3*x+2)->SetLineColor(colors.at(x%6));
                massFit.at(3*x+1)->SetLineStyle(7);
                massFit.at(3*x+2)->SetLineStyle(7);
            }
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
    // f1->SetParameter(0, 0.7);
    // f1->SetParLimits(0, 0, 2);
    // f1->SetParameter(1, 1e-2);
    // f1->SetLineColor(kBlack);
    // graphs.at(0)->Fit(f1, "ME", "", 4, 14);
    // cout << f1->GetParameter(0) << endl;

    can->SaveAs(name + ".png");
    // can->SaveAs("./Plots/" + name + ".cpp");

    //can->Clear();
}
