//c++
#include <map>
#include <vector>
#include <cmath>
#include <cassert>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <iomanip>
#include <math.h>
#include <stdlib.h>
#include <stdexcept>
#include <algorithm>
#include "stdio.h"

//ROOT
#include "TSystem.h"
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TMath.h"
#include "TMinuit.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TH1D.h"
#include "TF1.h"
#include "TPaveText.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TLorentzVector.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TVector3.h"
#include "TVectorD.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TDirectoryFile.h"
#include "TSystemDirectory.h"
#include "TString.h"
#include "TArc.h"
#include "TText.h"
#include "TMarker.h"
#include "TObjArray.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TGraphAsymmErrors.h"
#include "TRandom3.h"
#include "TGaxis.h"
#include "TPad.h"
#include "Math/DistFunc.h"
#include <cassert>

using namespace std;


Double_t fitcos(Double_t *x,Double_t *par) {

  Double_t fitval = par[0]*cos((x[0]-(2*par[1]))) + par[2];
  return fitval;
}

TVector3 calculateImpactParamVec(TLorentzVector pion, TVector3 primaryVertex, TVector3 tauDecayVertex) {
    TVector3 pionVector = pion.Vect();
    //cout << "pion Px " << pionVector[0] << " Py " << pionVector[1] << " Pz " << pionVector[2] << endl;
    TVector3 pv(primaryVertex.x(), primaryVertex.y(), primaryVertex.z());
    //cout << "PV x " << pv[0] << " y " << pv[1] << " z " << pv[2] << endl;
    TVector3 initialPoint(tauDecayVertex.x(), tauDecayVertex.y(), tauDecayVertex.z());
    //cout << "initialPoint x " << initialPoint[0] << " y " << initialPoint[1] << " z " << initialPoint[2] << endl;
    // solve equations in order to determine shortest distance from the PV to the extrapolated pi-track
    double k = (pv - initialPoint) * pionVector / pionVector.Mag2();
    // calculate impact parameter vector
    TVector3 IPVector = initialPoint + k * pionVector - pv;
    //cout << "IPVector x " << IPVector[0] << " y " << IPVector[1] << " z " << IPVector[2] << endl;
    return IPVector;
}

float Acoplanarity_IP(TLorentzVector piPlus, TLorentzVector ipVectorPlus, TLorentzVector piMinus, TLorentzVector ipVectorMinus,
                            const TLorentzVector& referenceFrame) {
    TVector3 boostIntoReferenceFrame = (-1) * referenceFrame.BoostVector();
    piPlus.Boost(boostIntoReferenceFrame);
    piMinus.Boost(boostIntoReferenceFrame);
    ipVectorPlus.Boost(boostIntoReferenceFrame);
    ipVectorMinus.Boost(boostIntoReferenceFrame);
    TVector3 ip3VectorPlus = ipVectorPlus.Vect();
    TVector3 ip3VectorMinus = ipVectorMinus.Vect();
    // parallel projections of impact vector onto pion direction
    TVector3 ipVectorParallelPlus = (ip3VectorPlus * piPlus.Vect().Unit()) * piPlus.Vect().Unit();
    TVector3 ipVectorParallelMinus = (ip3VectorMinus * piMinus.Vect().Unit()) * piMinus.Vect().Unit();
    // perpendicular component of impact vector
    TVector3 ipVectorPerpPlus = (ip3VectorPlus - ipVectorParallelPlus).Unit();
    TVector3 ipVectorPerpMinus = (ip3VectorMinus - ipVectorParallelMinus).Unit();
    // Triple-odd correlation
    float triplecorr = piMinus.Vect().Unit().Dot(ipVectorPerpPlus.Cross(ipVectorPerpMinus));

    float phistar = TMath::ACos(ipVectorPerpPlus * ipVectorPerpMinus);
    if (triplecorr < 0) { phistar = TMath::TwoPi() - phistar; }

    return phistar;
}


void phistarneu() {


  gROOT->Reset();
  //gStyle->SetOptStat(0);
  gStyle->SetOptStat(1111);
  gStyle->SetOptFit(1111);
//  gStyle->SetPalette(1);
  gStyle->SetPalette(55);
  gStyle->SetStatFontSize(0.03);
  gStyle->SetStatY(0.90);                // Set y-position (fraction of pad size)
  gStyle->SetStatX(0.90);                // Set x-position (fraction of pad size)
  gStyle->SetStatH(0.20);                // Set height of stat-box (fraction of pad size)
  gStyle->SetStatW(0.10);                // Set width of stat-box (fraction of pad size)

  std::ofstream out("quant_results/phistarneu.txt");  ///Name of output file with info

  auto pi = TMath::Pi();
  //int bins = 12;
  int bins = 7;

  TH1F* hphistar = new TH1F("hphistar", "#varphi* IP-IP Analysis", bins, 0.5, 5.9);
  hphistar->SetTitle(""); hphistar->SetXTitle("#varphi (rad)");


  TTree *treeU = new TTree () ;
  treeU->ReadFile("NeuNeuTauvert.tex");
  Int_t evtNum;
  Float_t px_plus, py_plus, pz_plus, e_plus;
  Float_t px_minus, py_minus, pz_minus, e_minus;
  Float_t x_tau_plus, y_tau_plus, z_tau_plus, t_tau_plus;
  Float_t x_tau_minus, y_tau_minus, z_tau_minus, t_tau_minus;
  Float_t Weight;
  treeU->SetBranchAddress("evtNum",&evtNum);
  treeU->SetBranchAddress("px_plus",&px_plus);
  treeU->SetBranchAddress("py_plus",&py_plus);
  treeU->SetBranchAddress("pz_plus",&pz_plus);
  treeU->SetBranchAddress("e_plus",&e_plus);
  treeU->SetBranchAddress("px_minus",&px_minus);
  treeU->SetBranchAddress("py_minus",&py_minus);
  treeU->SetBranchAddress("pz_minus",&pz_minus);
  treeU->SetBranchAddress("e_minus",&e_minus);
  treeU->SetBranchAddress("x_tau_plus",&x_tau_plus);
  treeU->SetBranchAddress("y_tau_plus",&y_tau_plus);
  treeU->SetBranchAddress("z_tau_plus",&z_tau_plus);
  treeU->SetBranchAddress("t_tau_plus",&t_tau_plus);
  treeU->SetBranchAddress("x_tau_minus",&x_tau_minus);
  treeU->SetBranchAddress("y_tau_minus",&y_tau_minus);
  treeU->SetBranchAddress("z_tau_minus",&z_tau_minus);
  treeU->SetBranchAddress("t_tau_minus",&t_tau_minus);
  treeU->SetBranchAddress("Weight",&Weight);
  Long64_t nentriesU = treeU->GetEntries();

  for (Long64_t i = 0; i < nentriesU; i++) {
    treeU->GetEntry(i);
    TLorentzVector piPlus, piMinus;
    piPlus.SetPxPyPzE(px_plus,py_plus,pz_plus,e_plus);
    piMinus.SetPxPyPzE(px_minus,py_minus,pz_minus,e_minus);
    float mass_plus = piPlus.M();
    float mass_minus = piMinus.M();
    out << "px_plus  " << px_plus  << " py_plus  " << py_plus  << " pz_plus  " << pz_plus  << " e_plus  " << e_plus  << " mass_plus  " << mass_plus << endl;
    out << "px_minus " << px_minus << " py_minus " << py_minus << " pz_minus " << pz_minus << " e_minus " << e_minus << " mass_minus " << mass_minus << endl;
    out << "x_tau_plus    " << x_tau_plus    << " y_tau_plus    " << y_tau_plus    << " z_tau_plus    " << z_tau_plus    << " t_tau_plus   " << t_tau_plus << endl;
    out << "x_tau_minus   " << x_tau_minus    << " y_tau_minus    " << y_tau_minus    << " z_tau_minus    " << z_tau_minus    << " t_tau_minus   " << t_tau_minus << endl;
    out << "weight" << Weight << endl;
    TVector3 primaryVertex(0.,0.,0.);
    TVector3 decayVertexPositive(x_tau_plus,y_tau_plus,z_tau_plus);
    TVector3 decayVertexNegative(x_tau_minus,y_tau_minus,z_tau_minus);
    TLorentzVector referenceFrame = piPlus + piMinus;
    TLorentzVector ipVectorPlus(calculateImpactParamVec(piPlus, primaryVertex, decayVertexPositive).Unit(), 0.);
    TLorentzVector ipVectorMinus(calculateImpactParamVec(piMinus, primaryVertex, decayVertexNegative).Unit(), 0.);
    float phistar_mc = Acoplanarity_IP(piPlus, ipVectorPlus, piMinus, ipVectorMinus, referenceFrame);
    hphistar->Fill(phistar_mc,Weight);
    float phistasr_deg = phistar_mc*TMath::RadToDeg();
    out << "phistar_mc " << phistar_mc << " in degrees " << phistasr_deg << endl;
  }

  hphistar->Sumw2();
  hphistar->SetMarkerStyle(20);

  hphistar->GetXaxis()->SetTitle("#varphi* (rad)");
  hphistar->GetYaxis()->SetTitle("Normalised Events");
  hphistar->Scale(1. / hphistar->Integral());

  TCanvas* c1 = new TCanvas("c1", "Phi Star IP-IP Analysis", 850, 550);
  c1->cd();

  hphistar->Draw();

  TF1 *f1 = new TF1 ("f1", "[0]*cos((x-(2*[1]))) + [2]",0,2.*pi);

  //f1->SetParNames ("u","#phi_{#tau}","w","bw");
  //f1->FixParameter(3,Plot->GetBinWidth(1));
  f1->SetParNames ("u","#phi_{#tau}","w");
  f1->SetParLimits(1,0.0,pi);

  f1->SetParameters(-0.1,0.,0.15);      // even
  //f1->SetParameters(-0.1,pi/2.,0.15);   // odd
  //f1->SetParameters(-0.1,pi/4.,0.15);   // mix

  hphistar->Fit("f1","R");
  hphistar->Fit("f1","R");

  auto legend = new TLegend(0.1,0.8,0.31,0.9);
  legend->AddEntry(f1,"u cos(#varphi*-2#phi_{#tau})+w","l");
  legend->AddEntry(hphistar,"#varphi*","lep");
  legend->SetTextSize(0.037);
  legend->Draw();

  c1->SaveAs("plots/neu_neu_phistar_mc.png");

} // phistar



