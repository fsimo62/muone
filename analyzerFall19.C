
#include <vector>
#include <string>
#include <math.h>
#include <iostream>
#include <sstream>
#include <iomanip>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TString.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TEfficiency.h"

// Theory and experiment parameters
const double r_ = 0.999996345272132;
const double eMass_ = 0.510999;
const double muMass_ = 105.658;
const double beamEnergy_ = 150000.;

const double minTheta = 0.;

// Crystal parameters
const double innerSize_ = 25. + 0.5*2;
const double outerSize_ = 97.8 + 0.5*2;
// Position measurement parameters
const int showerSize_ = 1;
const double w0_ = 4.;

// const double mpvMuLandau_ = 350;

double eFromTheta(double theta)
{
    return eMass_*(1+pow(r_,2)*pow(TMath::Cos(theta),2)/(1-pow(r_,2)*pow(TMath::Cos(theta),2) ) );
}

void analyzerFall19(TString file = "/lustre/cmswork/mpresill/MUonE/Fall2019/depth23/ntuD23T20E130.root"
,   double maxTheta = 5.
,   int nEvents = -1
    ){

    TString suffix = file(file.Index("ntu")+3, file.Index(".root") - file.Index("ntu")-3);

    cout<<"--- You choose:"<<endl;
    cout<<"----- file = "<<file<<endl;
    cout<<"----- nEvents = "<<nEvents<<endl;
    cout<<"----- suffix = "<<suffix<<endl;

    gStyle->SetOptStat(0);
    gStyle->SetOptFit(1);

    TFile *f = new TFile(file);

    cout<<"--- File opened"<<endl;

    TTree *trTracks  = (TTree*)f->Get("tracks");
    TTree *trHits    = (TTree*)f->Get("hits");
    TTree *trCalo    = (TTree*)f->Get("calo");
    TTree *trConfig  = (TTree*)f->Get("config");

    cout<<"--- Tree read"<<endl;

    // Read config
    int cGeometry, cTrackers, cModules, cSubstations;
    double cCalZentry;
    trConfig->SetBranchAddress( "geometry", &cGeometry);
    trConfig->SetBranchAddress( "trackers", &cTrackers);
    trConfig->SetBranchAddress( "modules", &cModules);
    trConfig->SetBranchAddress( "substastions", &cSubstations);
    trConfig->SetBranchAddress( "calZentry", &cCalZentry);

    // Get n events, max tracks/hits/calo
    int eventID, nTracks, nHits;
    trTracks->SetBranchAddress("eventID", &eventID);
    trTracks->SetBranchAddress("nTracks", &nTracks);
    if(nEvents<0) nEvents = trTracks->GetEntries();
    if(nEvents>trTracks->GetEntries()) nEvents = trTracks->GetEntries();

    //this numebrs are hadcoded in the G4 simulation
    int maxTracks = 500;
    int maxHits = 1000;
    int maxParticles = 5000;
    int maxcrystals = 5000;

    // Matching variables

    // Tracks
    vector<int> t_PDG, t_parentPDG, t_trackID;
    vector<double> t_kinEnergy, t_xVertex, t_yVertex, t_zVertex, t_theta;

    t_PDG.resize(maxTracks, 0);
    t_parentPDG.resize(maxTracks, 0);
    t_trackID.resize(maxTracks, 0);
    t_kinEnergy.resize(maxTracks, 0);
    t_xVertex.resize(maxTracks, 0);
    t_yVertex.resize(maxTracks, 0);
    t_zVertex.resize(maxTracks, 0);
    t_theta.resize(maxTracks, 0);

    trTracks->SetBranchAddress("PDG", &(t_PDG[0]));
    trTracks->SetBranchAddress("parentPDG", &(t_parentPDG[0]));
    trTracks->SetBranchAddress("trackID", &(t_trackID[0]));
    trTracks->SetBranchAddress("kinEnergy", &(t_kinEnergy[0]));
    trTracks->SetBranchAddress("xVertex", &(t_xVertex[0]));
    trTracks->SetBranchAddress("yVertex", &(t_yVertex[0]));
    trTracks->SetBranchAddress("zVertex", &(t_zVertex[0]));
    trTracks->SetBranchAddress("theta", &(t_theta[0]));

    // Calo
    float c_totalEnergy;
    int c_ncrystals;
    vector<double> c_deposit, c_calposx, c_calposy;
    vector<int> c_crystalType;

    c_deposit.resize(maxcrystals, 0);
    c_calposx.resize(maxcrystals, 0);
    c_calposy.resize(maxcrystals, 0);
    c_crystalType.resize(maxcrystals, 0);

    trCalo->SetBranchAddress("totalEnergy", &c_totalEnergy);
    trCalo->SetBranchAddress("ncrystals", &c_ncrystals);
    trCalo->SetBranchAddress("crystalType", &(c_crystalType[0]));
    trCalo->SetBranchAddress("deposit", &(c_deposit[0]));
    trCalo->SetBranchAddress("calposx", &(c_calposx[0]));
    trCalo->SetBranchAddress("calposy", &(c_calposy[0]));

    cout<<"--- Branches setted"<<endl;

    trConfig->GetEvent(0);

    // Define result histograms
    TH1D *measuredAngle = new TH1D("measuredAngle", "#theta_{calo} "+suffix+";#theta", 50, minTheta, maxTheta);
    TH1D *angleResolution = new TH1D("angleResolution", "#sigma(#theta_{calo}) "+suffix+";#sigma(#theta_{calo})", 50, -(maxTheta-minTheta)/2, +(maxTheta-minTheta)/2);
    TH1D *measuredEnergy = new TH1D("measuredEnergy", "E_{loss} "+suffix+";E_{loss}", 150, 0, 150);
    TH1D *measuredEnFraction = new TH1D("measuredEnFraction", "E_{loss}/E_{gen} "+suffix+";E_{loss}/E_{gen}", 100, 0, 2);

    // Define utility variables
    int nShowerCrystals = pow(2*showerSize_ + 1, 2);
    
    // Events loop begins here
    cout<<"--- Event loop begins"<<endl;

    for (int i=0; i<nEvents; ++i){
        
        trTracks->GetEvent(i);
        trCalo->GetEvent(i);

        if(i%1000 == 0) cout<<"event "<<i<<endl;

        // Event Selection and Primary Tracks
        int eleIndex = -1;

        for(int j=0; j<nTracks; ++j){
            if(t_trackID[j]==2){
                eleIndex = j;
                break;
            }
        }
        if(eleIndex == -1) continue;
        if(t_PDG[eleIndex] != 11) continue;

        // --- ELOSS ANALYSIS
        measuredEnergy->Fill(c_totalEnergy/1000); // convert MeV in GeV
        measuredEnFraction->Fill(c_totalEnergy/t_kinEnergy[eleIndex]);

        // --- CRYSTAL ANALYSIS
        double maxDeposit = 0;
        int maxCrystal = -1;

        // Search for crystal with the largest energy deposit
        for(int j=0; j<c_ncrystals; ++j){
            if(c_crystalType[j] == 0) continue;
            double en = c_deposit[j];
            if(en>maxDeposit){
                maxDeposit = en;
                maxCrystal = j;
            }
        }

        double maxCalX = c_calposx[maxCrystal];
        double maxCalY = c_calposy[maxCrystal];
        int    maxType = c_crystalType[maxCrystal]; // 1 == inner calo, 2 == outer calo

        // Define shower
        double size = (maxType == 1 ? innerSize_ : outerSize_);
        double showerSizeMM = (float)showerSize_*size + size*0.05; // add 5% tolerance

        // Measure total shower energy
        double xMean = 0.;
        double yMean = 0.;
        double norm  = 0.;
        double etot = 0.;
        int    nshower = 0;
        vector<int> showerCrystals;

        for(int j=0; j<c_ncrystals; ++j){
            if(c_deposit[j] == 0) continue;
            if(c_crystalType[j] == 0) continue;
            if(abs(c_calposx[j] - maxCalX) > showerSizeMM ) continue;
            if(abs(c_calposy[j] - maxCalY) > showerSizeMM ) continue;
            nshower++;
            etot += c_deposit[j];
            showerCrystals.push_back(j); // store shower crystals indices
        }

        if(nshower < nShowerCrystals) cout<<"nshower < nShowerCrystals: "<<nshower<<" < "<<nShowerCrystals<<endl; // shower not contained in inner calo
        if(nshower > nShowerCrystals) cout<<"nshower > nShowerCrystals: "<<nshower<<" > "<<nShowerCrystals<<endl; // shower not contained in outer calo

        for(auto it:showerCrystals){
            double en = c_deposit[it];
            double w = TMath::Max(0., w0_ + log(en/etot));
            norm  += w;
            xMean += w*c_calposx[it];
            yMean += w*c_calposy[it];
        }

        xMean /= norm;
        yMean /= norm;

        // Compute theta calo
        double r = sqrt(pow(xMean,2)+pow(yMean,2));
        double dist = cCalZentry - t_zVertex[eleIndex];
        double theta = TMath::ATan(r/dist);
        double reso = (theta-t_theta[eleIndex]);

        measuredAngle->Fill(theta*1000);
        angleResolution->Fill(reso*1000);

    }

    TCanvas *c1 = new TCanvas("c1","c1",1000,600);
    c1->Divide(2,2);
    c1->cd(1);
    gPad->SetGrid();
    measuredAngle->SetMarkerStyle(20);
    measuredAngle->SetMarkerSize(.5);
    measuredAngle->Draw("PE");
    c1->cd(2);
    gPad->SetGrid();
    angleResolution->SetMarkerStyle(20);
    angleResolution->SetMarkerSize(.5);
    angleResolution->Draw("P0E");
    TF1 *fit = new TF1("fit", "[3]*exp(-0.5*((x-[0])/[1])**2)+[4]*exp(-0.5*((x-[0])/[2])**2)", -(maxTheta-minTheta)/2, +(maxTheta-minTheta)/2);
    fit->SetParName(0,"Mean");
    fit->SetParName(1,"Sigma1");
    fit->SetParName(2,"Sigma2");
    fit->SetParName(3,"Constant1");
    fit->SetParName(4,"Constant2");
    fit->SetParameter(0,0.);
    fit->SetParameter(1,.1);
    fit->SetParameter(2,.1);
    fit->SetParameter(3,10.);
    fit->SetParameter(4,10.);

    angleResolution->Fit("fit", "Q");

    cout<<endl<<"Standard Deviation = "<<angleResolution->GetStdDev()<<endl;

    c1->cd(3);
    gPad->SetGrid();
    measuredEnergy->SetMarkerStyle(20);
    measuredEnergy->SetMarkerSize(.5);
    measuredEnergy->Draw("PE");
    c1->cd(4);
    gPad->SetGrid();
    measuredEnFraction->SetMarkerStyle(20);
    measuredEnFraction->SetMarkerSize(.5);
    measuredEnFraction->Draw("PE");

    c1->Print("plots/calo" + suffix + ".png");

    return;
}
