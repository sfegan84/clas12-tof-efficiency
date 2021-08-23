// ECAL cuts remove last 2 counters for FTOF1A and FTOF1B


//#include "DC_Fiducial_Cuts_CLAS12.cxx"
#include <iostream>
#include <cstdlib>
#include <chrono>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <vector>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TPaletteAxis.h>
#include <TBenchmark.h>
#include "clas12reader.h"
#include <stdlib.h>
#include "Riostream.h"
#include "TLine.h"
#include "TVirtualPad.h"
#include "TClass.h"
#include "TVirtualX.h"
#include "TMath.h"
#include "TStyle.h"

using namespace clas12;

void SetLorentzVector(TLorentzVector &p4,clas12::region_part_ptr rp){
  p4.SetXYZM(rp->par()->getPx(),rp->par()->getPy(),
  rp->par()->getPz(),p4.M());

}

//Provide an input filename. If none provided, a hard coded filename is used (void CTOF_eff() version below)
void CTOF_eff(TString inFileName){

  TString inputFile = inFileName;

  // Creating a TChain of all the input files
  TChain fake("hipo");
  // Adding the different input files to the TChain
  fake.Add(inputFile.Data());
  // fake.Add(inputFile2.Data());


  auto db=TDatabasePDG::Instance();
  TLorentzVector beam(0,0,10.6,10.6); // beam Px,Py,Pz,E
  TLorentzVector target(0,0,0,db->GetParticle(2212)->Mass()); // target Px,Py,Pz,E
  TLorentzVector el(0,0,0,db->GetParticle(11)->Mass()); // scattered e^- Px,Py,Pz,E
  TLorentzVector pr(0,0,0,db->GetParticle(2212)->Mass()); // proton Px,Py,Pz,E
  //TLorentzVector pip(0,0,0,db->GetParticle(211)->Mass()); // pi^+ Px,Py,Pz,E
  TLorentzVector pim(0,0,0,db->GetParticle(-211)->Mass()); // pi^+ Px,Py,Pz,E

  //TLorentzVector misspim;
TLorentzVector misspip;
// Negative particle set to pi^-
//TLorentzVector pim;
 TLorentzVector pip;

// Variables for 2 pi events
// Particle numbers for 2pi events
Int_t negative, positive, nonelectron;
// Creating variables for comparing detected and missing pion
Double_t DeltaP, DeltaTheta, DeltaPhi;

  TVector3 V3_q;

  Int_t part_pid;

  Double_t nu, W_var;

  // x,y,z positions for 3 DC layers
  Double_t part_DC_c1x,part_DC_c1y,part_DC_c1z;
  Double_t part_DC_c2x,part_DC_c2y,part_DC_c2z;
  Double_t part_DC_c3x,part_DC_c3y,part_DC_c3z;

  // Retrieving list of files
  auto files=fake.GetListOfFiles();
  // Gets total events in all files for run dependence binning
  Int_t Bins = files->GetEntries();
  // Output file location and name

  TString outFileName( inFileName(0,inFileName.Length()-5) + "_eff.root"); //trim '.hipo' and add '.root' for output file name

  cout << outFileName << endl;

  //TFile fileOutput1("/u/home/sfegan/CTOF_Efficiency_RGA_FALL2018_testMissPip5000s.root","recreate");
  //TFile fileOutput1("/home/stuart/CLAS/CTOF/CTOF_Efficiency_RGA_FALL2018_test5038MissPip.root","recreate");
  TFile fileOutput1(outFileName,"recreate");

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Create histograms here

  auto *h_z_vertex=new TH1D("h_z_vertex","z vertex",100,-15,15);
  auto *h_z_vertex_CD=new TH1D("h_z_vertex_CD","z vertex (central detector)",100,-15,15);

  auto *h_beta=new TH1D("h_beta","particle velocity",500,-1.0,1.5);
  auto *h_beta2=new TH1D("h_beta2","particle velocity",500,-1.0,1.5);
  auto *h_beta_mom=new TH2D("h_beta_mom","beta vs track momentum",200,0,10,100,-1.0,1.5);
  auto *h_beta_mom_cut=new TH2D("h_beta_mom_cut","beta vs track momentum (post cuts)",200,0,10,100,-1.0,1.5);
  auto *h_beta_mom_cut2=new TH2D("h_beta_mom_cut2","beta vs track momentum (post cuts)",200,0,10,100,-1.0,1.5);

  auto *h_verttime=new TH1D("h_verttime","vertex time",1000,-200,200);
  auto *h_time=new TH1D("h_time","hit time",900,-100,200);
  auto *h_time2=new TH1D("h_time2","hit time (cut)",900,-100,200);

  auto *h_Wvar=new TH1D("h_Wvar","W (CD)",400,-3,5);

  auto *h_trackPath=new TH1D("h_trackPath","Track Path Length (cm)",550,-10,100);

  auto *h_CTOF_E=new TH1D("h_CTOF_E","Energy deposition in CTOF",200,0.0,50.0);

  // Arrays of histograms [layer][charge][sector]
  //TH3F *h_Trajectories[3][2][6]; // Trajectories from DC
  //TH3F *h_Tracks[3][2][6]; // Trajectories from DC with energy deposited in FTOF
  //TH3F *h_Efficiency[3][2][6]; // Tracks divided by Trajectories

  TH2F *h_Traj_CD[2]; // Trajectories from Central Tracker
  TH2F *h_Trk_CD[2]; // Trajectories from Central Tracker with energy deposited in CND
  TH2F *h_Eff_CD[2]; // Tracks divided by Trajectories

  //paddle version
  TH2F *h_Traj_CD_paddle[2]; // Trajectories from Central Tracker
  TH2F *h_Trk_CD_paddle[2]; // Trajectories from Central Tracker with energy deposited in CND
  TH2F *h_Eff_CD_paddle[2]; // Tracks divided by Trajectories

  TH2F *h_Scint_XY[2];  //X,Y scintillator hit positions
  TH3F *h_Scint_XYZ[2]; //X,Y,Z scintillator hit positions

  TH2F *h_CND_XY[2];  //X,Y scintillator hit positions
  TH3F *h_CND_XYZ[2]; //X,Y,Z scintillator hit positions

  //Track NDF dependent histos
  TH2F *h_Traj_CD_NDF[2][10]; // Trajectories from Central Tracker
  TH2F *h_Trk_CD_NDF[2][10]; // Trajectories from Central Tracker with energy deposited in CND
  TH2F *h_Eff_CD_NDF[2][10]; // Tracks divided by Trajectories

  TH1F *h_Trk_NDF_CD[2]; // Trajectories from Central Tracker

// 2 pi event histograms
auto* hmass=new TH1F("pipmmass","Missing Mass e' p #pi^{-};MM(e'p#pi^{-}) [GeV];Counts",200,-1,1);
auto* hmass2=new TH1F("pipmmass2","Missing Mass e' p #pi^{-} (post cuts);MM(e'p#pi^{-}) [GeV];Counts",200,-1,1);
auto* hdeltaP=new TH1F("DeltaMomentum","Momentum difference of #pi^{+} detected and reconstructed;#Delta P [GeV];Counts",400,-2,2);
auto* hdeltaTheta=new TH1F("DeltaTheta","#theta difference of #pi^{+} detected and reconstructed;#Delta #theta [deg];Counts",360,-180,180);
auto* hdeltaPhi=new TH1F("DeltaPhi","#phi difference of #pi^{+} detected and reconstructed;#Delta #phi [deg];Counts",360,-180,180);

auto* h_el_thetaPhi=new TH2D("el_thetaPhi","#theta versus #phi, electron;#phi [deg];#theta [deg]",360,-180,180,180,0,180);
auto* h_prot_thetaPhi=new TH2D("prot_thetaPhi","#theta versus #phi, Proton;#phi [deg];#theta [deg]",360,-180,180,180,0,180);
auto* h_pipl_thetaPhi=new TH2D("pipl_thetaPhi","#theta versus #phi, #pi^{+};#phi [deg];#theta [deg]",360,-180,180,180,0,180);
auto* h_pimi_thetaPhi=new TH2D("pimi_thetaPhi","#theta versus #phi, #pi^{-};#phi [deg];#theta [deg]",360,-180,180,180,0,180);


  // Looping over negative and positive particles in Central Detector
  for(int i_charge=0;i_charge<2;i_charge++){
    
    //create a string which we can append integers to, which allows us to define a number of histograms in a for loop
    // Histogram names
    ostringstream TrajCD_name_stream;
    ostringstream TracksCD_name_stream;
    ostringstream TrajCD_paddle_name_stream;
    ostringstream TracksCD_paddle_name_stream;
    
    // Histogram Titles
    ostringstream TrajCD_title_stream;
    ostringstream TracksCD_title_stream;
    
    // Defining the histogram name strings
    TrajCD_name_stream<<"h_TrajCD_Charge_"<<i_charge;
    TracksCD_name_stream<<"h_TracksCD_Charge_"<<i_charge;
    TrajCD_paddle_name_stream<<"h_TrajCD_paddle_Charge_"<<i_charge;
    TracksCD_paddle_name_stream<<"h_TracksCD_paddle_Charge_"<<i_charge;
    
    // Defining the histogram title strings
    TrajCD_title_stream<<"Trajectories CTOF Charge "<<2*i_charge-1<<"; z position of hit [cm]; Polar Angle [degrees]";
    
    //convert the stringstream to a string and define our histograms in each element of the array
    h_Traj_CD[i_charge] = new TH2F(TrajCD_name_stream.str().c_str(),"", 70,-70,70, 56, -210, 210);
    h_Trk_CD[i_charge] = new TH2F(TracksCD_name_stream.str().c_str(),"", 70,-70,70, 56, -210 , 210);

    h_Traj_CD_paddle[i_charge] = new TH2F(TrajCD_paddle_name_stream.str().c_str(),"", 70, -70, 70, 50, -1, 49);
    h_Trk_CD_paddle[i_charge] = new TH2F(TracksCD_paddle_name_stream.str().c_str(),"", 70, -70, 70, 50, -1 , 49);

//    h_Trk_NDF_CD[i_charge] = new TH1F(Form("h_TracksNDF_CD_%d",((int)i_charge).c_str()),"",50,0,50);
    h_Trk_NDF_CD[i_charge] = new TH1F(Form("h_TracksNDF_CD_%d",i_charge),"",50,0,50);
    
    h_Scint_XY[i_charge] = new TH2F(Form("h_Scint_XY_%d",(2*i_charge-1)),"", 100,-50,50, 100, -50 , 50);
    h_Scint_XYZ[i_charge] = new TH3F(Form("h_Scint_XYZ_%d",(2*i_charge-1)),"", 100,-50, 50, 100, -50 , 50, 100, -50 , 50);

    h_CND_XY[i_charge] = new TH2F(Form("h_Scint_XY_%d",(2*i_charge-1)),"", 100,-50,50, 100, -50 , 50);
    h_CND_XYZ[i_charge] = new TH3F(Form("h_Scint_XYZ_%d",(2*i_charge-1)),"", 100,-50, 50, 100, -50 , 50, 100, -50 , 50);

    // Setting the title for each histogram as it did not work when put in the lines above
    h_Traj_CD[i_charge]->SetTitle(TrajCD_title_stream.str().c_str());
    h_Trk_CD[i_charge]->SetTitle(TracksCD_title_stream.str().c_str());
    h_Trk_NDF_CD[i_charge]->SetTitle(Form("h_TracksNDF_CD_%d",(2*i_charge-1)));

    for(int i_NDF=0;i_NDF<10;i_NDF++){
      //Track NDF dependent histos
      h_Traj_CD_NDF[i_charge][i_NDF] =  new TH2F(Form("h_TrajCD_Charge_%d_NDF_%d",(2*i_charge-1),i_NDF),"", 70,-70,70, 56, -210, 210); // Trajectories from Central Tracker
      h_Trk_CD_NDF[i_charge][i_NDF] =  new TH2F(Form("h_TracksCD_Charge_%d_NDF_%d",(2*i_charge-1),i_NDF),"", 70,-70,70, 56, -210, 210); // Trajectories from Central Tracker

      h_Traj_CD_NDF[i_charge][i_NDF]->SetTitle(Form("Trajectories CTOF Charge %d NDF %d",(2*i_charge-1),i_NDF));
      h_Trk_CD_NDF[i_charge][i_NDF]->SetTitle(Form("Tracks CTOF Charge %d NDF %d",(2*i_charge-1),i_NDF));

      //old binning
      //      h_Traj_CD_NDF[i_charge][i_NDF] =  new TH2F(Form("h_TrajCD_Charge_%d_NDF_%d",(2*i_charge-1),i_NDF),"", 160,-80,80, 200, -200, 200); // Trajectories from Central Tracker
      //      h_Trk_CD_NDF[i_charge][i_NDF] =  new TH2F(Form("h_TracksCD_Charge_%d_NDF_%d",(2*i_charge-1),i_NDF),"", 160,-80,80, 200, -200, 200); // Trajectories from Central Tracker
      //h_Eff_CD_NDF[i_charge][i_NDF]; // Tracks divided by Trajectories
    }

  }


//   // Looping over the FTOF layers   //only one layer in CTOF
//   for(Int_t i_detector=0;i_detector<3;i_detector++){
//     // Looping over negative and positive particles
//     for(Int_t i_charge=0;i_charge<2;i_charge++){
//       // Looping over the sectors   //no sectors in CTOF, start with scintillator bar?
//       for(Int_t i_sector=0;i_sector<6;i_sector++){

//         //create a string which we can append integers to, which allows us to define a number of histograms in a for loop
//         // Histogram names
//         ostringstream Traj_name_stream;
//         ostringstream Tracks_name_stream;

//         // Histogram Titles
//         ostringstream Traj_title_stream;
//         ostringstream Tracks_title_stream;

//         // Defining the histogram name strings
//         Traj_name_stream<<"h_Traj_Det_"<<i_detector<<"_Charge_"<<i_charge<<"_Sec_"<<i_sector;
//         Tracks_name_stream<<"h_Tracks_Det_"<<i_detector<<"_Charge_"<<i_charge<<"_Sec_"<<i_sector;

//         // Defining the histogram title strings
//         if (i_detector==0) Traj_title_stream<<"Trajectories FTOF1A Charge "<<2*i_charge-1<<" Sec "<<i_sector+1<<"; Run no.; L [cm]";
//         else if (i_detector==1) Traj_title_stream<<"Trajectories FTOF1B Charge "<<2*i_charge-1<<" Sec "<<i_sector+1<<"; Run no.; L [cm]";
//         else if (i_detector==2) Traj_title_stream<<"Trajectories FTOF2 Charge "<<2*i_charge-1<<" Sec "<<i_sector+1<<"; Run no.; L [cm]";

//         //convert the stringstream to a string and define our histograms in each element of the array
//         h_Trajectories[i_detector][i_charge][i_sector] = new TH3F(Traj_name_stream.str().c_str(),"", Bins,0,Bins,500,0,500, 200, -225, 225);
//         h_Tracks[i_detector][i_charge][i_sector] = new TH3F(Tracks_name_stream.str().c_str(),"", Bins,0,Bins,500,0,500, 200, -225 , 225);

//         // Setting the title for each histogram as it did not work when put in the lines above
//         h_Trajectories[i_detector][i_charge][i_sector]->SetTitle(Traj_title_stream.str().c_str());
//         h_Tracks[i_detector][i_charge][i_sector]->SetTitle(Tracks_title_stream.str().c_str());
//       }
//     }
//   }

  // Distance between trajectory point and scintillator hit
  auto* h_radia_residual_CD = new TH1D("h_radia_residual_CD","Distance Between CVT and CTOF hits",150,0,30);
  auto* h_radia_CTOF_CND = new TH1D("h_radia_CTOF_CND","Distance Between CND and CTOF hits",150,0,50);
  auto* h_path_CTOF_CND = new TH1D("h_path_CTOF_CND","Path Length Difference Between Track at CTOF and CND",150,0,50);

  auto* h_momentum_CD = new TH1D("h_momentum_CD","Particle Momentum (Central Detector)",200,0,10);
  auto* h_momentum_CD_cut = new TH1D("h_momentum_CD_cut","Particle Momentum (Central Detector)",200,0,10);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Creating variables for different layers of FTOF

  // Positions and angles
  Double_t x_1a, x_1b, x_2, y_1a, y_1b, y_2, z_1a, z_1b, z_2; // DC trajectory x,y,z positions
  Double_t d_1a, d_1b, d_2; // Distance to xy position (used for calculating perpendicular distance)
  Double_t x_FTOF, y_FTOF, z_FTOF; // Scintillator x,y,z positions
  Double_t L_1a, L_1b, L_2; // Perpendicular distance to sector in lab frame
  Double_t L_det_1a, L_det_1b, L_det_2; // Perpendicular distance to sector in detector frame

  Double_t alpha_1a, alpha_1b, alpha_2; // angle to hit (used to determine sector)
  Double_t L_Perp_1a, L_Perp_1b, L_Perp_2; // Distance along counter
  Double_t L_Theta; // Angle at middle of sector (check if hit is left or right of the middle of the sector)
  Double_t radia_residual; // Distance between trajectory bank and scintillator hit
  Double_t radia_CTOF_CND; // Distance between CTOF and CND hits

  Int_t TrackNDF; //Number of degrees of freedom in CD track

  Int_t paddleNo; //Arbitrary index to permit plotting by CTOF paddle

  Double_t x_CD, y_CD, z_CD; // Central Tracker trajectory x,y,z positions
  Double_t x_CND, y_CND, z_CND; // Central Detector x,y,z hit positions
  Double_t trackMom;  //Momentum of CD track
  Double_t trackBeta;  //beta of CD track
  Double_t trackBetaCalc;  //calculated beta of CD track
  Double_t stTime;  //start time of CD track
  Double_t trackTime;  //vertex time of CD track
  Double_t trk_px, trk_py, trk_pz; //Momentum components of track
  Double_t d_CD; // Distance to xy position (used for calculating perpendicular distance)
  Double_t x_CTOF, y_CTOF, z_CTOF; // Scintillator x,y,z positions

  Double_t alpha_CD, alpha_CTOF; // angle to hit (used to determine sector)
  Double_t L_Perp_CD; // Distance along counter

  Double_t CND_path, CTOF_path; //path lengths for CND and CTOF

  // Run Number
  Int_t runno; // runno as a integer
  vector<TString> v_Runno; // runno as a vector of string
  Int_t Binno=0; // Count the number of files this corresponds to the number of x bins
  // Status
  Int_t Status, Calorimeter_Hits; // Status is used to find out if there is a ECAL hit
  Int_t StatusCD, Cal_Hits; // StatusCD is used to find out if there is a CND hit
  Int_t StatusTrack;
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Looping over data

  //Loop over files
  for(Int_t i=0;i<files->GetEntries();i++){

    //create the event reader
    clas12reader c12(files->At(i)->GetTitle());

    Binno++; // Count the number of files, therefore the number of x bins

    //c12.setEntries(1E5);
    while(c12.next()==true){
      
      auto particles = c12.getDetParticles();
      
      negative = 0;
      positive = 0;
      nonelectron = 0;
      DeltaP = 0;
      DeltaTheta = 0;
      DeltaPhi = 0;
      
      auto electrons=c12.getByID(11);
      auto protons=c12.getByID(2212);
      //auto pips=c12.getByID(211);
      auto pims=c12.getByID(-211);
      
      for(auto& p : particles){

	//nonelectron reappropriated to refer to 'missing' particle	

        // Looking at negative particles
        if(p->par()->getCharge() < 0){
	  negative++;
          // negative particles not electron are set to pi^-
          if(p->par()->getPid() != 11){
            //nonelectron++;
            //pim.SetXYZM(p->par()->getPx(),p->par()->getPy(),p->par()->getPz(),db->GetParticle(-211)->Mass());
          }
        }
        // Looking at positive particles
        else if(p->par()->getCharge() > 0){
	  positive++;
          // positive particles not proton are set to pion
          if(p->par()->getPid() != 2212){
            nonelectron++;
            pip.SetXYZM(p->par()->getPx(),p->par()->getPy(),p->par()->getPz(),db->GetParticle(211)->Mass());
          }
	}
      }
      
      
      // Getting 2pi events
       if(nonelectron == 1 && electrons.size() == 1 && pims.size() == 1 && protons.size() == 1){
 	SetLorentzVector(el,electrons[0]);
         SetLorentzVector(pr,protons[0]);
         //SetLorentzVector(pip,pips[0]);
         SetLorentzVector(pim,pims[0]);
         misspip = beam + target - el - pim - pr;
         hmass->Fill(misspip.M2());
	
         // Cut on missing mass of the pi^-
         if((misspip.M2() > -0.1) && (misspip.M2() < 0.2)){
           DeltaP = misspip.Rho() - pip.Rho();
           DeltaTheta = TMath::RadToDeg()* (misspip.Theta() - pip.Theta());
           DeltaPhi = TMath::RadToDeg()* (misspip.Phi() - pip.Phi());
	   
           // Plotting pi^- variables
           hdeltaP->Fill(DeltaP);
           hdeltaTheta->Fill(DeltaTheta);
           hdeltaPhi->Fill(DeltaPhi);
	   
	   //(initially) loose cuts
	   if(fabs(DeltaP) > 0.3) continue;
	   if(fabs(DeltaTheta) > 10) continue;
	   if(fabs(DeltaPhi) > 10) continue;
	   
	   h_el_thetaPhi->Fill(TMath::RadToDeg()*el.Phi(), TMath::RadToDeg()*el.Theta());
	   h_prot_thetaPhi->Fill(TMath::RadToDeg()*pr.Phi(), TMath::RadToDeg()*pr.Theta());
	   h_pipl_thetaPhi->Fill(TMath::RadToDeg()*pip.Phi(), TMath::RadToDeg()*pip.Theta());
	   h_pimi_thetaPhi->Fill(TMath::RadToDeg()*pim.Phi(), TMath::RadToDeg()*pim.Theta());
	   
	   nu = -((el - beam).E());
	   V3_q = (beam-el).Vect();
	   W_var = TMath::Sqrt((0.938+nu)*(0.938+nu)-V3_q*V3_q);
	   
	   hmass2->Fill(misspip.M2());
//{
//{	  


	  //second loop
	  //Set the particle index to 0 and loop through the particles
	  int pindex=0;
	  //std::cout << "event" << std::endl;
	  for(auto& p : particles){
	    
	    //get information from the different detectors
	    switch(p->getRegion()) {

	      //add a skip of electrons in FT
	    case FD:
	      pindex++;
	      
	    case CD :
	      
	      //std::cout << p->par()->getCharge() << std::endl;

	      //std::cout << "central detector" << endl;
	      // Increase the particle index with each loop of the particles
	      pindex++;
	      
	      //if(pindex==1){
		//std::cout << "first CD particle" << p->par()->getPid() << std::endl;
	      //}

	      //Ignore the first particle (trigger, if it falls in the CTOF) and any neutrals
	      if (pindex==1 || p->par()->getCharge()==0) continue;
	      //just ignore neutrals
	      //if (p->par()->getCharge()==0) continue;

	      //if(p->par()->getPid() ==0) continue;
	      
	      runno = c12.runconfig()->getRun(); // Getting the run number
	      StatusCD = p->par()->getStatus(); // Getting the status
	      
	      // Getting the hits in CND
	      Cal_Hits = (StatusCD / 4000);//%10;
	      //track status
	      StatusTrack = ((StatusCD%4000)%100)/10; //Scintillator hits i.e. CND hits
	      //std::cout<< "Cal_Hits = " << Cal_Hits << " Scintillator Hits = " << StatusTrack <<std::endl;
	      
	      //beta, timing?
	      trackBeta = p->par()->getBeta();
	      trackBetaCalc = (p->sci(CTOF)->getPath()/((p->sci(CTOF)->getTime())-(p->par()->getVt())))/29.9792;
	      h_beta->Fill(trackBetaCalc);
	      
	      stTime = p->par()->getVt();
	      h_verttime->Fill(stTime);
	      trackTime = p->sci(CTOF)->getTime();
	      h_time->Fill(trackTime);
	      h_trackPath->Fill(p->sci(CTOF)->getPath());
	      
	      // Checking the z vertex before applying a cut
	      h_z_vertex_CD->Fill(p->par()->getVz());
	      if(p->par()->getVz() > 2 || p->par()->getVz() < -9)continue;  
	      
	      
	      //Compute track momentum from components
	      trk_px = p->par()->getPx();
	      trk_py = p->par()->getPy();
	      trk_pz = p->par()->getPz();
	      trackMom = TMath::Sqrt((trk_px*trk_px)+(trk_py*trk_py)+(trk_pz*trk_pz));
	      
	      h_beta2->Fill(trackBetaCalc);	  
	      h_beta_mom->Fill(trackMom,trackBetaCalc);
	      
	      //Fill Histogram with and without cut
	      h_momentum_CD->Fill(trackMom);
	      //this currently applies to all tracks, we don't need all tracks to reach CND, just the ones we're measuring efficiency of
	      //if(trackMom<0.3) continue;  //cut at particle momentum less than 300 MeV
	      //if(trackMom<0.4) continue;  //cut at particle momentum less than 400 MeV, so particle reaches CND
	      h_momentum_CD_cut->Fill(trackMom);
	      
	      
	      //if((trackBetaCalc<0.2)||(trackBetaCalc>1.2)) continue;
	      //if(trackBeta<0.2) continue;
	      //if(trackBeta>1.2) continue;
	      
	      
	      h_CTOF_E->Fill(p->sci(CTOF)->getEnergy());
	      
	      
	      //Skips any particle with no CD status
	      //if(Cal_Hits!=0){
	      //Skips any particle with no CD status and no CND hit
	      if((Cal_Hits!=0)&&(StatusTrack!=0)){
		
		h_beta_mom_cut->Fill(trackMom,trackBeta);
		h_beta_mom_cut2->Fill(trackMom,trackBetaCalc);
		
		h_Wvar->Fill(W_var);
		
		//if(W_var < 1) continue;
		
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//CTOF
		
		// Here you can put cuts on the particles you are looking at
		
		if((p->trk(CVT)->getDetector() ==5)&&(p->trk(CVT)->getNDF() >=0)){
		  
		  TrackNDF = (Int_t)(p->trk(CVT)->getNDF());
		  //std::cout << "Track NDF = "<< TrackNDF << endl;
		  
		  if(TrackNDF < 2) continue;  //skip events with track NDF less than 2
		  
		  if(p->par()->getCharge()>0){//positive tracks
		    h_Trk_NDF_CD[1]->Fill(p->trk(CVT)->getNDF());
		  }
		  else if(p->par()->getCharge()<0){//negative tracks
		    h_Trk_NDF_CD[0]->Fill(p->trk(CVT)->getNDF());
		  }
		  else{
		    continue;  //skip neutrals
		  }
		  //}
		  
		  h_time2->Fill(trackTime);
		  
		  //CND hits?
		  if(p->sci(CND)->getDetector()==3){
		    x_CND = p->sci(CND)->getX();
		    y_CND = p->sci(CND)->getY();
		    z_CND = p->sci(CND)->getZ();
		  }
		  
		  
		  //std::cout << "something " << trackMom << std::endl;
		  //CTOF,1 is inner layer, 2 middle, 3 outer
		  
		  if(p->traj(CTOF,1)->getDetector()==4 && p->traj(CTOF,1)->getLayer()==1){
		    // Getting the x-, y- and z- co-ordinates from CVT, i.e the track
		    x_CD =  p->traj(CTOF, 1)->getX();
		    y_CD =  p->traj(CTOF, 1)->getY();
		    z_CD = p->traj(CTOF, 1)->getZ();

// 		    //first order ``Fiducial'' cut, -25 < z_CD < 30 cm
// 		    if((z_CD < -25) || (z_CD > 30)){
// 		      continue;
// 		    }
		    
		    // Getting x-, y- and z- co-ordinates from CTOF hit, i.e. the scintillator hit
		    if(p->sci(CTOF)->getEnergy()>0){
		      x_CTOF =  p->sci(CTOF)->getX();
		      y_CTOF =  p->sci(CTOF)->getY();
		      z_CTOF =  p->sci(CTOF)->getZ();
		      
		      // Distance between 'track' (probably CVT) and CTOF co-ordinates
		      radia_residual = sqrt(pow(x_CD-x_CTOF,2) + pow(y_CD-y_CTOF,2) + pow(z_CD-z_CTOF,2));
		      h_radia_residual_CD->Fill(radia_residual); //do we want to cut on this?
		      
		      if(radia_residual>6) continue;
		      
		      // Distance between CND and CTOF co-ordinates
		      radia_CTOF_CND = sqrt(pow(x_CND-x_CTOF,2) + pow(y_CND-y_CTOF,2) + pow(z_CND-z_CTOF,2));
		      h_radia_CTOF_CND->Fill(radia_CTOF_CND);

		      //if(radia_CTOF_CND>30) continue;
		      
		      h_path_CTOF_CND->Fill((p->sci(CTOF)->getPath())-(p->sci(CND)->getPath()));
		      
		      //if((p->sci(CTOF)->getPath())-(p->sci(CND)->getPath())>33) continue;

		      //std::cout << ((p->sci(CTOF)->getPath())-(p->sci(CND)->getPath())) <<std::endl;
		      
		    }
		    
		    
// 		h_Scint_XY[0]->Fill(x_CTOF,y_CTOF);
// 		h_Scint_XYZ[0]->Fill(x_CTOF,y_CTOF,z_CTOF);

// 		h_Scint_XY[1]->Fill(x_CD,y_CD);
// 		h_Scint_XYZ[1]->Fill(x_CD,y_CD,z_CD);

		    //lets start by converting x,y,z to r,phi,z and visualising the hits
		    //r_CTOF =  sqrt(pow(x_CTOF,2) + pow(y_CTOF,2));
		    //phi_CTOF = TMath::RadToDeg()*atan2(y_CTOF,x_CTOF);
		    ////z_CTOF = z_CTOF;		
		    
		    // Calculating d, distance to hit from (0,0) to (x,y)
		    d_CD = sqrt(pow(x_CD,2) + pow(y_CD,2));   //this is r
		    
		    // Calculating alpha, angle from (0,0) to hit, this is phi
		    //alpha_CD = TMath::RadToDeg()*atan(y_CD/x_CD);  //Check, and fix if necessary, the use of atan, replacing with atan2(y_CD,x_CD)
		    alpha_CD = TMath::RadToDeg()*atan2(y_CD,x_CD);   //Replacing atan(y/x) with atan2(y_CD,x_CD)
		    alpha_CTOF = TMath::RadToDeg()*atan2(y_CTOF,x_CTOF);
		    
		    //attempt to calculate a "paddle index"
		    paddleNo = (int)(floor((alpha_CD+180)/7.5));
		    
		    
		    //z_CD = z_CD;
		    
		    // Positive particles
		    if(p->par()->getCharge()>0){
		      ////h_Traj_CD[1]->Fill(i,L_det_1a, L_Perp_1a);
		      if(p->par()->getPid()==2212){
			h_Traj_CD[1]->Fill(z_CTOF, alpha_CTOF);
		      }
		      if(p->par()->getPid()==211){
			h_Traj_CD[0]->Fill(z_CTOF, alpha_CTOF);
		      }

		      
		      h_Traj_CD_NDF[1][0]->Fill(z_CD, alpha_CD);
		      //h_Traj_CD_NDF[1][TrackNDF]->Fill(z_CD, alpha_CD);
		      
		      for(int ii=1;ii<10;ii++){
			if(TrackNDF>=ii){
			  h_Traj_CD_NDF[1][ii]->Fill(z_CD, alpha_CD);
			}
		      }
		      
		    }
		    
		    // Negative particles
		    else if(p->par()->getCharge()<0){
		      ////h_Traj_CD[0]->Fill(i,L_det_1a, L_Perp_1a);
		      //h_Traj_CD[0]->Fill(z_CTOF, alpha_CTOF);       
		      
		      h_Traj_CD_NDF[0][0]->Fill(z_CD, alpha_CD);           
		      
		      //h_Traj_CD_NDF[0][TrackNDF]->Fill(z_CD, alpha_CD);
		      
		      for(int ii=1;ii<10;ii++){
			if(TrackNDF>=ii){
			  h_Traj_CD_NDF[0][ii]->Fill(z_CD, alpha_CD);
			}
		      }
		      
		    }
		    
		    // Check if there is energy deposited on the scintillator
		    if(p->sci(CTOF)->getEnergy()>0){
		      // Positive particles
		      if(p->par()->getCharge()>0){
			////h_Trk_CD[1]->Fill(i,L_det_1a, L_Perp_1a);
			if(p->par()->getPid()==2212){
			  h_Trk_CD[1]->Fill(z_CTOF, alpha_CTOF);
			}
			if(p->par()->getPid()==211){
			  h_Trk_CD[0]->Fill(z_CTOF, alpha_CTOF);
			}
			//h_Trk_CD[1]->Fill(z_CTOF, alpha_CTOF);
			h_Trk_CD_paddle[1]->Fill(z_CTOF, paddleNo);
			
			h_Trk_CD_NDF[1][0]->Fill(z_CD, alpha_CD);
			
			//h_Trk_CD_NDF[1][TrackNDF]->Fill(z_CD, alpha_CD);
			
			for(int ii=1;ii<10;ii++){
			  if(TrackNDF>=ii){
			    h_Trk_CD_NDF[1][ii]->Fill(z_CD, alpha_CD);
			  }
			}
			
		      }
		      
		      // Negative particles
		      else if(p->par()->getCharge()<0){
			////h_Trk_CD[0]->Fill(i,L_det_1a, L_Perp_1a);
			//h_Trk_CD[0]->Fill(z_CTOF, alpha_CTOF);
			
			h_Trk_CD_NDF[0][0]->Fill(z_CD, alpha_CD);
			
			//h_Trk_CD_NDF[0][TrackNDF]->Fill(z_CD, alpha_CD);
			for(int ii=1;ii<10;ii++){
			  if(TrackNDF>=ii){
			    h_Trk_CD_NDF[0][ii]->Fill(z_CD, alpha_CD);
			  }
			}
		      }

		    }
		  }
		}
	      }
	    }
	  }	  
	}
      }
    }
    
    v_Runno.push_back(to_string(runno)); // Converting runno integer to a string
  }



  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Determining Efficiencies - Make it a function

  for(Int_t i_charge=0;i_charge<2;i_charge++){
    h_Eff_CD[i_charge] = (TH2F*)h_Trk_CD[i_charge]->Clone(Form("Efficiency_%d",i_charge));

    h_Eff_CD[i_charge]->Divide(h_Traj_CD[i_charge]);

    for(Int_t i_NDF=0;i_NDF<10;i_NDF++){
      h_Eff_CD_NDF[i_charge][i_NDF] = (TH2F*)h_Trk_CD_NDF[i_charge][i_NDF]->Clone(Form("Efficiency_%d_NDF_%d",i_charge,i_NDF));
      h_Eff_CD_NDF[i_charge][i_NDF]->Divide(h_Traj_CD_NDF[i_charge][i_NDF]);
      h_Eff_CD_NDF[i_charge][i_NDF]->SetTitle(Form("Efficiency %d NDF %d",(2*i_charge-1),i_NDF));
    }
  }

  // Looping over the FTOF layers
  for(Int_t i_detector=0;i_detector<3;i_detector++){
    // Looping over negative and positive particles
    for(Int_t i_charge=0;i_charge<2;i_charge++){

      // Looping over the sectors
      for(Int_t i_sector=0;i_sector<6;i_sector++){
        ostringstream Efficiency_name_stream;
        Efficiency_name_stream<<"h_Efficiency_Det_"<<i_detector<<"_Charge_"<<i_charge<<"_Sec_"<<i_sector;
        ostringstream Efficiency_title_stream;

        if (i_detector==0) Efficiency_title_stream<<"Efficiency FTOF1A Charge "<<2*i_charge-1<<" Sec "<<i_sector+1<<"; Run no.; L [cm]";
        else if (i_detector==1) Efficiency_title_stream<<"Efficiency FTOF1B Charge "<<2*i_charge-1<<" Sec "<<i_sector+1<<"; Run no.; L [cm]";
        else if (i_detector==2) Efficiency_title_stream<<"Efficiency FTOF2 Charge "<<2*i_charge-1<<" Sec "<<i_sector+1<<"; Run no.; L [cm]";


//        h_Efficiency[i_detector][i_charge][i_sector]=(TH3F*)h_Tracks[i_detector][i_charge][i_sector]->Clone(Efficiency_name_stream.str().c_str());


//        h_Efficiency[i_detector][i_charge][i_sector]->Divide(h_Trajectories[i_detector][i_charge][i_sector]);

//         // Looping over the x bins and changing them to string run numbers
//         for(Int_t i=0;i<Binno;i++){

//           h_Trajectories[i_detector][i_charge][i_sector]->GetXaxis()->SetBinLabel(h_Trajectories[i_detector][i_charge][i_sector]->GetXaxis()->FindBin(i),v_Runno.at(i));
//           h_Tracks[i_detector][i_charge][i_sector]->GetXaxis()->SetBinLabel(h_Tracks[i_detector][i_charge][i_sector]->GetXaxis()->FindBin(i),v_Runno.at(i));
// //          h_Efficiency[i_detector][i_charge][i_sector]->GetXaxis()->SetBinLabel(h_Efficiency[i_detector][i_charge][i_sector]->GetXaxis()->FindBin(i),v_Runno.at(i));

//         }
      }
    }
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //saving the file
  fileOutput1.Write();

}


//Function wrapper for hard coded filename
void CTOF_eff(){

  // Data files to process
  //TString inFile("/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass1/v0/dst/train/skim4/skim4_005*.hipo");
  //TString inFile("/volatile/clas12/rg-a/production/recon/fall2018/torus-1/pass1/v0/dst/train/dst/train/skim4/*.hipo");
  TString inFile("/home/stuart/CLAS/Data/skim4_00503*.hipo");
  // TString inputFile2("/volatile/clas12/rg-b/production/recon/spring2019/torus-1/pass1/v0/dst/train/inc/*.hipo");

  CTOF_eff(inFile); //call the analysis function with this filename 

}
