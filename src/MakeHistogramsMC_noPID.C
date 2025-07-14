#include <cstdlib>
#include <iostream>
#include <chrono>
#include <TFile.h>
#include <TNtuple.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TSystemDirectory.h>
#include <TH1.h>
#include <TH2.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLeaf.h>
#include "clas12reader.h"
#include "helonline.h"
#include "DCfiducialcuts.h"

using namespace clas12;

auto db=TDatabasePDG::Instance();

const int INBEND=1, OUTBEND=0;

int debug = 0;
double defval = 0;

void SetLorentzVector(TLorentzVector &p4, clas12::region_part_ptr rp, double mass){
  p4.SetXYZM(rp->par()->getPx(),rp->par()->getPy(),
             rp->par()->getPz(),mass);
  
}



double s30 = sin(TMath::Pi()/6);
double c30 = cos(TMath::Pi()/6);


bool pcalOK(clas12::region_part_ptr p){
  
  return p->cal(PCAL)->getLv()>9.0 && p->cal(PCAL)->getLw()>9.0;
  
}

//list of cuts
double cut_ECfracmin = 0.17;
double cut_pcalmin = 0.07;
double cut_evzmin_inb =  -13; // inbending    -8;//-13;
double cut_evzmax_inb =  12;  // inbending     1;//12;

double cut_evzmin_outb = -18;
double cut_evzmax_outb = 10;

double cut_evzmin = cut_evzmin_inb;
double cut_evzmax = cut_evzmax_inb;


//double cut_Wmin  = 2;
//double cut_Q2min = 1;
//double cut_Wmin = 2;
//double cut_ymax  = 0.85;

double cut_dvzmin = -20; //-5;//-20;
double cut_dvzmax = 20; //5;//20;

//removed this cut
double cut_dtime_corr = 10000;//0.75;//0.3;

double cut_HTCCmin = 2;

double sf1[7], sf2[7], sf3[7], sf4[7];
double sfs1[7], sfs2[7], sfs3[7], sfs4[7];


void toCM(TLorentzVector cm, TLorentzVector p,TLorentzVector& result){
  result = p;
  
  result.RotateZ(TMath::Pi()-cm.Phi());
  result.RotateY(cm.Theta());
  
  result.Boost(0,0,-cm.Beta());
  
}

double angle(double phi){
  double pi = TMath::Pi();
  while(phi>pi)
    phi-=2*pi;
  while(phi<-pi)
    phi+=2*pi;
  return phi;
}

int torus;

/*
bool acceptHadron(clas12::region_particle * h, clas12::region_particle * electron){
  int h_pid = h->getPid();
  
  TLorentzVector had;
  if(h_pid != 211 && h_pid != -211 && h_pid != 2212 && h_pid != -2212 && h_pid != 321 && h_pid != -321)
    return false;
  
  
  double dtime = electron->getTime()-h->getTime();

  //if(!dcOK(h,torus)) //|| had.Theta>35*3.14159/180)
  //  return false;
  
  //if(debug) cout << "pip" << endl;
  double mass = db->GetParticle(h_pid)->Mass();
  SetLorentzVector(had,h, mass);
  double h_p = had.P();
  double c = 29.9792458; //cm/ns
  auto dtime_corr =dtime-electron->getPath()/c+h->getPath()/(had.Beta()*c);
  if(abs(dtime_corr) > cut_dtime_corr)
    return false;
  
  
  //h_pid = h->par()->getPid();
  auto h_chi2pid = h->par()->getChi2Pid();
  
  //scaling constant for pid cut
  double C = .88*(h_pid == 211)+.93*(h_pid == -211)+ 1.30*(abs(h_pid)==2212)+1.0*(abs(h_pid) != 211 && abs(h_pid)!=2212);
  //if(abs(h_chi2pid) >3*C)
    //  return false;
  
  //tighter cut for pions
  //if(abs(h_pid) == 211 && h_p>2.44 && h_chi2pid >(0.00869+14.98587*exp(-h_p/1.18236)+1.81751*exp(-h_p/4.86394))*C)
  //  return false;
  
  auto dvz = electron->par()->getVz()-h->par()->getVz();
  
  //if(dvz < cut_dvzmin || dvz > cut_dvzmax)
  //  return false;
}*/



void MakeHistogramsMC_noPID(){
  // Record start time
  auto start = std::chrono::high_resolution_clock::now();
  //macro to create a single column
#define leaf(name) double name=0;  tree->Branch(#name,&name,#name+(TString)"/D");
  TTree* tree = new TTree("events","events");
  leaf(Q2);
  leaf(W);
  leaf(x);
  leaf(nu);
  leaf(y);
  
  leaf(e_p);
  leaf(e_px);
  leaf(e_py);
  leaf(e_pz);
  leaf(e_ph);
  leaf(e_th);
  leaf(e_E);
  leaf(e_vx);
  leaf(e_vy);
  leaf(e_vz);

  leaf(invmass);
  leaf(missmass);
  //0=K+pi-
  //1=K-pi+
  //2=K-pi+p
  //3=K+pi-p
  leaf(topo);

  //count the number of each of these types of particles
  leaf(n_pim);
  leaf(n_pip);
  leaf(n_p);
  leaf(n_Kp);
  leaf(n_Km);

//macro to create one column for each particle.  
#define leaf3(name) leaf(K_##name); leaf(pi_##name); leaf(prot_##name);

  leaf3(p);
  leaf3(px);
  leaf3(py);
  leaf3(pz);
  leaf3(E);
  leaf3(ph);
  leaf3(th);
  leaf3(vx);
  leaf3(vy);
  leaf3(vz);

#define fillVec(v,pref) pref##_p=v.P(); pref##_px=v.Px(); pref##_py=v.Py(); pref##_pz=v.Pz(); pref##_E=v.E(); pref##_th=v.Theta(); pref##_ph=v.Phi();
#define blankVec(pref) pref##_p=0; pref##_px=0; pref##_py=0; pref##_pz=0; pref##_E=0; pref##_th=0; pref##_ph=0;

  //time
  //leaf(e_t);
  leaf(prot_dt);
  leaf(K_dt);
  leaf(pi_dt);

  leaf(evt_num); leaf(run_num);
  /////////////////////////////////////
  //ignore this just getting file name!
  
  TString outputFile;
  TChain input("hipo");
  int skipEvents=0;
  int maxevents=-1;
  bool useCuts = 1;
  bool FTonly=0;
  bool DConly=0;
  bool isMC = 0;
  TString qadbPath="";
  
  for(Int_t ii=1;ii<gApplication->Argc();ii++){
    TString opt=gApplication->Argv(ii);
    if((opt.Contains("--in="))){
      TString inputFile=opt(5,opt.Sizeof());
      cout << inputFile << endl;
      if(inputFile.EndsWith("/")){
        cout << inputFile << endl;
        inputFile=inputFile(0,inputFile.Sizeof()-1);
        
        TSystemDirectory dir(inputFile, inputFile);
        TList *files = dir.GetListOfFiles();
        for (int j =0;j<files->GetEntries();j++){
          TString name = files->At(j)->GetName();
          if (name == "." || name == "..")
            continue;
          input.Add(inputFile+"/" +name);
          cout << files->At(j)->GetName() <<endl;
        }
      } else
        input.Add(inputFile);
    } else if (opt.Contains("--out=")){
      outputFile = opt(6,opt.Sizeof());
    } else if (opt.Contains("--N=")){
      cout << (TString)opt(4,opt.Sizeof())<< endl;
      maxevents= ((TString)opt(4,opt.Sizeof())).Atoi();
      cout << "maxevents set to " << maxevents << endl;
    } else if(opt.Contains("--skipEvents=")){
      skipEvents = ((TString)opt(13,opt.Sizeof())).Atoi();
      cout << "skipping the first " << skipEvents << " events" << endl;
    } else if (opt.Contains("--skipCuts")){
      useCuts=0;
    } else if (opt.EqualTo("--isMC")){  //not used
      isMC=1;
    } else if (opt.Contains("--FTonly")){ //only use electrons in the FT
      FTonly=1;
    } else if (opt.Contains("--DConly")){ //only use electrons in the DC
      DConly=1;
    } else if (opt.Contains("--debug")){
      debug=1;
    } else if (opt.Contains("--qadbPath=")){
      qadbPath=opt(11,opt.Sizeof());
      cout << "qadb path: "<< qadbPath <<endl;
    }
    
  }

  TLorentzVector had;

  //D0bar recon
  TH1* invmass_Kp_pim = new TH1D("Kp_pim_invmass", "pair mass K+pi-;m_{K^{+}#pi^{-}} [GeV];events",2000, 0, 10);
  TH1* invmass_Kp_pim_zoom = new TH1D("Kp_pim_invmass_zoom", "pair mass K+pi-;m_{K^{+}#pi^{-}} [GeV];events",500, 1.5, 2.0);
  TH1* missingmass_Kp_pim = new TH1D("Kp_pim_missmass", "missing mass K+pi-;m_{X}(ep,e'K^{+}#pi^{-}X) [GeV];events",2000, -10, 10);
  TH2* Kp_pim_2d = new TH2D("Kp_pim_2d", "invariant vs missing mass K+pi-;m_{K^{+}#pi^{-}} [GeV];m_{X}(ep,e'K^{+}#pi^{-}X) [GeV];events",2000,0, 10, 2000, -10, 10);
  TH2* Kp_pim_2d_zoom=new TH2D("Kp_pim_2d_zoom", "invariant vs missing mass K+pi-;m_{K^{+}#pi^{-}} [GeV];m_{X}(ep,e'K^{+}#pi^{-}X) [GeV];events", 500, 1.5, 2, 500, 2.0, 2.5);
  TH2* invmass_Kp_pim_vs_W = new TH2D("Kp_pim_invmass_vs_w", "pair mass K+pi-;m_{K^{+}#pi^{-}} [GeV];W [GeV]",500, 0, 10,100, 0, 10);


  //Lambda c recon
  TH1* invmass_Km_pip_p = new TH1D("Km_pip_p_invmass", "pair mass K-pi+p;m_{K^{-}#pi^{+}p} [GeV];events",2000, 0, 10);
  TH1* missingmass_Km_pip_p = new TH1D("Km_pip_p_missmass", "missing mass K-pi+p;m_{X}(ep,e'K^{-}#pi^{+}pX) [GeV];events",2000, -10, 10);
  TH2* Km_pip_p_2d = new TH2D("Km_pip_p_2d", "invariant vs missing mass K-pi+p;m_{K^{-}#pi^{+}p} [GeV];m_{X}(ep,e'K^{+}#pi^{-}p'X) [GeV];events",2000,0, 10, 2000, -10, 10);
  TH2* Km_pip_p_2d_zoom = new TH2D("Km_pip_p_2d_zoom", "invariant vs missing mass K-pi+p;m_{K^{-}#pi^{+}p} [GeV];m_{X}(ep,e'K^{+}#pi^{-}p'X) [GeV];events",500,2.0, 2.5, 500, 1.5, 2);
  TH2* invmass_Km_pip_p_vs_W = new TH2D("Km_pip_p_invmass_vs_w", "pair mass K+pi-;m_{K^{-}#pi^{+}p} [GeV];W [GeV]",2000, 0, 10,100, 0, 10);


  //Lambda recon
  TH1* invmass_pim_p = new TH1D("pim_p_invmass", "pair mass pi-p;m_{#pi{-}p} [GeV];events",1000, 0, 5);
  TH1* invmass_pim_p_zoom = new TH1D("pim_p_invmass_zoom", "pair mass pi-p;m_{#pi{-}p} [GeV];events",500, 1.05, 1.30);

  //D0 recon (not D0 bar)
  TH1* invmass_Km_pip = new TH1D("Km_pip_invmass", "pair mass K-pi+;m_{K^{-}#pi^{+}} [GeV];events",2000, 0, 10);
  TH1* missingmass_Km_pip = new TH1D("Km_pip_missmass", "missing mass K-pi+;m_{X}(ep,e'K^{-}#pi^{+}X) [GeV];events",2000, -10, 10);
  TH2* Km_pip_2d = new TH2D("Km_pip_2d", "invariant vs missing mass K-pi+;m_{K^{-}#pi^{+}} [GeV];m_{X}(ep,e'K^{-}#pi^{+}X) [GeV];events",2000,0, 10, 2000, -10, 10);
  TH2* Km_pip_2d_zoom=new TH2D("Km_pip_2d_zoom", "invariant vs missing mass K-pi+;m_{K^{-}#pi^{+}} [GeV];m_{X}(ep,e'K^{-}#pi^{+}X) [GeV];events", 500, 1.5, 2, 500, 2.0, 2.5);
  TH2* invmass_Km_pip_vs_W = new TH2D("Km_pip_invmass_vs_w", "pair mass K-pi+;m_{K^{-}#pi^{+}} [GeV];W [GeV]",2000, 0, 10,100, 0, 10);
  TH1* invmass_Km_pip_zoom = new TH1D("Km_pip_invmass_zoom", "pair mass K-pi+;m_{K^{-}#pi^{+}} [GeV];events",500, 1.5, 2.0);
  

  std::map<int, TH2*> singleParticle2d;
  std::map<int, TH1*> singleParticle_p;
  std::map<int, TH1*> singleParticle_theta;

  //std::map<int, TH1*> mcres_dp;
  //std::map<int, TH1*> mcres_dth;
  //std::map<int, TH1*> mcres_dph;
  TH1* hadron_pres = new TH1D("hadron p res", "hadron #Delta p/p;#Delta p/p [%];events",200, -50, 50);
  TH1* hadron_thetares = new TH1D("hadron theta res", "hadron #Delta#Theta;#Delta#Theta [mrad];events",200, -50, 50);
  TH1* hadron_phires = new TH1D("hadron phi res", "hadron #Delta#Phi;#sin#Theta #Delta#Phi [mrad];events",200, -50, 50);

  TH2* hadron_pres_vs_p = new TH2D("hadron p res vs p", "hadron #Delta p/p;p_{truth} [GeV];#Delta p/p [%]",200, 0, 10, 200, -50, 50);
  TH2* hadron_thetares_vs_p = new TH2D("hadron p res vs p", "hadron #Delta#Theta;p_{truth} [GeV];#Delta#Theta [mrad]",200, 0, 10, 200, -50, 50);
  TH2* hadron_phires_vs_p = new TH2D("hadron phi res vs p", "hadron #Delta#Phi vs p;p_{truth} [GeV];#sin#Theta #Delta#Phi [mrad];events",200, 0,10, 200, -50, 50);


  // Loop through the list using a range-based for loop
  for (const int pid : {11,211,2212,321,-11,-211,-321}){
    singleParticle2d[pid]=new TH2D(("pvstheta"+std::to_string(pid)).c_str(),("p vs theta "+std::to_string(pid)+";P [GeV];#Theta [deg]").c_str(), 100, 0, 10, 180, 0, 180);
    singleParticle_theta[pid]=new TH1D(("theta"+std::to_string(pid)).c_str(),("theta "+std::to_string(pid)+";#Theta [deg];events").c_str(), 180, 0, 180);
    singleParticle_p[pid]=new TH1D(("p"+std::to_string(pid)).c_str(),("p "+std::to_string(pid)+";P [GeV];events").c_str(), 100, 0, 10);
  }

  // leptonic variables
  TH1* h_nu = new TH1D("nu", "#nu;#nu [GeV];events", 100, 8, 11);
  TH1* h_nu_mc = new TH1D("nu MC", "#nu (MC):#nu [GeV];evetns",100, 8, 11);
  TH1* h_Q2 = new TH1D("Q2", "Q^{2};Q^{2} [GeV^2];events", 100,0, 11);
  TH1* h_W = new TH1D("W", "W;W [GeV];events", 100, 0, 11);
  TH2* h_Q2_nu = new TH2D("Q2_nu", "Q^{2} vs #nu;Q^{2} [GeV^2];#nu [GeV]",100, 0, 11, 100, 0, 11);
  
  clas12::clas12databases *dbc12;
  
  clas12::clas12databases::SetCCDBRemoteConnection();
  if(!isMC){
    clas12::clas12databases::SetRCDBRemoteConnection();
    /*if(!qadbPath.EqualTo("")){
     cout << "setting qadb connection" <<endl;
     clas12::clas12databases::SetQADBConnection((const std::string)qadbPath);
     cout << "successfully connected to qadb" <<endl;
     }*/
  }
  
  dbc12 = new clas12::clas12databases();
  
  
  auto files=input.GetListOfFiles();
  
  //some particles
  //static auto db=TDatabasePDG::Instance();
  
  
  //double recoilMass = 0;
  //by default the target mass is that of the proton.
  TLorentzVector target(0,0,0,db->GetParticle(2212)->Mass());
  TLorentzVector el(0,0,0,db->GetParticle(11)->Mass());
  TLorentzVector cm;
  cout << "proton mass is " << db->GetParticle(2212)->Mass() << endl;  
  
  gBenchmark->Start("timer");
  int count=0;
  
  
  for(Int_t filenum=0;filenum<files->GetEntries();filenum++){
    //create the event reader
    clas12reader c12(files->At(filenum)->GetTitle(),{0});
    
    c12.connectDataBases(dbc12);
    clas12::ccdb_reader *ccdb = c12.ccdb();
    
    const TableOfDoubles_t& electron_sf = ccdb->requestTableDoubles("/calibration/eb/electron_sf");
    
    
    for(int row = 0; row<=6; row++){
      sf1[row] = electron_sf[row][3];
      sf2[row] = electron_sf[row][4];
      sf3[row] = electron_sf[row][5];
      sf4[row] = electron_sf[row][6];
      sfs1[row] = electron_sf[row][7];
      sfs2[row] = electron_sf[row][8];
      sfs3[row] = electron_sf[row][9];
      sfs4[row] = electron_sf[row][10];
    }
    ccdb->close();
    
    
    double E;
    if(!isMC){
      
      
      
      if(!qadbPath.EqualTo("")){
        c12.db()->qadb_requireGolden(true);
        c12.applyQA();
      }
      clas12::rcdb_reader *rcdb = c12.rcdb();
      
      auto& current = rcdb->current();
      E = current.beam_energy/1000;
      
      
      
      double torus_curr = current.torus_current;
      double torus_scale = current.torus_scale;
      if(torus_scale <0){
        torus=INBEND;
        cut_evzmin = cut_evzmin_inb;
        cut_evzmax = cut_evzmax_inb;
      } else {
        torus=OUTBEND;
        cut_evzmin = cut_evzmin_outb;
        cut_evzmax = cut_evzmax_outb;
        cout << cut_evzmin << " " << cut_evzmax <<  endl;
      }
      
      
      cout << "torus current, scale = " << torus_curr << " " << torus_scale << endl;
      TString targetName(rcdb->current().target);
      cout << "target is " << targetName << "." << endl;
      std::cout << "beam energy is " << E << std::endl;
      rcdb->close();
      
    } else{
      E=0; //get the beam energy later
    }
    
    
    TLorentzVector beam(0,0,E,E);
    
    c12.useFTBased();

    while(c12.next()==true){
      if(isMC && E==0){
        auto dict = c12.getDictionary();
        if(dict.hasSchema("MC::Event")){
          auto schema = dict.getSchema("MC::Event");
          hipo::bank bank(schema);
          c12.getStructure(&bank);
          E = bank.getFloat(schema.getEntryOrder("ebeam"),0);
          beam = {0, 0, E, E};
        }
        if(dict.hasSchema("MC::Header")){
          auto schema = dict.getSchema("MC::Header");
          hipo::bank bank(schema);
          c12.getStructure(&bank);
          //helicity = bank.getFloat(schema.getEntryOrder("helicity"),0);
        }
      } else {
        auto dict = c12.getDictionary();
        if(dict.hasSchema("RUN::config")){
          auto schema = dict.getSchema("RUN::config");
          hipo::bank bank(schema);
          c12.getStructure(&bank);
          evt_num = bank.getInt(schema.getEntryOrder("event"),0);
          run_num = bank.getInt(schema.getEntryOrder("run"),0);
        }
        
      }
      //cout << "new event" << endl;
      //       c12.addARegionCDet();
      if(!((count-skipEvents) %10000) && count >= skipEvents)
        cout << (count-skipEvents) <<"events processed; " << maxevents << "requested"<< endl;
      //count ++;
      
      if(count-skipEvents > maxevents && maxevents >0)
        break;
      if(count < skipEvents){
        count ++;
        continue;
      }
      count ++;
      
      
      // get particles by type
      
      auto electrons=c12.getByID(11);
      
      TLorentzVector el(0,0,0,db->GetParticle(11)->Mass());
      

      
      auto parts=c12.getDetParticles();

      auto mcparts = c12.mcparts();
      //cout << "number of rows in mc parts" << mcparts->getRows() << endl;
      for(short i=0; i<mcparts->getRows(); i++){
	mcparts->setEntry(i);
	//cout << "mc particle #"<< i << ":  pid is "<< mcparts[i].getPid() <<endl ;
	if(mcparts->getPid()==11){
	  h_nu_mc->Fill(E-hypot(mcparts->getP(),0.000511));
	  break;
	}
      }
      //cout << "end electron loop"<<endl;
      int ntracks=0;

      //cout << "counting tracks"<<endl;                                                                                                  
      for(int kk = 0; kk<parts.size();kk++){
	auto part = parts[kk];
	int pid = part->getPid();
	if(db->GetParticle(pid) == NULL || db->GetParticle(pid)->Charge() == 0)
	  continue;
	ntracks++;
      }

      //if (ntracks>5)
      //  continue;
      int nelectrons = electrons.size();
      int electrons_passCuts = 0;
      vector<int> matchedMCindices = {};
      //if(debug) cout <<"CHECK 0: " << nelectrons << " electrons" << endl;
      for(int i=0; i<nelectrons; i++){
        
        bool Fdet_ok = 1;
        bool FT_ok = 1;
        if(debug) cout << "starting electron loop" << endl;
        int sector = electrons[i]->getSector();
        if(debug) cout << "electron in sector" << sector <<endl;
        

        if(debug) cout << "checking electron dc fid"<<endl;
        if(!dcOK(electrons[i],torus) || !pcalOK(electrons[i]))
          Fdet_ok = 0;
        if(debug) cout << "electron passes dc fid cuts" <<endl;
        
        SetLorentzVector(el,electrons[i], 0.000511);
        e_p = el.P();
        e_px = el.Px();
        e_py = el.Py();
        e_pz = el.Pz();
        e_th = el.Theta();
        e_ph = el.Phi();
	e_E = el.E();
        double ecal = electrons[i]->getDetEnergy();
        
        double e_ecalfrac = ecal/e_p;
        double e_pcal = electrons[i]->cal(PCAL)->getEnergy();
        if(ecal == 0 || e_p==0 || e_pcal == 0)
          Fdet_ok=0;
        double e_PCALx = electrons[i]->cal(PCAL)->getX();
        double e_PCALy = electrons[i]->cal(PCAL)->getY();
        
        
        //ecal fraction
        
        //cout << "sampling fraction check" << endl;
        double ecalfrac_mu = sf1[sector]*(sf2[sector]+sf3[sector]/ecal+sf4[sector]/pow(ecal,2));
        double ecalfrac_sigma = sfs1[sector]*(sfs2[sector]+sfs3[sector]/ecal+sfs4[sector]/pow(ecal,2));
        double nsig = 3.5;
        if(e_ecalfrac>ecalfrac_mu+nsig*ecalfrac_sigma || e_ecalfrac<ecalfrac_mu-nsig*ecalfrac_sigma)
          Fdet_ok=0;
        
        //cout << "passed sampling fraction check" << endl;
        
        if(e_pcal<cut_pcalmin)
          Fdet_ok=0;
        double e_ecalin = electrons[i]->cal(ECIN)->getEnergy();
        double e_ecalout = electrons[i]->cal(ECOUT)->getEnergy();
        //cout << e_ecalin << " " << e_ecalout << " " << e_p << " " << e_pcal << endl;
        //->cal(ECOUT)->getEnergy();
        
        //further cut to remove pions
        if(e_p> 4.5 && e_ecalin/e_p < 0.2 - e_pcal/e_p)
          Fdet_ok=0;
        //if(debug) cout << "/electron"<<endl;
        
	e_vx = electrons[i]->par()->getVx();
	e_vy = electrons[i]->par()->getVy();
        e_vz = electrons[i]->par()->getVz();
        
        if(useCuts && (e_vz<cut_evzmin || e_vz> cut_evzmax))
          Fdet_ok=0;
        
        double nHTCC = electrons[i]->che(HTCC)->getNphe();
        if(nHTCC<=cut_HTCCmin)
          Fdet_ok=0;
        
        if (electrons[i]->ft(0)==NULL)
          FT_ok=0;
        if (e_p<0.500 || e_th > 4.5*3.14159/180 || e_th<2.5*3.14159/180)
          FT_ok=0;
	int status = electrons[i]->getStatus();
	if (!(abs(status)>=1000 && abs(status)<2000))
          FT_ok=0;
        if (FT_ok)
          if(debug) cout << "electron from FT" << endl;
        if (!(Fdet_ok || FT_ok))
          continue;
        if (!Fdet_ok && DConly)
          continue;
        if (!FT_ok && FTonly)
          continue;
        
        //done with electron id cuts
        //cout << "check 2" <<endl;
        
        //if(debug) cout << "now for electron kinematics cuts"<<endl;
        // now for electron kinematics cuts
        Q2 = -(beam-el)*(beam-el);
        W = (target+beam-el).M();
        x = Q2/(2*target.M()*(beam.E()-el.E()));
        nu = (beam.E()-el.E());
        y = nu/E;
        

        
        //quasi-real cuts
        double mrho=.775;
        if (Q2>mrho*mrho)
          continue;
        if (nu<8) 
          continue;
        
	h_Q2->Fill(Q2);
	h_nu->Fill(nu);
	h_Q2_nu->Fill(Q2, nu);
	h_W->Fill(W);

        cm = beam+target-el;
        
        TLorentzVector target_cm;
        toCM(cm, target,target_cm);

        
        //double nHTCC = electrons[i]->che(HTCC)->getNphe();
        //if(nHTCC<=cut_HTCCmin)
        //  continue;
        
        auto parts=c12.getDetParticles();
        
        
        vector<int> accepted_indices;
        map<int, double> dtime_map;
	map<int, double> vx_map;
	map<int, double> vy_map;
	map<int, double> vz_map;
	n_pip=0; n_pim=0;n_Kp=0;n_Km=0; n_p=0;

        for(int j =0; j<parts.size();j++){
          if(debug) cout << "hadrons loop"<< endl;
          auto h = parts[j];
          
          auto h_pid = h->getPid();
          if(h_pid != 211 && h_pid != -211 && h_pid != 2212 && h_pid != -2212 && h_pid != 321 && h_pid != -321)
            continue;
          
          
          double dtime = electrons[i]->getTime()-h->getTime();
          bool useDCcuts=0;
          //if(debug) cout << "dc"<<endl;
          if(useDCcuts && !dcOK(h,torus))
            continue;
          double h_p = had.P();  
          //if(debug) cout << "pip" << endl;
	  
          double mass = db->GetParticle(h_pid)->Mass();
          SetLorentzVector(had,h, mass);
          
          double c = 29.9792458; //cm/ns
          auto dtime_corr =dtime-electrons[i]->getPath()/c+h->getPath()/(had.Beta()*c);
          if(abs(dtime_corr) > cut_dtime_corr)
            continue;
          dtime_map[j]=dtime_corr;
          
          
          auto h_chi2pid = h->par()->getChi2Pid();
          /*
          //scaling constant for pid cut
          double C = .88*(h_pid == 211)+.93*(h_pid == -211)+ 1.30*(abs(h_pid)==2212)+1.0*(abs(h_pid) != 211 && abs(h_pid)!=2212);
          if(abs(h_chi2pid) >3*C)
            continue;
          
          //tighter cut for pions
          if(abs(h_pid) == 211 && h_p>2.44 && h_chi2pid >(0.00869+14.98587*exp(-h_p/1.18236)+1.81751*exp(-h_p/4.86394))*C)
	  continue;*/
          vx_map[j]=h->par()->getVx();
	  vy_map[j]=h->par()->getVy();
	  vz_map[j]=h->par()->getVz();
          //auto dvz = electrons[i]->par()->getVz()-h->par()->getVz();
          
          //if(dvz < cut_dvzmin || dvz > cut_dvzmax)
          //  continue;
          if (debug) cout << "accepted hadron" << h_pid << endl;
          accepted_indices.push_back(j);
	  
	  if (h_pid==211) n_pip+=1;
	  if (h_pid==-211) n_pim+=1;
	  if (h_pid==321) n_Kp+=1;
          if (h_pid==-321) n_Km+=1;
	  if (h_pid==2212) n_p+=1;

	  /*
	  if (abs(h_pid)==211 or abs(h_pid)==321 or abs(h_pid)==11 or h_pid==2212){
	    singleParticle2d[h_pid]->Fill(h_p, had.Theta()*180/3.14159);
	    singleParticle_theta[h_pid]->Fill(had.Theta()*180/3.14159);
            singleParticle_p[h_pid]->Fill(had.P());
	  }
	  */
	  if (1) {
	    int bestmcindex=-1;
	    double bestdist=0.01;
	    for(int mcindex=0; mcindex<mcparts->getRows(); mcindex++){
	      mcparts->setEntry(mcindex);
	      //cout << "mc particle #" << mcindex <<"  pid=" << mcparts->getPid() << "  theta=" << mcparts->getTheta() <<endl;
	      double dist=hypot(h->getTheta()-mcparts->getTheta(), sin(h->getTheta())*angle(h->getPhi()-mcparts->getPhi()));
	      if (dist <bestdist && mcparts->getPid()!=11){
		bestmcindex=mcindex;
		bestdist=dist;
	      }
	    }
	    if (bestmcindex !=-1){
	      mcparts->setEntry(bestmcindex);
	      //cout <<"recon pid:" <<h_pid<< "  matched to true pid:" << mcparts->getPid() << endl;
	      hadron_pres->Fill(100*(h->getP()-mcparts->getP())/mcparts->getP());
	      hadron_thetares->Fill(1000*(h->getTheta()-mcparts->getTheta()));
	      hadron_pres_vs_p->Fill(mcparts->getP(), 100*(h->getP()-mcparts->getP())/mcparts->getP());
              hadron_thetares_vs_p->Fill(mcparts->getP(), 1000*(h->getTheta()-mcparts->getTheta()));
	      hadron_phires->Fill(1000*angle(h->getPhi()-mcparts->getPhi())*sin(mcparts->getTheta()));
	      hadron_phires_vs_p->Fill(mcparts->getP(),1000*angle(h->getPhi()-mcparts->getPhi())*sin(mcparts->getTheta()));
	    }
	  }
	}
      
        double Kmass=db->GetParticle(321)->Mass();
        double pimass=db->GetParticle(211)->Mass();
        double pmass=db->GetParticle(2212)->Mass();
	
	// D0 bar
        for(int j = 0; j< accepted_indices.size();j++){
          auto part = parts[accepted_indices[j]];
          if (part->getPid()>0)//positive charge (K+ candidate)
          {
            if (debug) cout << "accepted K+" << endl;
            for(int k = 0; k< accepted_indices.size();k++){
              auto part2 = parts[accepted_indices[k]];
              if(part2->getPid()<0){//negative charge (pi- candidate)   
                if (debug) cout << "accepted pi-" << endl;
                TLorentzVector Kp;
                TLorentzVector pim;
                SetLorentzVector(Kp,part, Kmass);
                SetLorentzVector(pim,part2,pimass);
                invmass = (Kp+pim).M();
                if(debug) cout << "pair mass "  << invmass <<  endl;
                invmass_Kp_pim->Fill(invmass);
		invmass_Kp_pim_zoom->Fill(invmass);
		invmass_Kp_pim_vs_W->Fill(invmass, W);
                missmass = (target+beam-el-Kp-pim).M();
                missingmass_Kp_pim->Fill(missmass);
                Kp_pim_2d->Fill(invmass, missmass);
                Kp_pim_2d_zoom->Fill(invmass, missmass);
		topo=0;
		fillVec(Kp,K);
		fillVec(pim,pi);
		blankVec(prot);
		
		
		prot_dt=0; K_dt=dtime_map[accepted_indices[j]];
		pi_dt=dtime_map[accepted_indices[k]];
		
		K_vx=vx_map[accepted_indices[j]];
		K_vy=vy_map[accepted_indices[j]];
		K_vz=vz_map[accepted_indices[j]];
		pi_vx=vx_map[accepted_indices[k]];
		pi_vy=vy_map[accepted_indices[k]];
		pi_vz=vz_map[accepted_indices[k]];
		prot_vx=0; prot_vy=0; prot_vz=0;
		
		tree->Fill();
              }
            }
          }
          
        }

	// Lambda c
        for(int j = 0; j< accepted_indices.size();j++){
          auto part = parts[accepted_indices[j]];
          if (part->getPid()<0) //K- candidate
          {
            if (debug) cout << "accepted K-" << endl;
            for(int k = 0; k< accepted_indices.size();k++){
              auto part2 = parts[accepted_indices[k]];
              if(part2->getPid()>0){ //pi+ candidate
                if (debug) cout << "accepted pi+" << endl;
                for(int l = 0; l< accepted_indices.size();l++){
                  auto part3 = parts[accepted_indices[l]];
                  if(part3->getPid()>0){ //p candidate
                    if (debug) cout << "accepted p" << endl;
                    TLorentzVector Km;
                    TLorentzVector pip;
                    TLorentzVector prot;
                    SetLorentzVector(Km,part, Kmass);
                    SetLorentzVector(pip,part2,pimass);
                    SetLorentzVector(prot,part3,pmass);
                    invmass = (Km+pip+prot).M();
                    if (debug) cout << "pair mass "  << invmass <<  endl;
                    invmass_Km_pip_p->Fill(invmass);
		    invmass_Km_pip_p_vs_W->Fill(invmass, W);
                    missmass = (target+beam-el-Km-pip-prot).M();
                    missingmass_Km_pip_p->Fill(missmass);
                    Km_pip_p_2d->Fill(invmass, missmass);
                    Km_pip_p_2d_zoom->Fill(invmass, missmass);
		    fillVec(Km,K);
		    fillVec(pip,pi);
		    fillVec(prot,prot);
		    prot_dt=dtime_map[accepted_indices[l]];
		    K_dt=dtime_map[accepted_indices[j]];
		    pi_dt=dtime_map[accepted_indices[k]];

		    K_vx=vx_map[accepted_indices[j]];
		    K_vy=vy_map[accepted_indices[j]];
		    K_vz=vz_map[accepted_indices[j]];
		    pi_vx=vx_map[accepted_indices[k]];
		    pi_vy=vy_map[accepted_indices[k]];
		    pi_vz=vz_map[accepted_indices[k]];
		    prot_vx=vx_map[accepted_indices[l]];
                    prot_vy=vy_map[accepted_indices[l]];
                    prot_vz=vz_map[accepted_indices[l]];
		    
		    topo=2;
		    tree->Fill();
                  }
                }
              }
            }
          }
          
        }
      }    
    }
  }
  gBenchmark->Stop("timer");
  gBenchmark->Print("timer");
  //hm2gCut->SetLineColor(2);
  //hm2gCut->DrawCopy("same");
  TFile *f = new TFile(outputFile,"RECREATE");
  tree->Write();
  //D0 bar
  missingmass_Kp_pim->Write();
  invmass_Kp_pim->Write();
  invmass_Kp_pim_vs_W->Write();
  invmass_Kp_pim_zoom->Write();
  Kp_pim_2d->Write();
  Kp_pim_2d_zoom->Write();
  //lambda c
  missingmass_Km_pip_p->Write();
  invmass_Km_pip_p->Write();
  invmass_Km_pip_p_vs_W->Write();
  Km_pip_p_2d->Write();
  Km_pip_p_2d_zoom->Write();
  //lambda
  invmass_pim_p->Write();
  invmass_pim_p_zoom->Write();
  //D0
  missingmass_Km_pip->Write();
  invmass_Km_pip_zoom->Write();
  invmass_Km_pip->Write();
  invmass_Km_pip_vs_W->Write();
  Km_pip_2d->Write();
  Km_pip_2d_zoom->Write();
  invmass_Km_pip_zoom->Write();
  for (const int pid : {11,211,2212,321,-11,-211,-321}){
    singleParticle2d[pid]->Write();
    singleParticle_p[pid]->Write();
    singleParticle_theta[pid]->Write();
  }

  //lepton variables
  h_nu->Write();
  h_W->Write();
  h_Q2->Write();
  h_Q2_nu->Write();

  
  h_nu_mc->Write();

  hadron_pres->Write();
  hadron_thetares->Write();
  hadron_phires->Write();

  hadron_pres_vs_p->Write();
  hadron_thetares_vs_p->Write();
  hadron_phires_vs_p->Write();
  f->Close();
  
  
  
  
  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish - start;
  std::cout << "Elapsed time: " << elapsed.count()<< " events = "<<count<< " s\n";
  
}
