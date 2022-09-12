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
double cut_dtime_corr = 100;//0.3;

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


void MakeHistograms(){
  // Record start time
  auto start = std::chrono::high_resolution_clock::now();
  
  
  /////////////////////////////////////
  //ignore this just getting file name!
  int torus;
  TString outputFile;
  TChain input("hipo");
  int skipEvents=0;
  int maxevents=-1;
  bool useCuts = 1;
  bool FTonly=0;
  bool DConly=0;
  bool isMC = 0;
  TString qadbPath="";

  TLorentzVector had;
  TH1* invmass_Kp_pim = new TH1D("Kp_pim_invmass", "pair mass K+pi-;m_{K^{+}#pi^{-}} [GeV];events",500, 0, 10);
  TH1* missingmass_Kp_pim = new TH1D("Kp_pim_missmass", "missing mass K+pi-;m_{X}(ep,e'K^{+}#pi^{-}X) [GeV];events",500, -10, 10);
  TH2* Kp_pim_2d = new TH2D("Kp_pim_2d", "invariant vs missing mass K+pi-;m_{K^{+}#pi^{-}} [GeV];m_{X}(ep,e'K^{+}#pi^{-}X) [GeV];events",500,0, 10, 500, -10, 10);
  TH2* Kp_pim_2d_zoom=new TH2D("Kp_pim_2d_zoom", "invariant vs missing mass K+pi-;m_{K^{+}#pi^{-}} [GeV];m_{X}(ep,e'K^{+}#pi^{-}X) [GeV];events", 500, 1.5, 2, 500, 2.0, 2.5);
  
  TH1* invmass_Km_pip_p = new TH1D("Km_pip_p_invmass", "pair mass K-pi+p;m_{K^{-}#pi^{+}p} [GeV];events",500, 0, 10);
  TH1* missingmass_Km_pip_p = new TH1D("Km_pip_p_missmass", "missing mass K-pi+p;m_{X}(ep,e'K^{-}#pi^{+}pX) [GeV];events",500, -10, 10);
  TH2* Km_pip_p_2d = new TH2D("Km_pip_p_2d", "invariant vs missing mass K-pi+p;m_{K^{-}#pi^{+}p} [GeV];m_{X}(ep,e'K^{+}#pi^{-}p'X) [GeV];events",500,0, 10, 500, -10, 10);
  TH2* Km_pip_p_2d_zoom = new TH2D("Km_pip_p_2d_zoom", "invariant vs missing mass K-pi+p;m_{K^{-}#pi^{+}p} [GeV];m_{X}(ep,e'K^{+}#pi^{-}p'X) [GeV];events",500,2.0, 2.5, 500, 1.5, 2);

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
  auto db=TDatabasePDG::Instance();
  
  
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
    
    //  clas12reader c12(files->At(i)->GetTitle(),{0});//add tags {tag1,tag2,tag3,...}
    
    //Add some event Pid based selections
    //////////c12.AddAtLeastPid(211,1); //at least 1 pi+
    //c12.addExactPid(11,1);    //exactly 1 electron
    //c12.addExactPid(211,1);    //exactly 1 pi+
    //c12.addExactPid(-211,1);    //exactly 1 pi-
    //c12.addExactPid(2212,1);    //exactly 1 proton
    //c12.addExactPid(22,2);    //exactly 2 gamma
    //////c12.addZeroOfRestPid();  //nothing else
    //////c12.useFTBased(); //and use the Pids from RECFT
    
    //can also access the integrated current at this point
    //c12.scalerReader();//must call this first
    //c12.getRunBeamCharge();
    
    //int max = 100000;
    
    
    // for data, this step can be done.
    // with MC, there are events with
    // electrons that fail recon,
    // so these must not be excluded
    if(!isMC)
      c12.addAtLeastPid(11,1);
    
    // for event mixing:  the cm frame fourvector two events ago
    TLorentzVector cm_mix = {0,0,0,0};
    // for event mixing:  the cm frame fourvector one events ago
    TLorentzVector cm_mix_tmp = {0,0,0,0};
    // for event mixing: the trigger hadron in the previous entry.
    TLorentzVector had_mix = {0,0,0,0};
    
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
          //evt_num = bank.getInt(schema.getEntryOrder("event"),0);
	  //run_num = bank.getInt(schema.getEntryOrder("run"),0);
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
      //can get an estimate of the beam current to this event
      //c12.getCurrApproxCharge();//if called c12.scalerReader();
      
      //c12.event()->getStartTime();
      
      
      //Loop over all particles to see how to access detector info.
      /*
       for(auto& p : c12.getDetParticles()){
       //  get predefined selected information
       p->getTime();
       p->getDetEnergy();
       p->getDeltaEnergy();
       
       //check trigger bits
       //	 if(c12.checkTriggerBit(25)) cout<<"MesonExTrigger"<<endl;
       //	 else cout<<"NOT"<<endl;
       
       // get any detector information (if exists for this particle)
       // there should be a get function for any entry in the bank
       switch(p->getRegion()) {
       case FD :
       p->cal(PCAL)->getEnergy();
       p->cal(ECIN)->getEnergy();
       p->cal(ECOUT)->getEnergy();
       p->sci(FTOF1A)->getEnergy();
       p->sci(FTOF1B)->getEnergy();
       p->sci(FTOF2)->getEnergy();
       p->trk(DC)->getSector();
       p->che(HTCC)->getNphe();
       p->che(LTCC)->getNphe();
       //trajectories
       p->traj(LTCC)->getX();
       // p->traj(DC,DC1)->getCx();; //First layer of DC, hipo4
       break;
       case FT :
       p->ft(FTCAL)->getEnergy();
       p->ft(FTHODO)->getEnergy();
       break;
       case CD:
       p->sci(CTOF)->getEnergy();
       p->sci(CND)->getEnergy();
       break;
       }
       //   covariance matrix (comment in to see!)
       // p->covmat()->print();
       p->cmat();
       }*/
      
      // get particles by type
      
      auto electrons=c12.getByID(11);
      
      if(c12.helonline() != NULL){
        //cout << "helicity is not null" << endl;
        //helicity = c12.helonline()->getHelicity();
        //helicity = c12.event()->getHelicity();
        //cout << "helicity is " << helicity << endl;
      }
      else{
        //cout << "helicity is NULL"<<endl;
      }
      //if(debug) cout << "helicity" << endl;
      TLorentzVector el(0,0,0,db->GetParticle(11)->Mass());
      
      /*mcpar_ptr mcparts;
      if(isMC){
        mcparts = c12.mcparts();
        nhtracks_truth =0;
        //skip the first two entries; these are the initial state particles
        //the intermediate state particles are removed by the pid cuts
        for(int kk = 2; kk<mcparts->getRows();kk++){
          int pid = mcparts->getPid(kk);
          if(debug) cout << pid << " ";
          if(abs(pid)==211 || abs(pid) == 321 || abs(pid)==2212)
            nhtracks_truth++;
        }
        if(debug) cout << endl;
      }
      */
      
      auto parts=c12.getDetParticles();
      /*if(electrons.size() == 2){
       if(debug) cout << "electrons" << endl;
       if(debug) cout << electrons[0]->par()->getP() <<endl;
       if(debug) cout << electrons[1]->par()->getP() <<endl;
       }*/
      int nelectrons = electrons.size();
      int electrons_passCuts = 0;
      vector<int> matchedMCindices = {};
      //if(debug) cout <<"CHECK 0: " << nelectrons << " electrons" << endl;
      for(int i=0; i<nelectrons; i++){

	bool Fdet_ok = 1;
	bool FT_ok = 1;
        if(debug) cout << "starting electron loop" << endl;
        //if(debug) cout << "CHECK 0.5" << endl;
        //if(electrons.size()>1) continue;
        int sector = electrons[i]->getSector();
        if(debug) cout << "electron in sector" << sector <<endl;
        
        //e_DC1x=electrons[i]->traj(DC,DC1)->getX();
        //e_DC1y=electrons[i]->traj(DC,DC1)->getY();
        //e_DC2x=electrons[i]->traj(DC,DC3)->getX();
        //e_DC2y=electrons[i]->traj(DC,DC3)->getY();
        //e_DC3x=electrons[i]->traj(DC,DC6)->getX();
        //e_DC3y=electrons[i]->traj(DC,DC6)->getY();
        if(debug) cout << "checking electron dc fid"<<endl;
        if(!dcOK(electrons[i],torus) || !pcalOK(electrons[i]))
          Fdet_ok = 0;
        if(debug) cout << "electron passes dc fid cuts" <<endl;
        //if(debug) cout << "dc and pcal ok"<<endl;
        //if(!dcok)
        //  continue;
        //if(debug) cout << "electron" << endl;
        
        //the electron mass is a myth.  I could set this to zero and nothing would change.
        //if(debug) cout << "CHECK 1" <<endl;
        SetLorentzVector(el,electrons[i], 0.000511);
        double e_p = el.P();
        double e_px = el.Px();
        double e_py = el.Py();
        double e_pz = el.Pz();
        //if(debug) cout<< "e_p: " << e_p <<" "<< el.P() << endl;
        double e_th = el.Theta();
        double e_ph = el.Phi();
        //if(theta*180/3.14159265 <=7)
        //continue;
        //if(p<0.01*E)
        //continue;
        
        //if(debug) cout << "CHECK 2"<<endl;
        double ecal = electrons[i]->getDetEnergy();
        
        double e_ecalfrac = ecal/e_p;
        double e_pcal = electrons[i]->cal(PCAL)->getEnergy();
        if(ecal == 0 || e_p==0 || e_pcal == 0)
          Fdet_ok=0;
        double e_PCALx = electrons[i]->cal(PCAL)->getX();
        double e_PCALy = electrons[i]->cal(PCAL)->getY();
        
        //double sf1 = 0, sf2 = 0, sf3 = 0, sf4 = 0;
        //double sfs1 =0, sfs2=0, sfs3 = 0, sfs4 = 0;
        
        if(debug) cout << sector << " " << sf1[sector] << " " <<sf2[sector] << " " <<sf3[sector] <<endl;
        
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
        
        double e_vz = electrons[i]->par()->getVz();
        
        if(useCuts && (e_vz<cut_evzmin || e_vz> cut_evzmax))
          Fdet_ok=0;
        
	double nHTCC = electrons[i]->che(HTCC)->getNphe();
        if(nHTCC<=cut_HTCCmin)
          Fdet_ok=0;
	
	if (electrons[i]->ft(0)==NULL)
	  FT_ok=0;
	if (e_p<0.500 || e_th > 4.5*3.14159/180 || e_th<2.5*3.14159/180)
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
        double Q2 = -(beam-el)*(beam-el);
        double W = (target+beam-el).M();
        double x = Q2/(2*target.M()*(beam.E()-el.E()));
        double nu = (beam.E()-el.E());
        double y = nu/E;
        
	// SIDIS cuts
        //if(Q2<cut_Q2min)
        //  continue;
        //if(W<cut_Wmin)
        //  continue;
        //if(y > cut_ymax)
        //  continue;
	
	//quasi-real cuts
        double mrho=.775;
	if (Q2>mrho*mrho)
	  continue;
	if (y < .70) //set this some margin away from the limit of 0.87
	  continue;

        cm = beam+target-el;
        
        TLorentzVector target_cm;
        toCM(cm, target,target_cm);
        //check that the cm outgoing electron has phi=0
        //TLorentzVector tmp;
        //toCM(cm, el,tmp);
        //cout << tmp.Phi() << endl;
        //toCM(cm, beam-el, tmp);
        //cout << "  " << tmp.Theta() << endl;
        
        //auto parts=c12.getDetParticles();
        
        //double nHTCC = electrons[i]->che(HTCC)->getNphe();
        //if(nHTCC<=cut_HTCCmin)
        //  continue;
        
        
        //if(debug) cout << "cm:"<<endl;
        //cm.Print();
        //npip = 0;
        //npim = 0;
        //npp = 0;
        //npm = 0;
        //nKp = 0;
        //nKm = 0;
        //nh = 0;

        
        
                
	auto parts=c12.getDetParticles();        
        
        //bool found_leader = 0;
        
        vector<int> accepted_indices;


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
	  
	  //if(debug) cout << "pip" << endl;
	  double mass = db->GetParticle(h_pid)->Mass();
	  SetLorentzVector(had,h, mass);
	  double h_p = had.P();
	  double c = 29.9792458; //cm/ns
	  auto dtime_corr =dtime-electrons[i]->getPath()/c+h->getPath()/(had.Beta()*c);
	  if(abs(dtime_corr) > cut_dtime_corr)
	    continue;
	    
            
	  //h_pid = h->par()->getPid();
	  auto h_chi2pid = h->par()->getChi2Pid();
            
	  //scaling constant for pid cut
	  double C = .88*(h_pid == 211)+.93*(h_pid == -211)+ 1.30*(abs(h_pid)==2212)+1.0*(abs(h_pid) != 211 && abs(h_pid)!=2212);
	  if(abs(h_chi2pid) >3*C)
	    continue;
	  
	  //tighter cut for pions
	  if(abs(h_pid) == 211 && h_p>2.44 && h_chi2pid >(0.00869+14.98587*exp(-h_p/1.18236)+1.81751*exp(-h_p/4.86394))*C)
	    continue;
	  
	  auto dvz = electrons[i]->par()->getVz()-h->par()->getVz();
          
	  if(dvz < cut_dvzmin || dvz > cut_dvzmax)
	    continue;
	  if (debug) cout << "accepted hadron" << h_pid << endl;
	  accepted_indices.push_back(j);
        }
        double Kmass=db->GetParticle(321)->Mass();
	double pimass=db->GetParticle(211)->Mass();
	double pmass=db->GetParticle(2212)->Mass();
	for(int j = 0; j< accepted_indices.size();j++){
	  auto part = parts[accepted_indices[j]];
	  if (part->getPid()==321)
	    {
	      if (debug) cout << "accepted K+" << endl;
	      for(int k = 0; k< accepted_indices.size();k++){
		auto part2 = parts[accepted_indices[k]];
		if(part2->getPid()==-211){
		  if (debug) cout << "accepted pi-" << endl;
		  TLorentzVector Kp;
		  TLorentzVector pim;
		  SetLorentzVector(Kp,part, Kmass);
		  SetLorentzVector(pim,part2,pimass);
		  double invmass = (Kp+pim).M();
		  if(debug) cout << "pair mass "  << invmass <<  endl;
		  invmass_Kp_pim->Fill(invmass);
		  double missmass = (target+beam-el-Kp-pim).M();
		  missingmass_Kp_pim->Fill(missmass);
		  Kp_pim_2d->Fill(invmass, missmass);
		  Kp_pim_2d_zoom->Fill(invmass, missmass);
		}
	      }
	    }
	  
	}
	//cout << target.M() << "xxxxx"<< target.E()<< endl;
	//cout << beam.Pz() << "yyyy" << beam.E() << endl;
	//cout << el.Pz()  << "yyyy"  << el.E() << endl;
	//cout <<  "zzzz"  << el.Theta()/3.14159*180 << endl;
	for(int j = 0; j< accepted_indices.size();j++){
          auto part = parts[accepted_indices[j]];
          if (part->getPid()==-321)
            {
              if (debug) cout << "accepted K-" << endl;
              for(int k = 0; k< accepted_indices.size();k++){
                auto part2 = parts[accepted_indices[k]];
                if(part2->getPid()==211){
                  if (debug) cout << "accepted pi+" << endl;
		  for(int l = 0; l< accepted_indices.size();l++){
		    auto part3 = parts[accepted_indices[l]];
		    if(part3->getPid()==2212){
		      if (debug) cout << "accepted p" << endl;
		      TLorentzVector Km;
		      TLorentzVector pip;
		      TLorentzVector prot;
		      SetLorentzVector(Km,part, Kmass);
		      SetLorentzVector(pip,part2,pimass);
                      SetLorentzVector(prot,part3,pmass);
		      double invmass = sqrt((Km+pip+prot)*(Km+pip+prot));
		      if (debug) cout << "pair mass "  << invmass <<  endl;
		      invmass_Km_pip_p->Fill(invmass);
		      double missmass = (target+beam-el-Km-pip-prot).M(); 
		      missingmass_Km_pip_p->Fill(missmass);
		      Km_pip_p_2d->Fill(invmass, missmass);
		      Km_pip_p_2d_zoom->Fill(invmass, missmass);
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
  
  missingmass_Kp_pim->Write();
  invmass_Kp_pim->Write();
  Kp_pim_2d->Write();
  Kp_pim_2d_zoom->Write();
  missingmass_Km_pip_p->Write();
  invmass_Km_pip_p->Write();
  Km_pip_p_2d->Write();
  Km_pip_p_2d_zoom->Write();
  f->Close();
  
  
  
  
  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish - start;
  std::cout << "Elapsed time: " << elapsed.count()<< " events = "<<count<< " s\n";
  
}
