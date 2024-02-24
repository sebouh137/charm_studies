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

int torus;


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
  if(abs(h_chi2pid) >3*C)
    return false;
  
  //tighter cut for pions
  if(abs(h_pid) == 211 && h_p>2.44 && h_chi2pid >(0.00869+14.98587*exp(-h_p/1.18236)+1.81751*exp(-h_p/4.86394))*C)
    return false;
  
  auto dvz = electron->par()->getVz()-h->par()->getVz();
  
  if(dvz < cut_dvzmin || dvz > cut_dvzmax)
    return false;
}



void TuplesForMassResolution(){
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

  leaf(invmass);
  leaf(missmass);

  //5=K+K-p
  //6=pi+pi-p
  leaf(topo);

  //count the number of each of these types of particles
  leaf(n_pim);
  leaf(n_pip);
  leaf(n_p);
  leaf(n_Kp);
  leaf(n_Km);

//macro to create one column for each particle.  
#define leaf3(name) leaf(h1_##name); leaf(h2_##name); leaf(prot_##name);

  leaf3(p);
  leaf3(px);
  leaf3(py);
  leaf3(pz);
  leaf3(E);
  leaf3(ph);
  leaf3(th);
#define fillVec(v,pref) pref##_p=v.P(); pref##_px=v.Px(); pref##_py=v.Py(); pref##_pz=v.Pz(); pref##_E=v.E(); pref##_th=v.Theta(); pref##_ph=v.Phi();
#define blankVec(pref) pref##_p=0; pref##_px=0; pref##_py=0; pref##_pz=0; pref##_E=0; pref##_th=0; pref##_ph=0;

  //time
  //leaf(e_t);
  leaf(prot_dt);
  leaf(h1_dt);
  leaf(h2_dt);
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
  
  TLorentzVector had;
  
  // leptonic variables
  TH1* h_nu = new TH1D("nu", "#nu;#nu [GeV];events", 100, 0, 11);
  TH1* h_Q2 = new TH1D("Q2", "Q^{2};Q^{2} [GeV^2];events", 100,0, 11);
  TH1* h_W = new TH1D("W", "W;W [GeV];events", 100, 0, 11);
  TH2* h_Q2_nu = new TH2D("Q2_nu", "Q^{2} vs #nu;Q^{2} [GeV^2];#nu [GeV]",100, 0, 11, 100, 0, 11);

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
      
      
      // get particles by type
      
      auto electrons=c12.getByID(11);
      
      TLorentzVector el(0,0,0,db->GetParticle(11)->Mass());
      

      
      auto parts=c12.getDetParticles();


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
          
          //if(debug) cout << "pip" << endl;
          double mass = db->GetParticle(h_pid)->Mass();
          SetLorentzVector(had,h, mass);
          double h_p = had.P();
          double c = 29.9792458; //cm/ns
          auto dtime_corr =dtime-electrons[i]->getPath()/c+h->getPath()/(had.Beta()*c);
          if(abs(dtime_corr) > cut_dtime_corr)
            continue;
          dtime_map[j]=dtime_corr;
          
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
	  
	  if (h_pid==211) n_pip+=1;
	  if (h_pid==-211) n_pim+=1;
	  if (h_pid==321) n_Kp+=1;
          if (h_pid==-321) n_Km+=1;
	  if (h_pid==2212) n_p+=1;

        }
        double Kmass=db->GetParticle(321)->Mass();
        double pimass=db->GetParticle(211)->Mass();
        double pmass=db->GetParticle(2212)->Mass();

	
	// p K+ K-
        for(int j = 0; j< accepted_indices.size();j++){
          auto part = parts[accepted_indices[j]];
          if (part->getPid()==321){
            if (debug) cout << "accepted K+" << endl;
            for(int k = 0; k< accepted_indices.size();k++){
              auto part2 = parts[accepted_indices[k]];
              if(part2->getPid()==-321){
                if (debug) cout << "accepted K-" << endl;
                for(int l = 0; l< accepted_indices.size();l++){
                  auto part3 = parts[accepted_indices[l]];
                  if(part3->getPid()==2212){
                    if (debug) cout << "accepted p" << endl;
                    TLorentzVector Km;
                    TLorentzVector Kp;
                    TLorentzVector prot;
                    SetLorentzVector(Kp,part, Kmass);
                    SetLorentzVector(Km,part2,Kmass);
                    SetLorentzVector(prot,part3,pmass);
                    invmass = (Km+Kp+prot).M();
                    missmass = (target+beam-el-Km-Kp-prot).M();
		    fillVec(Kp,h1);
		    fillVec(Km,h2);
		    fillVec(prot,prot);
		    prot_dt=dtime_map[accepted_indices[l]];
		    h1_dt=dtime_map[accepted_indices[j]];
		    h2_dt=dtime_map[accepted_indices[k]];
		    topo=5;
		    tree->Fill();
                  }
                }
              }
            }
          }  
        }
	// p pi+ pi-
        for(int j = 0; j< accepted_indices.size();j++){
          auto part = parts[accepted_indices[j]];
          if (part->getPid()==211){
	    if (debug) cout << "accepted K-" << endl;
	    for(int k = 0; k< accepted_indices.size();k++){
	      auto part2 = parts[accepted_indices[k]];
	      if(part2->getPid()==-211){
		if (debug) cout << "accepted pi+" << endl;
		for(int l = 0; l< accepted_indices.size();l++){
		  auto part3 = parts[accepted_indices[l]];
		  if(part3->getPid()==2212){
		    if (debug) cout << "accepted p" << endl;
		    TLorentzVector pim;
		    TLorentzVector pip;
		    TLorentzVector prot;
		    SetLorentzVector(pim,part, pimass);
		    SetLorentzVector(pip,part2,pimass);
		    SetLorentzVector(prot,part3,pmass);
		    invmass = (pim+pip+prot).M();
		    missmass = (target+beam-el-pim-pip-prot).M();
		    fillVec(pip,h1);
		    fillVec(pim,h2);
		    fillVec(prot,prot);
		    prot_dt=dtime_map[accepted_indices[l]];
		    h1_dt=dtime_map[accepted_indices[j]];
		    h2_dt=dtime_map[accepted_indices[k]];
		    topo=6;
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
  f->Close();  
  
  
  
  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish - start;
  std::cout << "Elapsed time: " << elapsed.count()<< " events = "<<count<< " s\n";
  
}
