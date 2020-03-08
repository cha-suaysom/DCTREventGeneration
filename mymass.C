/*
This macro shows how to compute jet energy scale.
root -l examples/Example4.C'("delphes_output.root", "plots.root")'
*/

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#else
class ExRootTreeReader;
class ExRootResult;
#endif

class ExRootResult;
class ExRootTreeReader;

#include<iostream>
using namespace std;


// Find the last index using recursion
int find_index(TClonesArray* branchParticle, int first_particle_index, int end_particle_PID) {
	if  (first_particle_index == -1) {
		return -1;
	}
	GenParticle *current_particle = (GenParticle*) branchParticle->At(first_particle_index); 
	
	if (current_particle->PID != end_particle_PID) {
		int d1_index = find_index(branchParticle, current_particle->D1, end_particle_PID);

		if (d1_index == -1) {
			return find_index(branchParticle, current_particle->D2, end_particle_PID);
		}
		else {
			return d1_index;
		}

		
	}
	else {
		return first_particle_index;
	}
}

void AnalyseEvents(ExRootTreeReader *treeReader, const char *outputFile_det, const char *outputFile_part, const char *outputFile_hadron,
	const char *outputFile_part_sub, const char *outputFile_det_sub)

{
  TClonesArray *branchGenJet = treeReader->UseBranch("GenJet");
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack");
  TClonesArray *branchEFlowPhoton = treeReader->UseBranch("EFlowPhoton");
  TClonesArray *branchEFlowNeutralHadron = treeReader->UseBranch("EFlowNeutralHadron");
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchFatJet = treeReader->UseBranch("FatJet");
  TClonesArray *branchTruthFatJet = treeReader->UseBranch("TruthFatJet");

  Long64_t allEntries = treeReader->GetEntries();
  ofstream myfile_det;
  ofstream myfile_part;
  ofstream myfile_all_stable_particles;
  ofstream myfile_part_sub;
  ofstream myfile_det_sub;

  cout << "** Chain contains " << allEntries << " events" << endl;

  Jet *jet, *genjet;
  //GenParticle *particle;
  Jet *particle;
  Jet *particle2;
  TObject *object;

  Track *track;
  Tower *tower;

  TLorentzVector jetMomentum, genJetMomentum, bestGenJetMomentum;

  Float_t deltaR;
  Float_t pt, eta;
  Long64_t entry;

  Int_t i, j;

  myfile_det.open (outputFile_det);
  myfile_part.open (outputFile_part);
  myfile_all_stable_particles.open (outputFile_hadron);
  myfile_part_sub.open (outputFile_part_sub);
  myfile_det_sub.open (outputFile_det_sub);

  // Define the csv header
  myfile_det  << "entry,px,py,pz,E,pid" << endl;
  myfile_part << "entry,index,pT,eta,phi,p,tau_0,tau_1,tau_2,tau_3,tau_4,soft_pT,soft_M,trimmed_pT,trimmed_M,n_trimmed" << endl;
  myfile_all_stable_particles << "entry,pT,eta,phi,pid" << endl;
  myfile_part_sub << "entry,jet_index,px,py,pz,E,pid" << endl;
  myfile_det_sub << "entry,pT/eT,eta,phi,pid" << endl;


  // Loop over all events
  for(entry = 0; entry < allEntries; ++entry)
  {
    if (entry > 1000) break;
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
    if(entry%100 == 0) cout << "Event number: "<< entry <<endl;
     
    // We want to find the truth bbar vdbar and l
    bool find_t = false;
	int t_PID = 6;
	int boson_plus_PID = 24;
	int electron_PID = -11;
	int muon_PID = -13;
	int b_PID = 5;

	bool find_tbar = false;
	int tb_PID = -6;
	int bb_PID = -5;


	 for(int k = 0; k < branchParticle ->GetEntriesFast(); ++k)
      {
      	GenParticle *gen = (GenParticle*) branchParticle->At(k); 

		// t
		if (gen->PID==t_PID && find_t == false) {
			int boson_index = find_index(branchParticle, k, boson_plus_PID);
			int b_index = find_index(branchParticle, k, b_PID);
			int electron_index = find_index(branchParticle, boson_index, electron_PID);
			int muon_index = find_index(branchParticle, boson_index, muon_PID);
			if (boson_index > 0 && b_index > 0 && (electron_index > 0 || muon_index > 0)) {
				cout << "Entry: " << entry << " " << k << " " << boson_index << " " << b_index << " " << electron_index << " " << muon_index << endl;
				find_t = true;
			}
			
		}

		//tbar
		if (gen->PID==tb_PID && find_tbar == false) {
			int boson_index = find_index(branchParticle, k, bb_PID);

			GenParticle *gen_boson = (GenParticle*) branchParticle->At(boson_index);
			
			if (boson_index > 0) {
				cout << "Entry: " << entry << " " << k << " " << boson_index << endl;
				find_tbar = true;
			}
		}

		if (find_t == true && find_tbar == true) {
			break;
		}
		
	  }
  }
}



//------------------------------------------------------------------------------

void mymass(const char *inputFile, const char *outputFile_det, const char *outputFile_part, const char *outputFile_hadron,
	const char *outputFile_part_sub, const char *outputFile_det_sub)
{
  gSystem->Load("libDelphes");

  TChain *chain = new TChain("Delphes");
  chain->Add(inputFile);

  ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
  ExRootResult *result = new ExRootResult();

  AnalyseEvents(treeReader, outputFile_det, outputFile_part, outputFile_hadron,
  		outputFile_part_sub, outputFile_det_sub);

  cout << "** Exiting..." << endl;

  // delete result;
  // delete treeReader;
  // delete chain;
}

//------------------------------------------------------------------------------
