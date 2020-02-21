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

//------------------------------------------------------------------------------

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
  ofstream myfile_hadron;
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
  myfile_hadron.open (outputFile_hadron);
  myfile_part_sub.open (outputFile_part_sub);
  myfile_det_sub.open (outputFile_det_sub);

  // Define the csv header
  myfile_part << "entry,index,pT,eta,phi,p,tau_0,tau_1,tau_2,tau_3,tau_4,soft_pT,soft_M,trimmed_pT,trimmed_M,n_trimmed" << endl;
  myfile_part_sub << "entry, pT,eta,phi,pid" << endl;
  myfile_det  << "entry,index,pT,eta,phi,p,tau_0,tau_1,tau_2,tau_3,tau_4,soft_pT,soft_M,trimmed_pT,trimmed_M,n_trimmed" << endl;
  myfile_det_sub << "entry,pT/eT,eta,phi,pid" << endl;
  myfile_hadron << "entry,status,pT,eta,phi,m,m1,m2,pid" << endl;

  // Loop over all events
  for(entry = 0; entry < allEntries; ++entry)
  {
    if (entry > 10000) break;
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    if(entry%500 == 0) cout << "Event number: "<< entry <<endl;
     

    //Loop over all reconstructed jets in event
    for(i = 0; i < branchFatJet->GetEntriesFast(); ++i)
      {

	//if (i > 0) continue; //let's just take the leading jet.
	
	jet = (Jet*) branchFatJet->At(i);
	jetMomentum = jet->P4();

	deltaR = 999;
	int whichjet = 0;
	



	// Loop over all hard partons in event
	for(j = 0; j < branchTruthFatJet->GetEntriesFast(); ++j)
	  {
	    particle = (Jet*) branchTruthFatJet->At(j);
	    
	    genJetMomentum = particle->P4();
	    
	    // take the closest parton candidate
	    if(genJetMomentum.DeltaR(jetMomentum) < deltaR)
	      {
		deltaR = genJetMomentum.DeltaR(jetMomentum);
		bestGenJetMomentum = genJetMomentum;
		particle2 = (Jet*) branchTruthFatJet->At(j);
		whichjet = j;
	      }
	  }
	


	if(deltaR < 0.3)
	  {
	    pt  = jetMomentum.Pt();
	    eta = TMath::Abs(jetMomentum.Eta());
	    
	    //For the part part

	    // Define csv header
	    myfile_part << entry << ",";

	    myfile_part << whichjet << "," 
	    			<< particle2->PT << "," 
	    			<< particle2->Eta << "," 
	    			<< particle2->Phi << "," 
	    			<< jetMomentum.M() << ",";

	    myfile_part << particle2->Tau[0] << "," 
	    			<< particle2->Tau[1] << "," 
	    			<< particle2->Tau[2] << "," 
	    			<< particle2->Tau[3] << "," 
	    			<< particle2->Tau[4] << ",";

	    myfile_part << particle2->SoftDroppedP4.Pt() << "," 
	    			<< particle2->SoftDroppedP4.M() << "," 
	    			<< particle2->TrimmedP4.Pt() << "," 
	    			<< particle2->TrimmedP4.M() << "," 
	    			<< particle2->NSubJetsTrimmed << endl;
	

	    

	    for(int j2 = 0; j2 < particle2->Constituents.GetEntriesFast(); ++j2)
	      {
			object = particle2->Constituents.At(j2);
			if(object == 0){

	                }
			if(object->IsA() == GenParticle::Class())
			  {
			    subparticle = (GenParticle*) object;
	   			myfile_part_sub << entry << ",";
			    myfile_part_sub << subparticle->PT << "," 
			    			<< subparticle->Eta << "," 
			    			<< subparticle->Phi << "," 
			    			<< subparticle->PID << endl;
			  }
	      }

	    

	    //Now for the det part.
	    myfile_det << entry << ",";

	    myfile_det << i << "," 
	    		   << jet->PT << "," 
	    		   << jet->Eta << "," 
	    		   << jet->Phi << "," 
	    		   << bestGenJetMomentum.M() << ",";

	    myfile_det << jet->Tau[0] << "," 
	    		   << jet->Tau[1] << "," 
	    		   << jet->Tau[2] << "," 
	    		   << jet->Tau[3] << "," 
	    		   << jet->Tau[4] << ",";
	    
	    myfile_det << jet->SoftDroppedP4.Pt() << "," 
	    		   << jet->SoftDroppedP4.M() << "," 
	    		   << jet->TrimmedP4.Pt() << "," 
	    		   << jet->TrimmedP4.M() << "," 
	    		   << jet->NSubJetsTrimmed << endl;

	    for(int j2 = 0; j2 < jet->Constituents.GetEntriesFast(); ++j2)
	      {
	      	
			object = jet->Constituents.At(j2);
			if(object == 0)
				{  
				}
		
		

		else if(object->IsA() == Track::Class()){
		  track = (Track*) object;
		  if (abs(track->PID)==11){
		  	myfile_det_sub << entry << ",";
		    myfile_det_sub << track->PT << "," 
		    		   << track->Eta << "," 
		    		   << track->Phi << "," 
		    		   << "11" << endl;  
		  }
		  else if (abs(track->PID)==13){
		  	myfile_det_sub << entry << ",";
		    myfile_det_sub << track->PT << "," 
		    		   << track->Eta << "," 
		    		   << track->Phi << "," 
		    		   << "13" << endl; 
		  }
		  else{
		  	myfile_det_sub << entry << ",";
		    myfile_det_sub << track->PT << "," 
		    		   << track->Eta << "," 
		    		   << track->Phi << "," 
		    		   << "211" << endl; 
		  }
		}
		else if(object->IsA() == Tower::Class()){

		  tower = (Tower*) object;
		  if (tower->Eem/(tower->Eem+tower->Ehad) == 1){
		  	myfile_det_sub << entry << ",";
		    myfile_det_sub << tower->ET << "," 
		    		   << tower->Eta << "," 
		    		   << tower->Phi << "," 
		    		   << "22" << endl; 
		  }
		  else{
		  	myfile_det_sub << entry << ",";
		    myfile_det_sub << tower->ET << "," 
		               << tower->Eta << "," 
		               << tower->Phi << ","
		               << "2112" << endl; 
		  }
		}
	      }
	    
	  }
      }

      
      // Output every single stable hadron
      for(int k = 0; k < branchParticle ->GetEntriesFast(); ++k)
	      {
	      	GenParticle *gen = (GenParticle*) branchParticle->At(k); 
	      	//Final state gen_particle
	      	if ( gen->Status!=1 ) continue;    
     		//skip electron type neutrinos
    		if (abs(gen->PID)==12 ) continue; 
			//skip muon type neutrinos
			if (abs(gen->PID)==14 ) continue; 
			//skip tau type neutrinos
			if (abs(gen->PID)==16 ) continue;
			//save 4-vector and PdgID
			
			myfile_hadron << entry << ",";

			myfile_hadron << gen->Status << "," 
						  << gen->PT << "," 
						  << gen->Eta << "," 
						  << gen->Phi << "," 
						  << gen->Mass << "," 
						  << gen->M1 <<"," 
						  << gen->M2 << "," 
						  << gen->PID << endl;
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

  delete result;
  delete treeReader;
  delete chain;
}

//------------------------------------------------------------------------------
