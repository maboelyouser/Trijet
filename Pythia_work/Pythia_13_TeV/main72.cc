  #include "Pythia8/Pythia.h"  
  // Include the Pythia8Rivet header file. 
  #include "Pythia8Plugins/Pythia8Rivet.h" 
  using namespace Pythia8;
 
  int main() { 
 Pythia pythia;
 pythia.readString("HardQCD:all = on");
 pythia.readString("PhaseSpace:pTHatMin = 15.");
 pythia.readString("PhaseSpace:pTHatMax = 1000.");
 pythia.readString("Tune:pp = 14");
 pythia.readString("PDF:pSet = LHAPDF6:NNPDF23_nlo_as_0119_qed");
  // Pick new random number seed for each run, based on clock.
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = 1");
  
  // No event record printout.
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");

  // LHC initialization.
  pythia.readString("Beams:eCM = 13000.");
 
    pythia.init(); 
 
    // Create a Pythia8Rivet object and add (one or several) analyses. 
 
    Pythia8Rivet rivet(pythia, "Pythia_13_TeV.yoda"); 

    rivet.addAnalysis("TriJet_Dist");
    //rivet.addAnalysis("TriJet_Cross");
    //rivet.addAnalysis("TriJet_Norm");
 
    for (int iEvent = 0; iEvent < 10000000; ++iEvent) { 
      if (!pythia.next()) continue; 
 
      // Push event to Rivet. 
      rivet(); 
 
      // Maybe do other non-Rivet analysis. 
    } 
 
    // Tell Rivet to finalise the run. 
    rivet.done(); 
 
  }  



