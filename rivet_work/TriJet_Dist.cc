// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// Generic analysis looking at various distributions of final state particles
  ///
  /// @deprecated Replaced by the better-named MC_FSPARTICLES
  class TriJet_Dist : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(TriJet_Dist);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Projections
      // Initialise and register projections

      // The basic final-state projection:
      // all final-state particles within
      // the given eta acceptance


      const FinalState fs;
      declare(fs,"FinalState");
      
      
       // The final-state particles declared above are clustered using FastJet with
      // the anti-kT algorithm and a jet-radius parameter 0.4
      // muons and neutrinos are excluded from the clustering
      FastJets jetfs(fs, FastJets::ANTIKT, 0.4, JetAlg::Muons::ALL, JetAlg::Invisibles::DECAY);
      declare(jetfs, "jets");

      // Histograms      
	  book(mult,"mult", 15, 0., 15.);
	  book(pT_1,"p_T_1",{35, 45, 55, 70, 85, 105, 125, 150, 175, 200, 250, 300, 375, 450});
	  book(pT_2,"p_T_2",{20, 30, 40, 52.5, 65, 80, 95, 120, 150, 200, 250, 300});
	  book(pT_3,"p_T_3",{20, 27.5, 35, 45, 57.5, 70, 85, 105, 150, 200, 250});
	  book(eta_1,"eta_1", {-4.7, -3.8, -3.4, -2.9, -2.5, -2.1, -1.8, -1.6, -1.4, -1.2, -1, -0.8, -0.6, -0.2, 0.2, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2.1, 2.5, 2.9, 3.4, 3.8, 4.7});
	  book(eta_2,"eta_2", {-4.7, -3.8, -3.4, -2.9, -2.5, -2.1, -1.8, -1.6, -1.4, -1.2, -1, -0.8, -0.6, -0.2, 0.2, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2.1, 2.5, 2.9, 3.4, 3.8, 4.7});
	  book(eta_3,"eta_3", {-4.7, -3.8, -3.4, -2.9, -2.5, -2.1, -1.8, -1.6, -1.4, -1.2, -1, -0.8, -0.6, -0.2, 0.2, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2.1, 2.5, 2.9, 3.4, 3.8, 4.7});
	  book(phi_1,"phi_1", 64, 0., 6.4);
	  book(phi_2,"phi_2", 64, 0., 6.4);
	  book(phi_3,"phi_3",  64, 0., 6.4);
	  book(delta_phi,"delta_phi",{0.0, 0.6, 1.3, 1.9, 2.4, 2.8,M_PI});
	  book(phi_dijet,"phi_dijet", 64, 0., 6.4);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
    
      // Retrieve clustered jets, sorted by pT, with a minimum pT cut
          const Jets& jet = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 10.0*GeV);
    double pT1, eta1, phi1;
    double pT2, eta2, phi2;
    double pT3, eta3, phi3;
    double deltaphi,phidijet;

    if(jet.size() > 2) {
		pT1=jet[0].pT()/GeV;
		pT2=jet[1].pT()/GeV;
		pT3=jet[2].pT()/GeV;
		eta1=jet[0].eta();
		eta2=jet[1].eta();
		eta3=jet[2].eta();
		if(pT1>=35&&pT2>=20&&pT3>=20){
				mult->fill( jet.size());	
				pT_1->fill(pT1);
				eta_1->fill(eta1);
				phi1=jet[0].phi();
				if (phi1 < 0) phi1 = 2*M_PI + phi1;
				phi_1->fill(phi1);	
				pT_2->fill(pT2);
				eta_2->fill(eta2);
				phi2=jet[1].phi();
				if (phi2 < 0) phi2 = 2*M_PI + phi2;
				phi_2->fill(phi2);
				pT_3->fill(pT3);
				eta_3->fill(eta3);
				phi3=jet[2].phi();
				if (phi3 < 0) phi3 = 2*M_PI + phi3;
				phi_3->fill(phi3);
				
			if(abs(eta1)>2.0&&abs(eta2)>2.&&abs(eta3)<2.&&eta1*eta2>0.0){	
				phidijet=FourMomentum(jet[0].momentum() + jet[1].momentum()).phi();
				if(phidijet<0.0) phidijet =2*M_PI+phidijet;
				phi_dijet->fill(phidijet);
				deltaphi=abs(phidijet-phi3);
				if(deltaphi>M_PI) deltaphi =2*M_PI-deltaphi;
				delta_phi->fill(deltaphi);

			}
			else if(abs(eta1)>2.0&&abs(eta3)>2.0&&abs(eta2)<2.0&&eta1*eta3>0.0){
				phidijet=FourMomentum(jet[0].momentum() + jet[2].momentum()).phi();
				if(phidijet<0.0) phidijet =2*M_PI+phidijet;
				phi_dijet->fill(phidijet);
				deltaphi=abs(phidijet-phi2);
				if(deltaphi>M_PI) deltaphi =2*M_PI-deltaphi;
				delta_phi->fill(deltaphi);

			}
			else if(abs(eta3)>2.0&&abs(eta2)>2.0&&abs(eta1)<2.0&&eta3*eta2>0.0){
				phidijet=FourMomentum(jet[2].momentum() + jet[1].momentum()).phi();
				if(phidijet<0.0) phidijet =2*M_PI+phidijet;
				phi_dijet->fill(phidijet);
				deltaphi=abs(phidijet-phi1);
				if(deltaphi>M_PI) deltaphi =2*M_PI-deltaphi;
				delta_phi->fill(deltaphi);

			}
		}
	}      
      // fill the jet pT spectra
    }
    /// Finalize
    void finalize() {

    }
    //@}
  private:
    /// @name Histograms

    //@{
    Histo1DPtr mult,pT_1,pT_2,pT_3,eta_1,eta_2,eta_3;
    Histo1DPtr phi_1,  phi_2,phi_3, delta_phi,phi_dijet;

    //@}
  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(TriJet_Dist);

}
