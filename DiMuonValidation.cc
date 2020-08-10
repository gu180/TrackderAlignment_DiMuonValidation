#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
 #include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"


#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.
#include "TH1D.h"
#include "TH3D.h"
#include "TLorentzVector.h"

using reco::TrackCollection;
using namespace std;
using namespace edm;
class DiMuonValidation : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit DiMuonValidation(const edm::ParameterSet& pset)
      {

        cout<<" DiMuonValidation constructor"<<endl;
        TkTag_ = pset.getParameter<string>("TkTag");
        theTrackCollectionToken = consumes<reco::TrackCollection>(TkTag_);
      }
      
      ~DiMuonValidation();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      const double Muon_mass=0.105658; //The invariant mass of muon 105.658MeV
    
      
      std::string TkTag_;
   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      edm::Service<TFileService> fs;
      edm::EDGetTokenT<reco::TrackCollection>  theTrackCollectionToken;

      //==================================================
      const static int eta_bins_number=10;
      const double eta_min=-5;
      const double eta_max=5;
      int Get_Idx_Eta(double eta)
      {
        int idx=-1;
        double r=(eta-eta_min)/(eta_max-eta_min);
        idx=int(r*eta_bins_number);
        cout<<"r="<<r<<"  idx="<<idx<<endl;
        return idx;

      }

      TH3D* th3d_mass_pt_phi[eta_bins_number];//Histogram sotre the invariant mass of mu+mu- pair for different eta, pt, phi bins.
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//



DiMuonValidation::~DiMuonValidation()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
DiMuonValidation::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   cout<<"xixihh"<<endl;
   edm::Handle<reco::TrackCollection> trackCollection;
   iEvent.getByToken(theTrackCollectionToken, trackCollection);
   const reco::TrackCollection tC = *(trackCollection.product());

   TLorentzVector TLVect_mother(0.,0.,0.,0.);
   for (reco::TrackCollection::const_iterator track1=tC.begin(); track1!=tC.end(); track1++)
   {
      TLorentzVector TLVect_track1(track1->px(),track1->py(),track1->pz(),sqrt((track1->p()*track1->p())+(Muon_mass*Muon_mass))); //old 106


      for (reco::TrackCollection::const_iterator track2=track1+1; track2!=tC.end(); track2++)
      {
        if(track1->charge()==track2->charge()) {continue;} // only reconstruct opposite charge pair

        TLorentzVector TLVect_track2(track2->px(),track2->py(),track2->pz(),sqrt((track2->p()*track2->p())+(Muon_mass*Muon_mass)));
        TLVect_mother=TLVect_track1+TLVect_track2;
        int idx_eta= Get_Idx_Eta(TLVect_mother.Eta());
        if(idx_eta<0 or idx_eta>9){continue;}
        th3d_mass_pt_phi[idx_eta]->Fill(TLVect_mother.M(), TLVect_mother.Pt(), TLVect_mother.Phi(), 1);
      }
    }
}


// ------------ method called once each job just before starting event loop  ------------
void
DiMuonValidation::beginJob()
{
  for(int idx_eta=0; idx_eta<eta_bins_number; idx_eta++)
  {
    cout<<"idx_eta="<<idx_eta<<endl;
    th3d_mass_pt_phi[idx_eta]=fs->make<TH3D>(Form("th3d_mass_pt_phi_eta%d",idx_eta),Form("th3d_mass_pt_phi_eta%d",idx_eta),120,60,120,20,0,10,24,-M_PI,M_PI);
    
  }

}

// ------------ method called once each job just after ending the event loop  ------------
void
DiMuonValidation::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DiMuonValidation::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

}

//define this as a plug-in
DEFINE_FWK_MODULE(DiMuonValidation);
