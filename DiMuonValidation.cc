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
#include "TH2D.h"
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
        pair_mass_min   = pset.getParameter<double>("Pair_mass_min");
        pair_mass_max   = pset.getParameter<double>("Pair_mass_max");
        pair_mass_bins  = pset.getUntrackedParameter<int>("Pair_mass_bins",120);
      }
      
      ~DiMuonValidation();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      const double Muon_mass=0.105658; //The invariant mass of muon 105.658MeV
    
      
      
      //===========Parameter==============================
      std::string TkTag_;
      double pair_mass_min=60;
      double pair_mass_max=120;
      int    pair_mass_bins=120;

      double pair_pt_min;
      double pair_pt_max;
      int    pair_pt_bins;

      double mu_p_pt_min;
      double mu_p_pt_max;
      int    mu_p_pt_bins;

      double mu_n_pt_min;
      double mu_n_pt_max;
      int    mu_n_pt_bins;

      


      //==================================================
   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      edm::Service<TFileService> fs;
      edm::EDGetTokenT<reco::TrackCollection>  theTrackCollectionToken;


     
      
      const static int variables_number=8;
      TH2D* th2d_mass_variables[variables_number];
      TString tstring_variables_name[variables_number]={"CosThetaCS","DeltaEta","EtaMinus","EtaPlus","PhiCS","PhiMinus","PhiPlus","Pt"};

      const int variables_bins_number[variables_number]={20, 20, 20  , 20  , 20  , 20   , 20   , 100};
      const double variables_min[variables_number]     ={-1, -5, -5  , -5  , -2  , -5   , -5   , 0  };
      const double variables_max[variables_number]     ={+1, +5, +5  , +5  , +2  , +5   , +5   , 400};

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
  // double variables_val[variables_number];
   for (reco::TrackCollection::const_iterator track1=tC.begin(); track1!=tC.end(); track1++)
   {
      TLorentzVector TLVect_track1(track1->px(),track1->py(),track1->pz(),sqrt((track1->p()*track1->p())+(Muon_mass*Muon_mass))); //old 106

      for (reco::TrackCollection::const_iterator track2=track1+1; track2!=tC.end(); track2++)
      {
        if(track1->charge()==track2->charge()) {continue;} // only reconstruct opposite charge pair

        TLorentzVector TLVect_track2(track2->px(),track2->py(),track2->pz(),sqrt((track2->p()*track2->p())+(Muon_mass*Muon_mass)));
        TLVect_mother=TLVect_track1+TLVect_track2;
        double mother_mass=TLVect_mother.M();
        double mother_pt=TLVect_mother.Pt();

        int charge1   = track1->charge();
        double etaMu1 = track1->eta();
        double phiMu1 = track1->phi();
        double ptMu1  = track1->pt();

        int charge2   = track2->charge();
        double etaMu2 = track2->eta();
        double phiMu2 = track2->phi();
        double ptMu2  = track2->pt();


        if(charge1<0){ // use Mu+ for charge1, Mu- for charge2
        swap(charge1,charge2);
        swap(etaMu1,etaMu2);
        swap(phiMu1,phiMu2);
        swap(ptMu1,ptMu2);
        }
        double delta_eta=etaMu1-etaMu2;
        
        double muplus = 1.0 / sqrt(2.0) * (TLVect_track1.E() + TLVect_track1.Z());
        double muminus = 1.0 / sqrt(2.0) * (TLVect_track1.E() - TLVect_track1.Z());
        double mubarplus = 1.0 / sqrt(2.0) * (TLVect_track2.E() + TLVect_track2.Z());
        double mubarminus = 1.0 / sqrt(2.0) * (TLVect_track2.E() - TLVect_track2.Z());
        //double costheta = 2.0 / Q.Mag() / sqrt(pow(Q.Mag(), 2) + pow(Q.Pt(), 2)) * (muplus * mubarminus - muminus * mubarplus);
        double costhetaCS = 2.0 / TLVect_mother.Mag() / sqrt(pow(TLVect_mother.Mag(), 2) + pow(TLVect_mother.Pt(), 2)) * (muplus * mubarminus - muminus * mubarplus);
        
        TLorentzVector Pbeam(0., 0., 3500., 3500.);
        TVector3 R = Pbeam.Vect().Cross(TLVect_mother.Vect());
        TVector3 Runit = R.Unit();
        TVector3 Qt = TLVect_mother.Vect();
        Qt.SetZ(0);
        TVector3 Qtunit = Qt.Unit();
        TLorentzVector D(TLVect_track1 - TLVect_track2);
        TVector3 Dt = D.Vect();
        Dt.SetZ(0);
        double tanphi = sqrt(pow(TLVect_mother.Mag(), 2) + pow(TLVect_mother.Pt(), 2)) / TLVect_mother.Mag() * Dt.Dot(Runit) / Dt.Dot(Qtunit);
        double  phiCS = atan(tanphi);

        if(mother_mass> pair_mass_min and mother_mass< pair_mass_max)
        {
          th2d_mass_variables[0]->Fill(mother_mass, costhetaCS,1);
          th2d_mass_variables[1]->Fill(mother_mass, delta_eta ,1);
          th2d_mass_variables[2]->Fill(mother_mass, etaMu2,1);
          th2d_mass_variables[3]->Fill(mother_mass, etaMu1,1);
          th2d_mass_variables[4]->Fill(mother_mass, phiCS,1);
          th2d_mass_variables[5]->Fill(mother_mass, phiMu2,1);
          th2d_mass_variables[6]->Fill(mother_mass, phiMu1,1);
          th2d_mass_variables[7]->Fill(mother_mass, mother_pt,1);
        }
       


      }
    }
}


// ------------ method called once each job just before starting event loop  ------------
void
DiMuonValidation::beginJob()
{
  for(int i=0; i<variables_number; i++)
  {
    TString th2d_name=Form("th2d_mass_%s",tstring_variables_name[i].Data());
    th2d_mass_variables[i]=fs->make<TH2D>(th2d_name, th2d_name,pair_mass_bins,pair_mass_min,pair_mass_max, variables_bins_number[i], variables_min[i], variables_max[i]);
    //cout<<"variables_bins_number[i], variables_min[i], variables_max[i]: "<< variables_bins_number[i]<<" "<< variables_min[i]<<" "<< variables_max[i]<<endl;
    //th2d_mass_variables[i]=fs->make<TH2D>(th2d_name, th2d_name,120,60,120, variables_bins_number[i], variables_min[i], variables_max[i]);
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
