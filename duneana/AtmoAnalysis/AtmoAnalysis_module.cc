////////////////////////////////////////////////////////////////////////
// Class:       atmoAnalysis
// Plugin Type: analyzer (art v3_05_01)
// File:        AtmoAnalysis_module.cc
//
// Generated at Fri May  3 15:01:35 2024 by Pierre Granger using cetskelgen
// from cetlib version v3_10_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art_root_io/TFileService.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/ArtDataHelper/TrackUtils.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "larsim/Utils/TruthMatchUtils.h"
#include "larsim/MCCheater/BackTrackerService.h"

#include "dunereco/AnaUtils/DUNEAnaEventUtils.h"
#include "dunereco/AnaUtils/DUNEAnaPFParticleUtils.h"
#include "dunereco/AnaUtils/DUNEAnaTrackUtils.h"
#include "dunereco/AnaUtils/DUNEAnaShowerUtils.h"
#include "dunereco/FDSensOpt/FDSensOptData/AngularRecoOutput.h"
#include "dunereco/FDSensOpt/FDSensOptData/EnergyRecoOutput.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include <TTree.h>
#include <vector>
#include <string>
#include <TH1D.h>
#include <TLorentzVector.h>


namespace test {
  class atmoAnalysis;
}


class test::atmoAnalysis : public art::EDAnalyzer {
public:
  explicit atmoAnalysis(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  atmoAnalysis(atmoAnalysis const&) = delete;
  atmoAnalysis(atmoAnalysis&&) = delete;
  atmoAnalysis& operator=(atmoAnalysis const&) = delete;
  atmoAnalysis& operator=(atmoAnalysis&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.

  TTree *fTree;

  unsigned int fEventID;
  unsigned int fRunID;
  unsigned int fSubRunID;

  unsigned int fNMCParticles;
  unsigned int fNPFParticles;

  double fTrueVtx_x;
  double fTrueVtx_y;
  double fTrueVtx_z;
  double fTrueNuP_x;
  double fTrueNuP_y;
  double fTrueNuP_z;
  double fTrueNuPdg;
  double fTrueNuE;
  bool fIsCC;
  int fInterMode;

  double fRecoVtx_x;
  double fRecoVtx_y;
  double fRecoVtx_z;
  double fDirectionRecNuE_x;
  double fDirectionRecNuE_y;
  double fDirectionRecNuE_z;
  double fDirectionRecNuMu_x;
  double fDirectionRecNuMu_y;
  double fDirectionRecNuMu_z;
  double fDirectionRecNuEPfps_x;
  double fDirectionRecNuEPfps_y;
  double fDirectionRecNuEPfps_z;
  double fDirectionRecNuMuPfps_x;
  double fDirectionRecNuMuPfps_y;
  double fDirectionRecNuMuPfps_z;
  double fDirectionRecHits_x;
  double fDirectionRecHits_y;
  double fDirectionRecHits_z;
  double fCVNScoreNuMu;
  double fCVNScoreNuE;
  double fCVNScoreNC;
  double fCVNScoreProton0;
  double fCVNScoreProton1;
  double fCVNScoreProton2;
  double fCVNScoreProton3;
  double fCVNScorePion0;
  double fCVNScorePion1;

  double fErecNuMu;
  double fErecNuE;
  double fErecNuMuRange;
  double fErecNuMuMCS;
  double fErecNC;
  double fVisibleEnergy;

  std::string fMCTruthLabel;
  std::string fPandoraNuVertexModuleLabel;
  std::string fDirectionRecoLabelNuMu;
  std::string fDirectionRecoLabelNuE;
  std::string fDirectionRecoLabelNuMuPfps;
  std::string fDirectionRecoLabelNuEPfps;
  std::string fDirectionRecoLabelHits;
  std::string fEnergyRecoNuELabel;
  std::string fEnergyRecoNuMuLabel;
  std::string fEnergyRecoNuMuMCSLabel;
  std::string fEnergyRecoNuMuRangeLabel;
  std::string fEnergyRecoNCLabel;
  std::string fCVNLabel;
  std::string fEdepLabel;

  const geo::Geometry* fGeom;
};


test::atmoAnalysis::atmoAnalysis(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}
{
  fMCTruthLabel = p.get<std::string>("MCTruthLabel");

  fPandoraNuVertexModuleLabel = p.get<std::string>("PandoraNuVertexModuleLabel");
  fDirectionRecoLabelNuMu = p.get<std::string>("DirectionRecoLabelNuMu");
  fDirectionRecoLabelNuE = p.get<std::string>("DirectionRecoLabelNuE");
  fDirectionRecoLabelNuMuPfps = p.get<std::string>("DirectionRecoLabelNuMuPfps");
  fDirectionRecoLabelNuEPfps = p.get<std::string>("DirectionRecoLabelNuEPfps");
  fDirectionRecoLabelHits = p.get<std::string>("DirectionRecoLabelHits");

  fEnergyRecoNuELabel = p.get<std::string>("EnergyRecoNuELabel");
  fEnergyRecoNuMuLabel = p.get<std::string>("EnergyRecoNuMuLabel");
  fEnergyRecoNuMuMCSLabel = p.get<std::string>("EnergyRecoNuMuMCSLabel");
  fEnergyRecoNuMuRangeLabel = p.get<std::string>("EnergyRecoNuMuRangeLabel");
  fEnergyRecoNCLabel = p.get<std::string>("EnergyRecoNCLabel");

  fCVNLabel = p.get<std::string>("CVNLabel");
  fEdepLabel = p.get<std::string>("EdepLabel");
  fGeom    = &*art::ServiceHandle<geo::Geometry>();
} 

void test::atmoAnalysis::analyze(art::Event const& evt)
{
  const art::ServiceHandle<cheat::BackTrackerService> btServ;
  fEventID = evt.id().event();
  fRunID = evt.id().run();
  fSubRunID = evt.id().subRun();
  // auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);
  std::cout << "=============== EVENT ID " << fEventID << " == RUN ID " << fRunID << " == SUBRUN ID " << fSubRunID << " ================" << std::endl;

  //True neutrino infos
  art::Handle<std::vector<simb::MCTruth>> mct = evt.getHandle< std::vector<simb::MCTruth> >(fMCTruthLabel);

  const simb::MCNeutrino &neutrino = mct->at(0).GetNeutrino();
  fTrueNuE = neutrino.Nu().E();
  fTrueNuPdg = neutrino.Nu().PdgCode();
  fIsCC = !(neutrino.CCNC()); // ccnc is 0=CC 1=NC
  fInterMode = neutrino.Mode();
  fTrueVtx_x = neutrino.Nu().Vx();
  fTrueVtx_y = neutrino.Nu().Vy();
  fTrueVtx_z = neutrino.Nu().Vz();
  fTrueNuP_x = neutrino.Nu().Momentum().X();
  fTrueNuP_y = neutrino.Nu().Momentum().Y();
  fTrueNuP_z = neutrino.Nu().Momentum().Z();


  //Reco infos

  lar_pandora::PFParticleVector particleVector;
  lar_pandora::LArPandoraHelper::CollectPFParticles(evt, fPandoraNuVertexModuleLabel, particleVector);
  lar_pandora::VertexVector vertexVector;
  lar_pandora::PFParticlesToVertices particlesToVertices;
  lar_pandora::LArPandoraHelper::CollectVertices(evt, fPandoraNuVertexModuleLabel, vertexVector, particlesToVertices);

  fRecoVtx_x = -999;
  fRecoVtx_y = -999;
  fRecoVtx_z = -999;

  for (unsigned int n = 0; n < particleVector.size(); ++n) {
      const art::Ptr<recob::PFParticle> particle = particleVector.at(n);
      if(particle->IsPrimary()){
        //Retrieving the reco vertex
        lar_pandora::PFParticlesToVertices::const_iterator vIter = particlesToVertices.find(particle);
        if (particlesToVertices.end() != vIter) {
          const lar_pandora::VertexVector &vertexVector = vIter->second;
          if (vertexVector.size() == 1) {
            const art::Ptr<recob::Vertex> vertex = *(vertexVector.begin());
            double xyz[3] = {0.0, 0.0, 0.0} ;
            vertex->XYZ(xyz);
            fRecoVtx_x = xyz[0];
            fRecoVtx_y = xyz[1];
            fRecoVtx_z = xyz[2];
          }
        }

        //Fill direction infos
        art::Handle<dune::AngularRecoOutput> dirReco = evt.getHandle<dune::AngularRecoOutput>(fDirectionRecoLabelNuMu);
        if(!dirReco.failedToGet()){
          fDirectionRecNuMu_x = dirReco->fRecoDirection.X();
          fDirectionRecNuMu_y = dirReco->fRecoDirection.Y();
          fDirectionRecNuMu_z = dirReco->fRecoDirection.Z();
        }
        else{
          mf::LogWarning("CAFMaker") << "No AngularRecoOutput found with label '" << fDirectionRecoLabelNuMu << "'";
        }

        dirReco = evt.getHandle<dune::AngularRecoOutput>(fDirectionRecoLabelNuE);
        if(!dirReco.failedToGet()){
          fDirectionRecNuE_x = dirReco->fRecoDirection.X();
          fDirectionRecNuE_y = dirReco->fRecoDirection.Y();
          fDirectionRecNuE_z = dirReco->fRecoDirection.Z();
        }
        else{
          mf::LogWarning("CAFMaker") << "No AngularRecoOutput found with label '" << fDirectionRecoLabelNuE << "'";
        }

        dirReco = evt.getHandle<dune::AngularRecoOutput>(fDirectionRecoLabelNuMuPfps);
        if(!dirReco.failedToGet()){
          fDirectionRecNuMuPfps_x = dirReco->fRecoDirection.X();
          fDirectionRecNuMuPfps_y = dirReco->fRecoDirection.Y();
          fDirectionRecNuMuPfps_z = dirReco->fRecoDirection.Z();
        }
        else{
          mf::LogWarning("CAFMaker") << "No AngularRecoOutput found with label '" << fDirectionRecoLabelNuMuPfps << "'";
        }

        dirReco = evt.getHandle<dune::AngularRecoOutput>(fDirectionRecoLabelNuEPfps);
        if(!dirReco.failedToGet()){
          fDirectionRecNuEPfps_x = dirReco->fRecoDirection.X();
          fDirectionRecNuEPfps_y = dirReco->fRecoDirection.Y();
          fDirectionRecNuEPfps_z = dirReco->fRecoDirection.Z();
        }
        else{
          mf::LogWarning("CAFMaker") << "No AngularRecoOutput found with label '" << fDirectionRecoLabelNuEPfps << "'";
        }

        dirReco = evt.getHandle<dune::AngularRecoOutput>(fDirectionRecoLabelHits);
        if(!dirReco.failedToGet()){
          fDirectionRecHits_x = dirReco->fRecoDirection.X();
          fDirectionRecHits_y = dirReco->fRecoDirection.Y();
          fDirectionRecHits_z = dirReco->fRecoDirection.Z();
        }
        else{
          mf::LogWarning("CAFMaker") << "No AngularRecoOutput found with label '" << fDirectionRecoLabelNuEPfps << "'";
        }

        fCVNScoreNuE = -999;
        fCVNScoreNuMu = -999;
        fCVNScoreNC = -999;
        //Filling CVN
        art::Handle<std::vector<cvn::Result>> cvnin = evt.getHandle<std::vector<cvn::Result>>(fCVNLabel);
        if( !cvnin.failedToGet() && !cvnin->empty()) {
          const std::vector<std::vector<float>> &scores = (*cvnin)[0].fOutput;

          fCVNScoreNC = scores[0][0];
          fCVNScoreNuE = scores[0][1];
          fCVNScoreNuMu = scores[0][2];
          fCVNScoreProton0 = scores[1][0];
          fCVNScoreProton1 = scores[1][1];
          fCVNScoreProton2 = scores[1][2];
          fCVNScoreProton3 = scores[1][3];
          fCVNScorePion0 = scores[2][0];
          fCVNScorePion1 = scores[2][1];
        }
         else{
          mf::LogWarning("CAFMaker") << "No CVNResult found with label '" << fCVNLabel << "'";
        }

        //Neutrino energy hypothese
        art::Handle<dune::EnergyRecoOutput> ereco = evt.getHandle<dune::EnergyRecoOutput>(fEnergyRecoNuELabel);
        if(ereco.failedToGet()){
          mf::LogWarning("CAFMaker") << fEnergyRecoNuELabel << " does not correspond to a valid EnergyRecoOutput product";
        }
        else{
          fErecNuE = ereco->fNuLorentzVector.E();
        }

        ereco = evt.getHandle<dune::EnergyRecoOutput>(fEnergyRecoNuMuLabel);
        if(ereco.failedToGet()){
          mf::LogWarning("CAFMaker") << fEnergyRecoNuMuLabel << " does not correspond to a valid EnergyRecoOutput product";
        }
        else{
          fErecNuMu = ereco->fNuLorentzVector.E();
        }

        ereco = evt.getHandle<dune::EnergyRecoOutput>(fEnergyRecoNuMuMCSLabel);
        if(ereco.failedToGet()){
          mf::LogWarning("CAFMaker") << fEnergyRecoNuMuMCSLabel << " does not correspond to a valid EnergyRecoOutput product";
        }
        else{
          fErecNuMuMCS = ereco->fNuLorentzVector.E();
        }

        ereco = evt.getHandle<dune::EnergyRecoOutput>(fEnergyRecoNuMuRangeLabel);
        if(ereco.failedToGet()){
          mf::LogWarning("CAFMaker") << fEnergyRecoNuMuRangeLabel << " does not correspond to a valid EnergyRecoOutput product";
        }
        else{
          fErecNuMuRange = ereco->fNuLorentzVector.E();
        }

        ereco = evt.getHandle<dune::EnergyRecoOutput>(fEnergyRecoNCLabel);
        if(ereco.failedToGet()){
          mf::LogWarning("CAFMaker") << fEnergyRecoNCLabel << " does not correspond to a valid EnergyRecoOutput product";
        }
        else{
          fErecNC = ereco->fNuLorentzVector.E();
        }

        //Compute visible energy
        fVisibleEnergy = 0;
        auto theseProds = evt.getHandle<std::vector<sim::SimEnergyDeposit>>(fEdepLabel);
        if(!theseProds.isValid()){
          mf::LogWarning("CAFMaker") << "No SimEnergyDeposit found with label '" << fEdepLabel << "'";
          
        }

        std::vector<art::Ptr<sim::SimEnergyDeposit>> edeps;
        art::fill_ptr_vector(edeps,theseProds);
        std::cout << "Number of edeps: " << edeps.size() << std::endl;
        for(const art::Ptr<sim::SimEnergyDeposit> &edep : edeps){
          fVisibleEnergy += edep->Energy();
        }
      }

    }

  fTree->Fill();
}

void test::atmoAnalysis::beginJob()
{

  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;
  fTree = tfs->make<TTree>("atmoOutput","Atmo Output Tree");

  //Event branches
  fTree->Branch("eventID",&fEventID,"eventID/i");
  fTree->Branch("runID",&fRunID,"runID/i");
  fTree->Branch("subrunID",&fSubRunID,"subrunID/i");

  fTree->Branch("TrueVtx_x", &fTrueVtx_x);
  fTree->Branch("TrueVtx_y", &fTrueVtx_y);
  fTree->Branch("TrueVtx_z", &fTrueVtx_z);
  fTree->Branch("TrueNuP_x", &fTrueNuP_x);
  fTree->Branch("TrueNuP_y", &fTrueNuP_y);
  fTree->Branch("TrueNuP_z", &fTrueNuP_z);
  fTree->Branch("TrueNuPdg", &fTrueNuPdg);
  fTree->Branch("TrueNuE", &fTrueNuE);
  fTree->Branch("RecoVtx_x", &fRecoVtx_x);
  fTree->Branch("RecoVtx_y", &fRecoVtx_y);
  fTree->Branch("RecoVtx_z", &fRecoVtx_z);
  fTree->Branch("CVNScoreNuMu", &fCVNScoreNuMu);
  fTree->Branch("CVNScoreNuE", &fCVNScoreNuE);
  fTree->Branch("CVNScoreNC", &fCVNScoreNC);
  fTree->Branch("fCVNScoreProton0", &fCVNScoreProton0);
  fTree->Branch("fCVNScoreProton1", &fCVNScoreProton1);
  fTree->Branch("fCVNScoreProton2", &fCVNScoreProton2);
  fTree->Branch("fCVNScoreProton3", &fCVNScoreProton3);
  fTree->Branch("fCVNScorePion0", &fCVNScorePion0);
  fTree->Branch("fCVNScorePion1", &fCVNScorePion1);
  fTree->Branch("ErecNuMu", &fErecNuMu);
  fTree->Branch("ErecNuE", &fErecNuE);
  fTree->Branch("ErecNuMuRange", &fErecNuMuRange);
  fTree->Branch("ErecNuMuMCS", &fErecNuMuMCS);
  fTree->Branch("ErecNC", &fErecNC);
  fTree->Branch("VisibleEnergy", &fVisibleEnergy);
  fTree->Branch("DirectionRecNuE_x", &fDirectionRecNuE_x);
  fTree->Branch("DirectionRecNuE_y", &fDirectionRecNuE_y);
  fTree->Branch("DirectionRecNuE_z", &fDirectionRecNuE_z);
  fTree->Branch("DirectionRecNuMu_x", &fDirectionRecNuMu_x);
  fTree->Branch("DirectionRecNuMu_y", &fDirectionRecNuMu_y);
  fTree->Branch("DirectionRecNuMu_z", &fDirectionRecNuMu_z);
  fTree->Branch("DirectionRecNuEPfps_x", &fDirectionRecNuEPfps_x);
  fTree->Branch("DirectionRecNuEPfps_y", &fDirectionRecNuEPfps_y);
  fTree->Branch("DirectionRecNuEPfps_z", &fDirectionRecNuEPfps_z);
  fTree->Branch("DirectionRecNuMuPfps_x", &fDirectionRecNuMuPfps_x);
  fTree->Branch("DirectionRecNuMuPfps_y", &fDirectionRecNuMuPfps_y);
  fTree->Branch("DirectionRecNuMuPfps_z", &fDirectionRecNuMuPfps_z);
  fTree->Branch("DirectionRecHits_x", &fDirectionRecHits_x);
  fTree->Branch("DirectionRecHits_y", &fDirectionRecHits_y);
  fTree->Branch("DirectionRecHits_z", &fDirectionRecHits_z);
  fTree->Branch("IsCC", &fIsCC);
  fTree->Branch("InterMode", &fInterMode);
}
void test::atmoAnalysis::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(test::atmoAnalysis)
