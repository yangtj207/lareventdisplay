////////////////////////////////////////////////////////////////////////
// Class:       EventDisplay
// Module Type: analyzer
// File:        EventDisplay_module.cc
//
// Generated at Fri Aug 28 14:45:48 2015 by Matt Robinson using artmod
// from cetpkgsupport v1_08_06.
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
#include "canvas/Persistency/Common/Assns.h"

//root
#include "TDatabasePDG.h"
#include "THashList.h"
#include "TCollection.h"

//larsoft
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"


//geant4
#include "Geant4/G4UImanager.hh"
#include "Geant4/G4RunManager.hh"
#include "Geant4/G4VisExecutive.hh"
#include "Geant4/QGSP_BERT.hh"

//nusoft
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nutools/G4Base/DetectorConstruction.h"

//system
#include <iostream>
#include <fstream>
#include <sqlite3.h>
#include <zlib.h>

//local
#include "LandedSocket.h"

namespace evd {
  class Landed;
}

class evd::Landed : public art::EDAnalyzer {
public:
  explicit Landed(fhicl::ParameterSet const & p);

  // Plugins should not be copied or assigned.
  Landed(Landed const &) = delete;
  Landed(Landed &&) = delete;
  Landed & operator = (Landed const &) = delete;
  Landed & operator = (Landed &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;
  
  //optional functions
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.
  art::ServiceHandle<geo::Geometry>  geom_;
  art::ServiceHandle<LandedSocket> sock_;
  //string containing gdml having ben simplified through root
  std::string vrml_;
  std::string outputFilename_;
  std::string which_;
  sqlite3* ppDb_;
  std::vector<int> events_;
};


evd::Landed::Landed(fhicl::ParameterSet const & p):
  EDAnalyzer(p),
 // More initializers here.
  geom_(art::ServiceHandle<geo::Geometry>()),
  sock_(art::ServiceHandle<LandedSocket>()),
  outputFilename_(p.get<std::string>("outputFilename", "")),
  which_(p.get<std::string>("which", "truth")),
  ppDb_(nullptr),
  events_(p.get<std::vector<int> >("events", std::vector<int>()))
{
  
  if (outputFilename_.empty() && !sock_->connect())
    {
      throw cet::exception("Landed") << "no output filename specified and unable to connect to LANDED app using local socket\n" << 
	"please either specify a filename, or check LANDED is running and that you have started the server\n";	  

    }
  //  sock_->hello();

  std::string gdml=geom_->GDMLFile();

  // g4b::G4Helper* g4helper=new g4b::G4Helper();
  // g4helper->ConstructDetector(gdml_);
  // g4helper->InitPhysics();
  G4RunManager* runManager=new G4RunManager;
  runManager->SetUserInitialization(new g4b::DetectorConstruction(gdml));
  runManager->SetUserInitialization(new QGSP_BERT);
  // Visualization manager construction
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();

  char tmpdir[1024];
  strcpy(tmpdir, "/tmp/vcXXXXXX");
  mkdtemp(tmpdir);
  strcat(tmpdir, "/");
  char* vrmlviewer=getenv("G4VRMLFILE_VIEWER");
  if (vrmlviewer)
    unsetenv("G4VRMLFILE_VIEWER");
  char* oldvrmldir=getenv("G4VRMLFILE_DEST_DIR");
  setenv("G4VRMLFILE_DEST_DIR", tmpdir, 1);
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  UImanager->ApplyCommand("/run/initialize");
  UImanager->ApplyCommand("/vis/open VRML2FILE");
  UImanager->ApplyCommand("/vis/drawVolume");
  UImanager->ApplyCommand("/vis/viewer/flush");
  if (oldvrmldir)
    setenv("G4VRMLFILE_DEST_DIR", oldvrmldir, 1);
  else
    unsetenv("G4VRMLFILE_DEST_DIR");
  if (vrmlviewer)
    setenv("G4VRMLFILE_VIEWER", vrmlviewer, 1);

  char tmpname[1024];
  strcpy(tmpname, tmpdir);
  strcat(tmpname, "g4_00.wrl");
  std::ifstream geofile(tmpname);
  std::stringstream geodata;
  geodata << geofile.rdbuf();
  geofile.close();
  ::remove(tmpname);
  ::remove(tmpdir);
  vrml_=geodata.str();

  //  std::cout << "!!!!!!!  " << geodata.str().length() << std::endl;
  //  std::cout << geodata.str() << std::endl;

}

void evd::Landed::beginJob()
{

  if (!outputFilename_.empty())
    {
      //compress geometry
      uLongf complen=compressBound(vrml_.length());
      Bytef* comp=new Bytef[complen];
      if (compress(comp, &complen, (const Bytef*)vrml_.c_str(), vrml_.length())!=Z_OK)
	throw cet::exception("Landed") << "failed to compress geometry vrml\n";

      //create sqlite output file, compress and store gdml geometry
      //create sqlite output file
      ::remove(outputFilename_.c_str());
      if (sqlite3_open_v2(outputFilename_.c_str(), &ppDb_, SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE, "unix")!=SQLITE_OK)
	throw cet::exception("Landed") << "failed to create sqlite file\n";
      sqlite3_stmt* stmt;
      char * sErrMsg = 0;
      sqlite3_exec(ppDb_, "PRAGMA synchronous = OFF", NULL, NULL, &sErrMsg);
      
      
      //create geometry table
      if (sqlite3_prepare(ppDb_, "CREATE TABLE Geometry ( vrml BLOB, size INT );", 100, &stmt, nullptr)!=SQLITE_OK)
	throw cet::exception("Landed") << "failed to create geometry table on prepare\n";
      if (sqlite3_step(stmt)!=SQLITE_DONE)
	throw cet::exception("Landed") << "failed to create geometry table on step\n";
      if (sqlite3_finalize(stmt)!=SQLITE_OK)
	throw cet::exception("Landed") << "failed to create geometry table on finalize\n";
      //store geometry
      if (sqlite3_prepare(ppDb_, "INSERT INTO Geometry ( vrml, size ) VALUES ( ?, ? );", -1, &stmt, nullptr)!=SQLITE_OK)
	throw cet::exception("Landed") << "failed to fill geometry table on prepare\n";
      if (sqlite3_bind_blob(stmt, 1, comp, complen, SQLITE_STATIC)!=SQLITE_OK ||
	  sqlite3_bind_int(stmt, 2, vrml_.length())!=SQLITE_OK)
	throw cet::exception("Landed") << "failed to fill geometry table on bind\n";
      if (sqlite3_step(stmt)!=SQLITE_DONE)
	throw cet::exception("Landed") << "failed to fill geometry table on step\n";
      if (sqlite3_finalize(stmt)!=SQLITE_OK)
	throw cet::exception("Landed") << "failed to fill geometry table on finalize\n";
      delete[] comp;
      
      //create events table
      if (sqlite3_prepare(ppDb_, "CREATE TABLE Events ( RunNumber INT, SubRunNumber INT, EventNumber INT);", 200, &stmt, nullptr)!=SQLITE_OK)
	throw cet::exception("Landed") << "failed to create events table on prepare\n";
      if (sqlite3_step(stmt)!=SQLITE_DONE)
	throw cet::exception("Landed") << "failed to create events table on step\n";
      if (sqlite3_finalize(stmt)!=SQLITE_OK)
	throw cet::exception("Landed") << "failed to create events table on finalize\n";

      //create simchannel ide table
      if (sqlite3_prepare(ppDb_, "CREATE TABLE Hits ( RunNumber INT, SubRunNumber INT, EventNumber INT, NumElectrons REAL, Energy REAL, X REAL, Y REAL, Z REAL, Track INT);", 200, &stmt, nullptr)!=SQLITE_OK)
	throw cet::exception("Landed") << "failed to create hits table on prepare\n";
      if (sqlite3_step(stmt)!=SQLITE_DONE)
	throw cet::exception("Landed") << "failed to create hits table on step\n";
      if (sqlite3_finalize(stmt)!=SQLITE_OK)
	throw cet::exception("Landed") << "failed to create hits table on finalize\n";

      //create mcparticle table
      if (sqlite3_prepare(ppDb_, "CREATE TABLE Tracks (RunNumber INT, SubRunNumber INT, EventNumber INT, ID INT, PDGCode INT, Segment INT, Energy REAL, X REAL, Y REAL, Z REAL);", 200, &stmt, nullptr)!=SQLITE_OK)
	throw cet::exception("Landed") << "failed to create tracks table on prepare\n";
      if (sqlite3_step(stmt)!=SQLITE_DONE)
	throw cet::exception("Landed") << "failed to create tracks table on step\n";
      if (sqlite3_finalize(stmt)!=SQLITE_OK)
	throw cet::exception("Landed") << "failed to create tracks table on finalize\n";
    }
  else
    {
      //send geometry through landed service
      sock_->sendGeometry(vrml_);
    }

}
void evd::Landed::endJob()
{
  sqlite3_close(ppDb_);
}
void evd::Landed::analyze(art::Event const & event)
{
  if (!events_.size() || outputFilename_.empty() || std::find(events_.begin(), events_.end(), event.event())!=events_.end())
    {
      sqlite3_stmt* stmt;
      char * sErrMsg = 0;

      //mc truth
      art::Handle<std::vector<simb::MCParticle> > mcparticles;
      art::Handle<std::vector<sim::SimChannel> > simchannels;
      //reconstruction
      art::Handle<std::vector<anab::Calorimetry> > calorimetry;
      //art::Handle<std::vector<recob::Track> > tracks;
      //art::Handle<std::vector<recob::PFParticle> > particles;
      art::Handle<art::Assns<recob::Track,anab::Calorimetry,void> > track2calo;
      art::Handle<art::Assns<recob::PFParticle,recob::Track,void> > particle2track;
      std::map<size_t, int> pdgcodes;

      //get unique event id numbers  
      art::RunNumber_t runnum=event.run();
      art::SubRunNumber_t subrunnum=event.subRun();
      art::EventNumber_t eventnum=event.event();

      std::string whichcalo;
      if (which_.substr(which_.length()-2, 2)=="dc")
	whichcalo=which_.substr(0, which_.length()-2)+std::string("calodc");
      else
	whichcalo=which_+std::string("calo");
  
      if (which_=="truth")
	{
	  event.getByLabel("largeant", mcparticles);
	  event.getByLabel("largeant", simchannels);
	}
      else
	{
	  event.getByLabel(whichcalo, calorimetry);
	  //event.getByLabel(which_, tracks);
	  //event.getByLabel(which_, particles);
	  event.getByLabel(whichcalo, track2calo);
	  event.getByLabel(which_, particle2track);

	  for (size_t i=0; i<particle2track->size(); i++)
	    {
	      //int pdgcode=particle2track->at(i).first->PdgCode();
	      for (size_t j=0; j<track2calo->size(); j++)
		{
		  if (track2calo->at(j).first==particle2track->at(i).second)
		    {
		      for (size_t k=0; k<calorimetry->size(); k++)
			{
			  if (track2calo->at(j).second==art::Ptr<anab::Calorimetry>(calorimetry, k))
			    {
			      pdgcodes[k]=particle2track->at(i).first->PdgCode();
			    }
			}
		    }
		}
	    }
	}
      //  //determine available reconstructions
      //  std::cout << "looking for reco space point sets\n";
      //  std::vector<art::Handle<std::vector<recob::SpacePoint> > > recohits;
      //  event.getManyByType(recohits);
      //  std::cout << "reco hit sets: " << recohits.size() << std::endl;
      //  for (size_t i=0; i<recohits.size(); i++)
      //    {
      //      const art::Provenance* prov=recohits[i].provenance();
      //      std::cout << i << ": ";
      //      std::cout << prov->moduleLabel();
      //      std::cout << std::endl;
      //    }
      // std::cout << "looking for reco track sets\n";
      //  std::vector<art::Handle<std::vector<recob::Track> > > recotracks;
      //  event.getManyByType(recotracks);
      //  std::cout << "reco track sets: " << recotracks.size() << std::endl;
      //  for (size_t i=0; i<recotracks.size(); i++)
      //    {
      //      const art::Provenance* prov=recotracks[i].provenance();
      //      std::cout << i << ": ";
      //      std::cout << prov->moduleLabel();
      //      std::cout << std::endl;
      //    }
 
      //  std::cout << "Event: " << event.id() << std::endl;
  
      int hitcount=0;
      int trajtotal=0;
      if (which_=="truth")
	{
	  //store sim channels and their ide's
	  for (auto const& simchannel:*simchannels)
	    if (geom_->SignalType(simchannel.Channel()) == geo::kCollection)
	      //ides
	      for (auto const& tdcide:simchannel.TDCIDEMap())
		hitcount+=tdcide.second.size();
	  if (mcparticles->size())
	    {
	      //put first mctruth particle into database
	      for (auto const& part:*mcparticles)
		{
		  for (size_t i=0; i<part.NumberTrajectoryPoints(); i++)
		    {
		      trajtotal++;
		    }
		}
	    }
	}
      else
	{
	  for (size_t tracknum=0; tracknum<calorimetry->size(); tracknum++)
	    {
	      hitcount+=calorimetry->at(tracknum).XYZ().size();
	    }
	}

      if (!outputFilename_.empty())
	{
	  if (sqlite3_prepare(ppDb_, "INSERT INTO Events (RunNumber, SubRunNumber, EventNumber) VALUES ( ?, ?, ? );", 200, &stmt, nullptr)!=SQLITE_OK)
	    throw cet::exception("Landed") << "failed to fill events table on prepare\n";
	  if (sqlite3_bind_int(stmt, 1, runnum)!=SQLITE_OK)
	    throw cet::exception("Landed") << "failed to fill events table on bind\n";
	  if (sqlite3_bind_int(stmt, 2, subrunnum)!=SQLITE_OK)
	    throw cet::exception("Landed") << "failed to fill events table on bind\n";
	  if (sqlite3_bind_int(stmt, 3, eventnum)!=SQLITE_OK)
	    throw cet::exception("Landed") << "failed to fill events table on bind\n";
	  if (sqlite3_step(stmt)!=SQLITE_DONE)
	    throw cet::exception("Landed") << "failed to fill events table on step\n";
	  if (sqlite3_finalize(stmt)!=SQLITE_OK)
	    throw cet::exception("Landed") << "failed to fill events table on finalize\n";

	  if (which_=="truth")
	    {
	      //store sim channels and their ide's
	      for (auto const& simchannel:*simchannels)
		if (geom_->SignalType(simchannel.Channel()) == geo::kCollection)
		  {
		    //ides
		    sqlite3_exec(ppDb_, "BEGIN TRANSACTION", NULL, NULL, &sErrMsg);
		    for (auto const& tdcide:simchannel.TDCIDEMap())
		      {
			//	      unsigned int tdc=tdcide.first;
			auto const& idevec=tdcide.second;
			for (auto const& ide:idevec)
			  {
			    // SimChannelID INT, TDC INT, NumElectrons REAL, Energy REAL, X REAL, Y REAL, Z REAL
			    if (sqlite3_prepare(ppDb_, "INSERT INTO Hits (RunNumber, SubRunNumber, EventNumber, NumElectrons, Energy, X, Y, Z, Track ) VALUES ( ?, ?, ?, ?, ?, ?, ?, ?, ? );", 200, &stmt, nullptr)!=SQLITE_OK)
			      throw cet::exception("Landed") << "failed to fill hits table on prepare\n";
			    if (sqlite3_bind_int(stmt, 1, runnum)!=SQLITE_OK)
			      throw cet::exception("Landed") << "failed to fill hits table on bind\n";
			    if (sqlite3_bind_int(stmt, 2, subrunnum)!=SQLITE_OK)
			      throw cet::exception("Landed") << "failed to fill hits table on bind\n";
			    if (sqlite3_bind_int(stmt, 3, eventnum)!=SQLITE_OK)
			      throw cet::exception("Landed") << "failed to fill hits table on bind\n";
			    if (sqlite3_bind_double(stmt, 4, ide.numElectrons)!=SQLITE_OK)
			      throw cet::exception("Landed") << "failed to fill hits table on bind\n";
			    if (sqlite3_bind_double(stmt, 5, ide.energy)!=SQLITE_OK)
			      throw cet::exception("Landed") << "failed to fill hits table on bind\n";
			    if (sqlite3_bind_double(stmt, 6, ide.x)!=SQLITE_OK)
			      throw cet::exception("Landed") << "failed to fill hits table on bind\n";
			    if (sqlite3_bind_double(stmt, 7, ide.y)!=SQLITE_OK)
			      throw cet::exception("Landed") << "failed to fill hits table on bind\n";
			    if (sqlite3_bind_double(stmt, 8, ide.z)!=SQLITE_OK)
			      throw cet::exception("Landed") << "failed to fill hits table on bind\n";
			    if (sqlite3_bind_int(stmt, 9, ide.trackID)!=SQLITE_OK)
			      throw cet::exception("Landed") << "failed to fill hits table on bind\n";		  
			    if (sqlite3_step(stmt)!=SQLITE_DONE)
			      throw cet::exception("Landed") << "failed to fill hits table on step\n";
			    if (sqlite3_finalize(stmt)!=SQLITE_OK)
			      throw cet::exception("Landed") << "failed to fill hits table on finalize\n";

			  }
		      }
		    sqlite3_exec(ppDb_, "END TRANSACTION", NULL, NULL, &sErrMsg);
		  }

	      //store mctruth tracks
	      sqlite3_exec(ppDb_, "BEGIN TRANSACTION", NULL, NULL, &sErrMsg);
	      for (auto const& part:*mcparticles)
		{
		  for (size_t i=0; i<part.NumberTrajectoryPoints(); i++)
		    {
		      if (sqlite3_prepare(ppDb_, "INSERT INTO Tracks ( RunNumber, SubRunNumber, EventNumber, ID, PDGCode, Segment, Energy, X, Y, Z) VALUES ( ?, ?, ?, ?, ?, ?, ?, ?, ?, ? );", 200, &stmt, nullptr)!=SQLITE_OK)
			throw cet::exception("Landed") << "failed to fill tracks table on prepare\n";
		      if (sqlite3_bind_int(stmt, 1, runnum)!=SQLITE_OK)
			throw cet::exception("Landed") << "failed to fill tracks table on bind\n";
		      if (sqlite3_bind_int(stmt, 2, subrunnum)!=SQLITE_OK)
			throw cet::exception("Landed") << "failed to fill tracks table on bind\n";
		      if (sqlite3_bind_int(stmt, 3, eventnum)!=SQLITE_OK)
			throw cet::exception("Landed") << "failed to fill tracks table on bind\n";
		      if (sqlite3_bind_int(stmt, 4, part.TrackId())!=SQLITE_OK)
			throw cet::exception("Landed") << "failed to fill tracks table on bind\n";
		      if (sqlite3_bind_int(stmt, 5, part.PdgCode())!=SQLITE_OK)
			throw cet::exception("Landed") << "failed to fill tracks table on bind\n";
		      if (sqlite3_bind_int(stmt, 6, i)!=SQLITE_OK)
			throw cet::exception("Landed") << "failed to fill tracks table on bind\n";
		      if (sqlite3_bind_double(stmt, 7, part.E(i))!=SQLITE_OK)
			throw cet::exception("Landed") << "failed to fill tracks table on bind\n";
		      if (sqlite3_bind_double(stmt, 8, part.Vx(i))!=SQLITE_OK)
			throw cet::exception("Landed") << "failed to fill tracks table on bind\n";
		      if (sqlite3_bind_double(stmt, 9, part.Vy(i))!=SQLITE_OK)
			throw cet::exception("Landed") << "failed to fill tracks table on bind\n";
		      if (sqlite3_bind_double(stmt, 10, part.Vz(i))!=SQLITE_OK)
			throw cet::exception("Landed") << "failed to fill tracks table on bind\n";
		      if (sqlite3_step(stmt)!=SQLITE_DONE)
			throw cet::exception("Landed") << "failed to fill tracks table on step\n";
		      if (sqlite3_finalize(stmt)!=SQLITE_OK)
			throw cet::exception("Landed") << "failed to fill tracks table on finalize\n";
		    }
		}
	      sqlite3_exec(ppDb_, "END TRANSACTION", NULL, NULL, &sErrMsg);
	    }
	  //reconstruction rather than truth
	  else
	    {
	      sqlite3_exec(ppDb_, "BEGIN TRANSACTION", NULL, NULL, &sErrMsg);  
	  
	      for (size_t tracknum=0; tracknum<calorimetry->size(); tracknum++)
		{
		  auto const& calo=calorimetry->at(tracknum);
		  for (size_t hitnum=0; hitnum<calo.dEdx().size() && hitnum<calo.dQdx().size() && hitnum<calo.XYZ().size(); hitnum++)
		    {
		      if (sqlite3_prepare(ppDb_, "INSERT INTO Hits (RunNumber, SubRunNumber, EventNumber, NumElectrons, Energy, X, Y, Z, Track ) VALUES ( ?, ?, ?, ?, ?, ?, ?, ?, ? );", 200, &stmt, nullptr)!=SQLITE_OK)
			throw cet::exception("Landed") << "failed to fill hits table on prepare\n";
		      if (sqlite3_bind_int(stmt, 1, runnum)!=SQLITE_OK)
			throw cet::exception("Landed") << "failed to fill hits table on bind\n";
		      if (sqlite3_bind_int(stmt, 2, subrunnum)!=SQLITE_OK)
			throw cet::exception("Landed") << "failed to fill hits table on bind\n";
		      if (sqlite3_bind_int(stmt, 3, eventnum)!=SQLITE_OK)
			throw cet::exception("Landed") << "failed to fill hits table on bind\n";
		      if (sqlite3_bind_double(stmt, 4, calo.TrkPitchVec().at(hitnum)*calo.dQdx().at(hitnum))!=SQLITE_OK)
			throw cet::exception("Landed") << "failed to fill hits table on bind\n";
		      if (sqlite3_bind_double(stmt, 5, calo.TrkPitchVec().at(hitnum)*calo.dEdx().at(hitnum))!=SQLITE_OK)
			throw cet::exception("Landed") << "failed to fill hits table on bind\n";
		      if (sqlite3_bind_double(stmt, 6, calo.XYZ().at(hitnum).X())!=SQLITE_OK)
			throw cet::exception("Landed") << "failed to fill hits table on bind\n";
		      if (sqlite3_bind_double(stmt, 7, calo.XYZ().at(hitnum).Y())!=SQLITE_OK)
			throw cet::exception("Landed") << "failed to fill hits table on bind\n";
		      if (sqlite3_bind_double(stmt, 8, calo.XYZ().at(hitnum).Z())!=SQLITE_OK)
			throw cet::exception("Landed") << "failed to fill hits table on bind\n";
		      if (sqlite3_bind_int(stmt, 9, tracknum+1)!=SQLITE_OK)
			throw cet::exception("Landed") << "failed to fill hits table on bind\n";
		      if (sqlite3_step(stmt)!=SQLITE_DONE)
			throw cet::exception("Landed") << "failed to fill hits table on step\n";
		      if (sqlite3_finalize(stmt)!=SQLITE_OK)
			throw cet::exception("Landed") << "failed to fill hits table on finalize\n";
		    }
		}
	      sqlite3_exec(ppDb_, "END TRANSACTION", NULL, NULL, &sErrMsg);
	  
	      sqlite3_exec(ppDb_, "BEGIN TRANSACTION", NULL, NULL, &sErrMsg);  
	      for (size_t tracknum=0; tracknum<calorimetry->size(); tracknum++)
		{
		  auto const& calo=calorimetry->at(tracknum);
		  if (calo.dEdx().size()>1)
		    {
		      int pdgcode=0;
		      if (pdgcodes.find(tracknum)!=pdgcodes.end())
			pdgcode=pdgcodes[tracknum];
		      double energy=calo.KineticEnergy();
		      for (size_t hitnum=0; hitnum<calo.dEdx().size() && hitnum<calo.dQdx().size() && hitnum<calo.XYZ().size(); hitnum++)
			{
			  if (sqlite3_prepare(ppDb_, "INSERT INTO Tracks ( RunNumber, SubRunNumber, EventNumber, ID, PDGCode, Segment, Energy, X, Y, Z) VALUES ( ?, ?, ?, ?, ?, ?, ?, ?, ?, ? );", 200, &stmt, nullptr)!=SQLITE_OK)
			    throw cet::exception("Landed") << "failed to fill tracks table on prepare\n";
			  if (sqlite3_bind_int(stmt, 1, runnum)!=SQLITE_OK)
			    throw cet::exception("Landed") << "failed to fill tracks table on bind\n";
			  if (sqlite3_bind_int(stmt, 2, subrunnum)!=SQLITE_OK)
			    throw cet::exception("Landed") << "failed to fill tracks table on bind\n";
			  if (sqlite3_bind_int(stmt, 3, eventnum)!=SQLITE_OK)
			    throw cet::exception("Landed") << "failed to fill tracks table on bind\n";
			  if (sqlite3_bind_int(stmt, 4, tracknum+1)!=SQLITE_OK)
			    throw cet::exception("Landed") << "failed to fill tracks table on bind\n";
			  if (sqlite3_bind_int(stmt, 5, pdgcode)!=SQLITE_OK)
			    throw cet::exception("Landed") << "failed to fill tracks table on bind\n";
			  if (sqlite3_bind_int(stmt, 6, hitnum+1)!=SQLITE_OK)
			    throw cet::exception("Landed") << "failed to fill tracks table on bind\n";
			  if (sqlite3_bind_double(stmt, 7, energy)!=SQLITE_OK)
			    throw cet::exception("Landed") << "failed to fill tracks table on bind\n";
			  if (sqlite3_bind_double(stmt, 8, calo.XYZ().at(hitnum).X())!=SQLITE_OK)
			    throw cet::exception("Landed") << "failed to fill tracks table on bind\n";
			  if (sqlite3_bind_double(stmt, 9, calo.XYZ().at(hitnum).Y())!=SQLITE_OK)
			    throw cet::exception("Landed") << "failed to fill tracks table on bind\n";
			  if (sqlite3_bind_double(stmt, 10, calo.XYZ().at(hitnum).Z())!=SQLITE_OK)
			    throw cet::exception("Landed") << "failed to fill tracks table on bind\n";
			  if (sqlite3_step(stmt)!=SQLITE_DONE)
			    throw cet::exception("Landed") << "failed to fill tracks table on step\n";
			  if (sqlite3_finalize(stmt)!=SQLITE_OK)
			    throw cet::exception("Landed") << "failed to fill tracks table on finalize\n"; 
			  energy-=calo.TrkPitchVec().at(hitnum)*calo.dEdx().at(hitnum);
			}
		    }
		}
	    }
	}	
      else
	{
	  sock_->sendEvent(hitcount, trajtotal, runnum, subrunnum, eventnum);
	  if (which_=="truth")
	    {

	      //store sim channels and their ide's
	      for (auto const& simchannel:*simchannels)
		if (geom_->SignalType(simchannel.Channel()) == geo::kCollection)
		  //ides
		  for (auto const& tdcide:simchannel.TDCIDEMap())
		    {
		      auto const& idevec=tdcide.second;
		      for (auto const& ide:idevec)
			{
			  sock_->sendHit(ide.x, ide.y, ide.z, ide.energy, ide.numElectrons, ide.trackID);
			}
		    }
	      //store mctruth tracks
	      for (auto const& part:*mcparticles)
		{
		  for (size_t i=0; i<part.NumberTrajectoryPoints(); i++)
		    {
		      sock_->sendVertex(part.Vx(i), part.Vy(i), part.Vz(i), part.E(i), part.TrackId(), part.PdgCode());
		    }
		}
	    }
	  //reco rather than truth
	  else
	    {
	      for (size_t tracknum=0; tracknum<calorimetry->size(); tracknum++)
		{
		  auto const& calo=calorimetry->at(tracknum);
		  for (size_t hitnum=0; hitnum<calo.dEdx().size() && hitnum<calo.dQdx().size() && hitnum<calo.XYZ().size(); hitnum++)
		    {
		      sock_->sendHit(calo.XYZ().at(hitnum).X(), calo.XYZ().at(hitnum).Y(), calo.XYZ().at(hitnum).Z(), calo.TrkPitchVec().at(hitnum)*calo.dEdx().at(hitnum), calo.TrkPitchVec().at(hitnum)*calo.dQdx().at(hitnum), tracknum+1);
		    }
		}
	      for (size_t tracknum=0; tracknum<calorimetry->size(); tracknum++)
		{
		  int pdgcode=0;
		  if (pdgcodes.find(tracknum)!=pdgcodes.end())
		    pdgcode=pdgcodes[tracknum];
		  auto const& calo=calorimetry->at(tracknum);
		  if (calo.dEdx().size()>1)
		    {
		      double energy=calo.KineticEnergy();
		      for (size_t hitnum=0; hitnum<calo.dEdx().size() && hitnum<calo.dQdx().size() && hitnum<calo.XYZ().size(); hitnum++)
			{
			  sock_->sendVertex(calo.XYZ().at(hitnum).X(), calo.XYZ().at(hitnum).Y(), calo.XYZ().at(hitnum).Z(), energy, tracknum+1, pdgcode);
			  energy-=calo.TrkPitchVec().at(hitnum)*calo.dEdx().at(hitnum);
			}
		    }
		}
	    }
	}
    }
}
DEFINE_ART_MODULE(evd::Landed)
