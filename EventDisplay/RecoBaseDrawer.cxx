/// \file    RecoBaseDrawer.cxx
/// \brief   Class to aid in the rendering of RecoBase objects
/// \author  brebel@fnal.gov
/// \version $Id: RecoBaseDrawer.cxx,v 1.3 2010/11/11 22:47:19 p-novaart Exp $
#include <cmath>
#include <map>
#include <stdint.h>

#include "TMath.h"
#include "TMarker.h"
#include "TBox.h"
#include "TH1.h"
#include "TLine.h"
#include "TPolyLine.h"
#include "TPolyLine3D.h"
#include "TPolyMarker.h"
#include "TPolyMarker3D.h"
#include "TVector3.h"
#include "TText.h"

#include "EventDisplayBase/View2D.h"
#include "EventDisplayBase/View3D.h"
#include "EventDisplayBase/EventHolder.h"
#include "EventDisplay/eventdisplay.h"
#include "EventDisplay/RecoBaseDrawer.h"
#include "EventDisplay/RecoDrawingOptions.h"
#include "EventDisplay/ColorDrawingOptions.h"
#include "EventDisplay/RawDrawingOptions.h"
#include "RecoBase/Wire.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/SpacePoint.h"
#include "RecoBase/Track.h"
#include "RecoBase/Shower.h"
#include "RecoBase/Event.h"
#include "RecoBase/EndPoint2D.h"
#include "RecoBase/Seed.h"
#include "RecoObjects/BezierTrack.h"
#include "RecoBase/Vertex.h"
#include "RecoBase/OpFlash.h"
#include "Filters/ChannelFilter.h"
#include "Geometry/Geometry.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/WireGeo.h"
#include "Utilities/LArProperties.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/AssociationUtil.h"

#include "cetlib/exception.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Event.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Framework/Principal/Handle.h" 
#include "art/Framework/Core/FindMany.h" 
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace {
  // Utility function to make uniform error messages.
  void writeErrMsg(const char* fcn,
		   cet::exception const& e)
  {
    mf::LogWarning("RecoBaseDrawer") << "RecoBaseDrawer::" << fcn
				     << " failed with message:\n"
				     << e;
  }
}

namespace evd{

  //......................................................................
  RecoBaseDrawer::RecoBaseDrawer() 
  {
    art::ServiceHandle<geo::Geometry>            geo;
    art::ServiceHandle<util::LArProperties>      larp;
    art::ServiceHandle<util::DetectorProperties> detp;
    art::ServiceHandle<evd::RawDrawingOptions>   rawopt;

    fWireMin.resize(0);   
    fWireMax.resize(0);    
    fTimeMin.resize(0);    
    fTimeMax.resize(0);    
    fRawCharge.resize(0);   
    fConvertedCharge.resize(0);
    
    // set the list of bad channels in this detector
    filter::ChannelFilter cf;
    for(size_t t = 0; t < geo->NTPC(); ++t){
      unsigned int nplanes = geo->Nplanes(t);
      fWireMin.resize(nplanes,-1);   
      fWireMax.resize(nplanes,-1);    
      fTimeMin.resize(nplanes,-1);    
      fTimeMax.resize(nplanes,-1);    
      fRawCharge.resize(nplanes,0);   
      fConvertedCharge.resize(nplanes,0);
      for(size_t p = 0; p < geo->Nplanes(t); ++p){
	fWireMin[p] = 0;
	fWireMax[p] = geo->TPC(t).Plane(p).Nwires();
	fTimeMin[p] = 0;
	fTimeMax[p] = rawopt->fTicks;
	for(size_t w = 0; w < geo->TPC(t).Plane(p).Nwires(); ++w){
	  uint32_t channel = geo->PlaneWireToChannel(p, w, t);
	  if( cf.BadChannel(channel) ) fBadChannels.push_back(channel);
	} // end loop over wires
      }// end loop over planes
    }// end loop over TPCs
    
  }

  //......................................................................
  RecoBaseDrawer::~RecoBaseDrawer() 
  {
 
  }

  //......................................................................
  void RecoBaseDrawer::Wire2D(const art::Event& evt,
			      evdb::View2D*     view,
			      unsigned int      plane)
  {
    art::ServiceHandle<evd::RawDrawingOptions>   rawOpt;
    art::ServiceHandle<evd::RecoDrawingOptions>  recoOpt;
    art::ServiceHandle<geo::Geometry>            geo;
    art::ServiceHandle<evd::ColorDrawingOptions> cst;

    if(rawOpt->fDrawRawDataOrCalibWires < 1)    return;
  
    int ticksPerPoint = rawOpt->fTicksPerPoint;
    int ticks         = 0; // fill this below

    // to make det independent later:
    double mint = 5000;
    double maxt = 0;
    double minw = 5000;
    double maxw = 0;

    geo::PlaneID pid(rawOpt->fCryostat, rawOpt->fTPC, plane);
    
    for(size_t imod = 0; imod < recoOpt->fWireLabels.size(); ++imod) {
      std::string const which = recoOpt->fWireLabels[imod];
    
      art::PtrVector<recob::Wire> wires;
      this->GetWires(evt, which, wires);

      if(wires.size() < 1) return;

      ticks = wires[0]->NSignal();
      
      for(size_t i = 0; i < wires.size(); ++i) {
      
	uint32_t channel = wires[i]->Channel();
	std::vector<geo::WireID> wireids = geo->ChannelToWire(channel);
      
	geo::SigType_t sigType = geo->SignalType(channel);

	for (auto const& wid : wireids){
	  if (wid.planeID() != pid) continue;

	  double wire = 1.*wid.Wire;
	  double tick = 0;
	  // get the unpacked ROIs
	  std::vector<float> wirSig = wires[i]->Signal();
	  if(wirSig.size() == 0) continue;
	  // get an iterator over the adc values
	  std::vector<float>::const_iterator itr = wirSig.begin();
	  while( itr != wirSig.end() ){
	    int ticksUsed = 0;
	    double tdcsum = 0.;
	    double adcsum = 0.;
	    while(ticksUsed < ticksPerPoint && itr != wirSig.end()){
	      tdcsum  += tick;
	      adcsum  += (1.*(*itr));
	      ++ticksUsed;
	      tick += 1.;
	      itr++; // this advance of the iterator is sufficient for the external loop too
	    }
	    double adc = adcsum/ticksPerPoint;
	    double tdc = tdcsum/ticksPerPoint;
	    
	    if(TMath::Abs(adc) < rawOpt->fMinSignal) continue;
	    
	    int    co = 0;
	    double sf = 1.;
	    double q0 = 1000.0;
	    
	    co = cst->CalQ(sigType).GetColor(adc);
	    if (rawOpt->fScaleDigitsByCharge) {
	      sf = sqrt(adc/q0);
	      if (sf>1.0) sf = 1.0;
	    }
	    
	    if(wire < minw) minw = wire;
	    if(wire > maxw) maxw = wire;  
	    if(tdc  < mint) mint = tdc;
	    if(tdc  > maxt) maxt = tdc;
	    
	    if(rawOpt->fAxisOrientation < 1){	
	      TBox& b1 = view->AddBox(wire-sf*0.5,
				      tdc-sf*0.5*ticksPerPoint,
				      wire+sf*0.5,
				      tdc+sf*0.5*ticksPerPoint);
	      b1.SetFillStyle(1001);
	      b1.SetFillColor(co);    
	      b1.SetBit(kCannotPick);
	    }
	    else{
	      TBox& b1 = view->AddBox(tdc-sf*0.5*ticksPerPoint,
				      wire-sf*0.5,
				      tdc+sf*0.5*ticksPerPoint,
				      wire+sf*0.5);
	      b1.SetFillStyle(1001);
	      b1.SetFillColor(co);    
	      b1.SetBit(kCannotPick);
	    }
	  }// end loop over samples 
	}//end loop over wire segments
      }//end loop over wires
    }// end loop over wire module labels

    fWireMin[plane] = minw;   
    fWireMax[plane] = maxw;    
    fTimeMin[plane] = mint;    
    fTimeMax[plane] = maxt; 
    
    // now loop over all the bad channels and set them to 0 adc
    for(size_t bc = 0; bc < fBadChannels.size(); ++bc){

      geo::SigType_t sigType = geo->SignalType(fBadChannels[bc]);

      std::vector<geo::WireID> wireids = geo->ChannelToWire(fBadChannels[bc]);

      for( auto const& wid : wireids){
	if(wid.planeID() != pid) continue;

	if(rawOpt->fMinSignal > 0) continue;
	
	
	int    co   = cst->CalQ(sigType).GetColor(0);
	double wire = 1.*wid.Wire;
	
	for(int i = 0; i < ticks; i += ticksPerPoint){
	  double tdc = i + 0.5*ticksPerPoint;
	  
	  if(rawOpt->fAxisOrientation < 1){	
	    TBox& b1 = view->AddBox(wire-0.5,
				    tdc-0.5*ticksPerPoint,
				    wire+0.5,
				    tdc+0.5*ticksPerPoint);
	    b1.SetFillStyle(1001);
	    b1.SetFillColor(co);    
	    b1.SetBit(kCannotPick);
	  }
	  else{
	    TBox& b1 = view->AddBox(tdc-0.5*ticksPerPoint,
				    wire-0.5,
				    tdc+0.5*ticksPerPoint,
				    wire+0.5);
	    b1.SetFillStyle(1001);
	    b1.SetFillColor(co);    
	    b1.SetBit(kCannotPick);
	  }
	}
      }
    }// end loop over bad channels    

    return;
  }

  //......................................................................
  ///
  /// Render Hit objects on a 2D viewing canvas
  ///
  /// @param evt    : Event handle to get data objects from
  /// @param view   : Pointer to view to draw on
  /// @param plane  : plane number of view
  ///
  void RecoBaseDrawer::Hit2D(const art::Event& evt,
			     evdb::View2D*     view,
			     unsigned int      plane)
  {
    art::ServiceHandle<evd::RawDrawingOptions>   rawOpt;
    art::ServiceHandle<evd::RecoDrawingOptions>  recoOpt;
    art::ServiceHandle<geo::Geometry>            geo;
    art::ServiceHandle<util::LArProperties>      larp;
    art::ServiceHandle<util::DetectorProperties> detp;
    
    if (recoOpt->fDrawHits == 0)                return;
    if (rawOpt->fDrawRawDataOrCalibWires < 1)   return;
    
    fRawCharge[plane]       = 0;
    fConvertedCharge[plane] = 0;

    for (size_t imod = 0; imod < recoOpt->fHitLabels.size(); ++imod) {
      std::string const which = recoOpt->fHitLabels[imod];
    
      std::vector<const recob::Hit*> hits;
      this->GetHits(evt, which, hits, plane);
  
      // Display all hits on the two 2D views provided
      for(auto itr : hits){
	// Try to get the "best" charge measurement, ie. the one last in
	// the calibration chain      
	fRawCharge[itr->WireID().Plane]    += itr->PeakAmplitude();
	double dQdX = itr->PeakAmplitude()/geo->WirePitch()/detp->ElectronsToADC();
	fConvertedCharge[itr->WireID().Plane] += larp->BirksCorrection(dQdX);	
      } // loop on hits

      this->Hit2D(hits, kBlack, view);

    } // loop on imod folders
      
    return;
  }

  //......................................................................
  ///
  /// Render Hit objects on a 2D viewing canvas
  ///
  /// @param hits   : vector of hits for the veiw
  /// @param color  : color of associated cluster/prong
  /// @param view   : Pointer to view to draw on
  ///
  /// assumes the hits are all from the correct plane for the given view
  void RecoBaseDrawer::Hit2D(std::vector<const recob::Hit*> hits,
			     int                            color,
			     evdb::View2D*                  view,
                             bool              drawConnectingLines) 
  {
    art::ServiceHandle<evd::RawDrawingOptions>   rawOpt;
    art::ServiceHandle<evd::RecoDrawingOptions>  recoOpt;
    art::ServiceHandle<geo::Geometry>            geo;

    unsigned int w  = 0;
    unsigned int wold = 0;
    float timeold = 0.;
    
    if(color==-1)
      color=recoOpt->fSelectedHitColor;
    
    for(unsigned int c = 0; c < hits.size(); ++c){

      // check that we are in the correct TPC
      // the view should tell use we are in the correct plane
      if(hits[c]->WireID().TPC      != rawOpt->fTPC || 
	 hits[c]->WireID().Cryostat != rawOpt->fCryostat) continue;
 
      w = hits[c]->WireID().Wire;

      // Try to get the "best" charge measurement, ie. the one last in
      // the calibration chain
      float time = hits[c]->PeakTime();

      if(rawOpt->fAxisOrientation < 1){
	TBox& b1 = view->AddBox(w-0.5, time-0.5, w+0.5, time+0.5);
        if(drawConnectingLines && c > 0) {
          TLine& l = view->AddLine(w, time, wold, timeold);
          l.SetLineColor(color);
	  l.SetBit(kCannotPick);
        }
	b1.SetFillStyle(0);
	b1.SetLineColor(color);
	b1.SetBit(kCannotPick);
      }
      else{
	TBox& b1 = view->AddBox(time-0.5, w-0.5, time+0.5, w+0.5);
        if(drawConnectingLines && c > 0) {
          TLine& l = view->AddLine(time, w, timeold, wold);
          l.SetLineColor(color);
	  l.SetBit(kCannotPick);
        }
	b1.SetFillStyle(0);
	b1.SetLineColor(color);
	b1.SetBit(kCannotPick);
      }
      wold = w;
      timeold = time;

    } // loop on hits
    
    return;
  }

  //........................................................................
  int RecoBaseDrawer::GetRegionOfInterest(int plane,
					  int& minw,
					  int& maxw,
					  int& mint,
					  int& maxt)
  {
    art::ServiceHandle<evd::RawDrawingOptions>   rawOpt;
    art::ServiceHandle<geo::Geometry>            geo;

    if((unsigned int)plane > fWireMin.size()){
      mf::LogWarning  ("RecoBaseDrawer") << " Requested plane " 
					 << plane 
					 << " is larger than those available ";
      return -1;
    }
  
    minw = fWireMin[plane];
    maxw = fWireMax[plane];
    mint = fTimeMin[plane];
    maxt = fTimeMax[plane];
  
    //make values a bit larger, but make sure they don't go out of bounds 
    minw = (minw-30<0) ? 0 : minw-30;
    mint = (mint-10<0) ? 0 : mint-10;

    int fTicks = rawOpt->fTicks;
  
    maxw = (maxw+10 > (int)geo->Nwires(plane)) ? geo->Nwires(plane) : maxw+10;
    maxt = (maxt+10 > fTicks) ? fTicks : maxt+10;
    
    return 0;
  }

  //......................................................................
  void RecoBaseDrawer::GetChargeSum(int plane,
				    double& charge,
				    double& convcharge)
  {
   charge     = fRawCharge[plane]; 
   convcharge = fConvertedCharge[plane];        
   
   return;
 }

  //......................................................................
  void RecoBaseDrawer::EndPoint2D(const art::Event& evt,
				  evdb::View2D*     view,
				  unsigned int      plane) 
  {
    art::ServiceHandle<evd::RawDrawingOptions>   rawOpt;
    art::ServiceHandle<evd::RecoDrawingOptions>  recoOpt;
    art::ServiceHandle<geo::Geometry>            geo;

    if (rawOpt->fDrawRawDataOrCalibWires < 1) return;
    if (recoOpt->fDraw2DEndPoints == 0)       return;

    geo::View_t gview = geo->TPC(rawOpt->fTPC).Plane(plane).View();

    for(size_t imod = 0; imod < recoOpt->fEndPoint2DLabels.size(); ++imod) {
      std::string const which = recoOpt->fEndPoint2DLabels[imod];
    
      art::PtrVector<recob::EndPoint2D> ep2d;
      this->GetEndPoint2D(evt, which, ep2d);
	
      for (size_t iep = 0; iep < ep2d.size(); ++iep) {

	// only worry about end points with the correct view
	if(ep2d[iep]->View() != gview) continue;

	///\todo - have to verify that we are in the right TPC, but to do that we
	// need to be sure that all EndPoint2D objects have filled the required information

	// draw cluster with unique marker  
	// Place this cluster's unique marker at the hit's location
	int color  = evd::kColor[ep2d[iep]->ID()%evd::kNCOLS];
	
	double x = ep2d[iep]->WireID().Wire;
	double y = ep2d[iep]->DriftTime();

	if(rawOpt->fAxisOrientation > 0){
	  x = ep2d[iep]->DriftTime();
	  y = ep2d[iep]->WireID().Wire;
	}

	TMarker& strt = view->AddMarker(x, y, color, 30, 2.0);
	strt.SetMarkerColor(color);
	
      } // loop on iep end points
    } // loop on imod folders

    return;
  }

  //......................................................................
  void RecoBaseDrawer::OpFlash2D(const art::Event& evt,
				 evdb::View2D*     view,
				 unsigned int      plane)
  {
    art::ServiceHandle<evd::RawDrawingOptions>   rawOpt;
    art::ServiceHandle<evd::RecoDrawingOptions>  recoOpt;
    art::ServiceHandle<geo::Geometry>            geo;

    if (rawOpt->fDrawRawDataOrCalibWires < 1) return;
    if (recoOpt->fDrawOpFlashes == 0)         return;

    for(size_t imod = 0; imod < recoOpt->fOpFlashLabels.size(); ++imod) {
      std::string const which = recoOpt->fOpFlashLabels[imod];

      art::PtrVector<recob::OpFlash> opflashes;
      this->GetOpFlashes(evt, which, opflashes);

      if(opflashes.size() < 1) continue;


      int NFlashes = opflashes.size();
      double TopCoord = 1000;

      // project each seed into this view
      for (size_t iof = 0; iof < opflashes.size(); ++iof) {
	std::vector<double> WireCenters = opflashes[iof]->WireCenters();
	std::vector<double> WireWidths = opflashes[iof]->WireWidths();
	
	LOG_VERBATIM("RecoBaseDrawer") << "Flash t: " 
				       << opflashes[iof]->Time() 
				       << "\t y,z : " 
				       << opflashes[iof]->YCenter() 
				       << ", "
				       << opflashes[iof]->ZCenter()
				       << " \t"
				       << opflashes[iof]->TotalPE();
		
        double LineTop = TopCoord * float(iof) / NFlashes;

        double x1 =  WireCenters.at(plane)+WireWidths.at(plane);
        double x2 =  WireCenters.at(plane)-WireWidths.at(plane);
        double y1 =  LineTop;
        double y2 =  LineTop;

	double s1x1 = x1;
	double s1x2 = x1;
	double s2x1 = x2;
        double s2x2 = x2;
        double s1y1 = 0;
	double s1y2 = LineTop;
        double s2y1 = 0;
        double s2y2 = LineTop;

        if(opflashes[iof]->OnBeamTime()==1)
          {
            s1y2 = TopCoord*1.1;
            s2y2 = TopCoord*1.1;
            y2 = TopCoord*1.1;
            y1 = TopCoord*1.1;
	  }

	if(rawOpt->fAxisOrientation > 0)
          {
            y1 =   WireCenters.at(plane)+WireWidths.at(plane);
            y2 =   WireCenters.at(plane)-WireWidths.at(plane);
            x1 =   LineTop;
	    x2 =   LineTop;

            s1y1 = x1;
            s1y2 = x1;
            s2y1 = x2;
            s2y2 = x2;
            s1x1 = 0;
            s1x2 = LineTop;
            s2x1 = 0;
            s2x2 = LineTop;
	    if(opflashes[iof]->OnBeamTime()==1)
              {
                s1x2 = TopCoord*1.1;
                s2x2 = TopCoord*1.1;
                x2 = TopCoord*1.1;
		x1 = TopCoord*1.1;
              }
          }
	
        TLine&   line = view->AddLine(x1, y1, x2, y2);


	if(opflashes[iof]->OnBeamTime()==1)
	  {
	    
	    line.SetLineColor(kRed);
	    TLine&   s1 = view->AddLine(s1x1, s1y1, s1x2, s1y2);
	    TLine&   s2 = view->AddLine(s2x1, s2y1, s2x2, s2y2);
	    s1.SetLineColor(kRed);
	    s2.SetLineColor(kRed);
	    s1.SetLineStyle(kDashed);
	    s2.SetLineStyle(kDashed);
	  }
	
	else
	  {
	    //      line.SetLineColor(kGray);
	  }
	line.SetLineStyle(kSolid);
	line.SetLineWidth(2.0);
	
      } // loop on opflashes
    } // loop on imod folders
    
    return;
  }
  
  
  //......................................................................
  void RecoBaseDrawer::Seed2D(const art::Event& evt,
			      evdb::View2D*     view,
			      unsigned int      plane) 
  {
    art::ServiceHandle<evd::RawDrawingOptions>   rawOpt;
    art::ServiceHandle<evd::RecoDrawingOptions>  recoOpt;
    art::ServiceHandle<geo::Geometry>            geo;

    if (rawOpt->fDrawRawDataOrCalibWires < 1) return;
    if (recoOpt->fDrawSeeds == 0)             return;

    for(size_t imod = 0; imod < recoOpt->fSeedLabels.size(); ++imod) {
      std::string const which = recoOpt->fSeedLabels[imod];
    
      art::PtrVector<recob::Seed> seeds;
      this->GetSeeds(evt, which, seeds);

      if(seeds.size() < 1) continue;
	
      // project each seed into this view
      for (size_t isd = 0; isd < seeds.size(); ++isd) {

	double SeedPoint[3];
	double SeedDir[3];
	double SeedPointErr[3];
	double SeedDirErr[3];
	double SeedEnd1[3];
	double SeedEnd2[3];

	seeds[isd]->GetPoint( SeedPoint,  SeedPointErr);
	seeds[isd]->GetDirection( SeedDir, SeedDirErr);
	
	SeedEnd1[0] = SeedPoint[0] + SeedDir[0];
	SeedEnd1[1] = SeedPoint[1] + SeedDir[1];
	SeedEnd1[2] = SeedPoint[2] + SeedDir[2];

	SeedEnd2[0] = SeedPoint[0] - SeedDir[0] ;
	SeedEnd2[1] = SeedPoint[1] - SeedDir[1] ;
	SeedEnd2[2] = SeedPoint[2] - SeedDir[2] ;
	
	// Draw seed on evd
	// int color  = kColor[seeds[isd]->ID()%kNCOLS];
	int color  = evd::kColor[0];
	unsigned int wirepoint = 0;
	unsigned int wireend1  = 0;
	unsigned int wireend2  = 0;
	try{
	  wirepoint = geo->NearestWire(SeedPoint, plane, rawOpt->fTPC, rawOpt->fCryostat);
	}
	catch(cet::exception &e){
	  wirepoint = atoi(e.explain_self().substr(e.explain_self().find("#")+1,5).c_str());
	}
	try{
	  wireend1  = geo->NearestWire(SeedEnd1,  plane, rawOpt->fTPC, rawOpt->fCryostat);
	}
	catch(cet::exception &e){
	  wireend1 = atoi(e.explain_self().substr(e.explain_self().find("#")+1,5).c_str());
	}
	try{
	  wireend2  = geo->NearestWire(SeedEnd2,  plane, rawOpt->fTPC, rawOpt->fCryostat);
	}
	catch(cet::exception &e){
	  wireend2 = atoi(e.explain_self().substr(e.explain_self().find("#")+1,5).c_str());
	}

	art::ServiceHandle<util::DetectorProperties> det;

	double x =  wirepoint;
	double y =  det->ConvertXToTicks(SeedPoint[0], plane, rawOpt->fTPC, rawOpt->fCryostat);
	double x1 = wireend1;
	double y1 = det->ConvertXToTicks(SeedEnd1[0],  plane, rawOpt->fTPC, rawOpt->fCryostat);
	double x2 = wireend2;
	double y2 = det->ConvertXToTicks(SeedEnd2[0],  plane, rawOpt->fTPC, rawOpt->fCryostat);

	if(rawOpt->fAxisOrientation > 0){
	  x  = det->ConvertXToTicks(SeedPoint[0], plane, rawOpt->fTPC, rawOpt->fCryostat);
	  y  = wirepoint;
	  x1 = det->ConvertXToTicks(SeedEnd1[0],  plane, rawOpt->fTPC, rawOpt->fCryostat);
	  y1 = wireend1;		
	  x2 = det->ConvertXToTicks(SeedEnd2[0],  plane, rawOpt->fTPC, rawOpt->fCryostat);
	  y2 = wireend2;		
	}

	TMarker& strt = view->AddMarker(x, y, color, 4, 1.5);
	TLine&   line = view->AddLine(x1, y1, x2, y2);
	strt.SetMarkerColor(color);
	line.SetLineColor(color);
	line.SetLineWidth(2.0);
	
      } // loop on seeds
    } // loop on imod folders

    return;
  }


  //......................................................................
  void RecoBaseDrawer::BezierTrack2D(const art::Event& evt,
                                     evdb::View2D*     view,
                                     unsigned int      plane)
  {
    art::ServiceHandle<evd::RawDrawingOptions>   rawOpt;
    art::ServiceHandle<evd::RecoDrawingOptions>  recoOpt;
    art::ServiceHandle<geo::Geometry>            geo;

    if (rawOpt->fDrawRawDataOrCalibWires < 1) return;
    if (recoOpt->fDrawBezierTracks == 0)      return;

    for(size_t imod = 0; imod < recoOpt->fBezierTrackLabels.size(); ++imod) {
      std::string const which = recoOpt->fBezierTrackLabels[imod];

      art::PtrVector<recob::Track> btrackbase;
      this->GetBezierTracks(evt, which, btrackbase);

      int N=100;

      // project each seed into this view
      for (size_t ibtb = 0; ibtb < btrackbase.size(); ++ibtb) {

	trkf::BezierTrack BTrack(*btrackbase.at(ibtb));

	std::vector<std::vector<double> > ProjPtUVWs(3);
	std::vector<std::vector<double> > ProjTimes(3);

	double projpt[3], ticks[3];
        int c=0, t=0;
        
        for(int i = 0; i != N; ++i){
	  try{
	    BTrack.GetProjectedPointUVWT(float(i)/N,projpt,ticks,c,t );
	    for(size_t n = 0; n != 3; ++n){
	      ProjPtUVWs[n].push_back(projpt[n]);
	      ProjTimes[n].push_back(ticks[n]);
	    }
	  }
	  catch(...){
	    continue;
	  }
	}

	TPolyLine& pl = view->AddPolyLine(ProjPtUVWs[plane].size(),kColor[ibtb%kNCOLS],1,0);
	
	for(size_t i = 0; i != ProjPtUVWs[plane].size(); ++i){
	  double x = ProjPtUVWs[plane][i];
	  double y = ProjTimes[plane][i];
	  if(rawOpt->fAxisOrientation > 0){
                y =  ProjPtUVWs[plane][i];
                x =  ProjTimes[plane][i];
	  }
	  pl.SetPoint(i,x,y);
	}
      }
    }
  }

  //......................................................................
  void RecoBaseDrawer::BezierTrack3D(const art::Event& evt,
                                     evdb::View3D*       view)
  {
    art::ServiceHandle<evd::RawDrawingOptions>   rawOpt;
    art::ServiceHandle<evd::RecoDrawingOptions>  recoOpt;
    art::ServiceHandle<geo::Geometry>            geo;

    if (rawOpt->fDrawRawDataOrCalibWires < 1) return;
    if (recoOpt->fDrawBezierTracks == 0)      return;

    for(size_t imod = 0; imod < recoOpt->fBezierTrackLabels.size(); ++imod) {
      std::string const which = recoOpt->fBezierTrackLabels[imod];

      art::PtrVector<recob::Track> btrackbase;
      art::PtrVector<recob::Vertex> btrackvertices;
      this->GetBezierTracks(evt, which, btrackbase);
      this->GetVertices(evt, which, btrackvertices);

      int N=100;



      // Draw bezier track lines
      for (size_t ibtb = 0; ibtb < btrackbase.size(); ++ibtb) {

	trkf::BezierTrack BTrack(*btrackbase.at(ibtb));
        TPolyLine3D& pl = view->AddPolyLine3D(N,kColor[ibtb%kNCOLS],2,0);
        double xyzpt[3] ;

        for(int i = 0; i != N; ++i){
	  BTrack.GetTrackPoint(float(i)/N,xyzpt );
	  double x = xyzpt[0];
	  double y = xyzpt[1];
	  double z = xyzpt[2];
	  
	  pl.SetPoint(i,x,y,z);
	  
	}
      }
      
      // Draw bezier track vertices

      TPolyMarker3D&  pmrk = view->AddPolyMarker3D(btrackvertices.size(), kYellow, 4, 1);
      for(size_t ivtx = 0; ivtx < btrackvertices.size(); ++ivtx){
	double xyz[3];
	btrackvertices.at(ivtx)->XYZ(xyz);
	pmrk.SetPoint(ivtx, xyz[0], xyz[1], xyz[2]);	  	  
      }
    }

  }

  //......................................................................
  void RecoBaseDrawer::Cluster2D(const art::Event& evt,
				 evdb::View2D*     view,
				 unsigned int      plane) 
  {
    art::ServiceHandle<evd::RawDrawingOptions>   rawOpt;
    art::ServiceHandle<evd::RecoDrawingOptions>  recoOpt;
    art::ServiceHandle<geo::Geometry>            geo;

    if (rawOpt->fDrawRawDataOrCalibWires < 1) return;
    if (recoOpt->fDrawClusters == 0)          return;

    geo::View_t gview = geo->TPC(rawOpt->fTPC).Plane(plane).View();

    // if user sets "DrawClusters" to 2, draw the clusters differently:
    bool drawAsMarkers = (recoOpt->fDrawClusters == 1 || 
                          recoOpt->fDrawClusters == 3);
    // draw connecting lines between cluster hits?
    bool drawConnectingLines = (recoOpt->fDrawClusters == 3);

    for(size_t imod = 0; imod < recoOpt->fClusterLabels.size(); ++imod) {
      std::string const which = recoOpt->fClusterLabels[imod];
    
      art::PtrVector<recob::Cluster> clust;
      this->GetClusters(evt, which, clust);

      if(clust.size() < 1) continue;

      art::FindMany<recob::Hit> fmh(clust, evt, which);

      for (size_t ic = 0; ic < clust.size(); ++ic) {

	// only worry about clusters with the correct view
	if(clust[ic]->View() != gview) continue;

 	std::vector<const recob::Hit*> hits;
	
 	hits = fmh.at(ic);

	// check for correct tpc, the view check done above
	// ensures we are in the correct plane
	if((*hits.begin())->WireID().TPC      != rawOpt->fTPC || 
	   (*hits.begin())->WireID().Cryostat != rawOpt->fCryostat) continue;

	if (drawAsMarkers) {
	  // draw cluster with unique marker  
	  // Place this cluster's unique marker at the hit's location
	  int color  = evd::kColor[clust[ic]->ID()%evd::kNCOLS];
	  this->Hit2D(hits, color, view, drawConnectingLines);
          if(recoOpt->fDrawClusters > 3) {
            // BB: draw the cluster ID 
            std::string s = std::to_string(clust[ic]->ID());
            char const* txt = s.c_str();
            double wire = clust[ic]->StartWire();
            double tick = 20 + clust[ic]->StartTick();
            TText& clID = view->AddText(wire, tick, txt);
            clID.SetTextColor(color);
          } // recoOpt->fDrawClusters > 3
	}
	else {

	  // default "outline" method:
	  std::vector<double> tpts, wpts;
      
	  this->GetClusterOutlines(hits, tpts, wpts, plane);
      
	  int lcolor = 9; // line color
	  int fcolor = 9; // fill color
	  int width  = 2; // line width
	  int style  = 1; // 1=solid line style
	  if (view != 0) {
	    TPolyLine& p1 = view->AddPolyLine(wpts.size(), 
					      lcolor,
					      width,
					      style);
	    TPolyLine& p2 = view->AddPolyLine(wpts.size(),
					      lcolor,
					      width,
					      style);
	    p1.SetOption("f");
	    p1.SetFillStyle(3003);
	    p1.SetFillColor(fcolor);
	    for (size_t i = 0; i < wpts.size(); ++i) {
	      if(rawOpt->fAxisOrientation < 1){
		p1.SetPoint(i, wpts[i], tpts[i]);
		p2.SetPoint(i, wpts[i], tpts[i]);
	      }
	      else{
		p1.SetPoint(i, tpts[i], wpts[i]);
		p2.SetPoint(i, tpts[i], wpts[i]);
	      }
	    } // loop on i points in ZX view
	  } // if we have a cluster in the ZX view
	}// end if outline mode

	// draw the direction cosine of the cluster as well as it's starting point
	// (average of the start and end angle -- by default they are the same value)
	this->Draw2DSlopeEndPoints(clust[ic]->StartWire(), 
				   clust[ic]->StartTick(), 
				   std::tan((clust[ic]->StartAngle() + clust[ic]->EndAngle())/2.),
				   evd::kColor[clust[ic]->ID()%evd::kNCOLS], view);

      } // loop on ic clusters
    } // loop on imod folders

    return;
  }

  //......................................................................
  void RecoBaseDrawer::Draw2DSlopeEndPoints(double        x,
					    double        y,
					    double        slope,
					    int           color,
					    evdb::View2D* view)
  {
    art::ServiceHandle<evd::RawDrawingOptions>   rawOpt;
    art::ServiceHandle<evd::RecoDrawingOptions>  recoOpt;
    
    if(recoOpt->fDraw2DSlopeEndPoints < 1) return;

    double x1 = x;
    double y1 = y;
    double slope1 = slope;

    if(rawOpt->fDrawRawDataOrCalibWires < 1) return;
    if(rawOpt->fAxisOrientation > 0){
      x1 = y;
      y1 = x;
      if(std::abs(slope) > 0.) slope1 = 1./slope;
      else slope1 = 1.e6;
    }		

    TMarker& strt = view->AddMarker(x1, y1, color, 29, 2.0);
    strt.SetMarkerColor(color); // stupid line to shut up compiler warning

    //    double stublen = 50.0 ;
    double stublen = 300.0;
    TLine& l = view->AddLine(x1, y1, x1+stublen, y1 + slope1*stublen);
    l.SetLineColor(color);
    l.SetLineWidth(2);

    return;
  }

  //......................................................................
  ///
  /// Make a set of points which outline a cluster
  ///
  /// @param c      : Reco base cluster to outline
  /// @param wpts   : wire values of the outlines
  /// @param tpts   : tdc values of the outlines
  /// @param plane  : plane number 
  ///
  void RecoBaseDrawer::GetClusterOutlines(std::vector<const recob::Hit*>& hits,
					  std::vector<double>&            wpts, 
					  std::vector<double>&      	  tpts, 
					  unsigned int              	  plane)
  {
    art::ServiceHandle<evd::RawDrawingOptions>   rawOpt;
        
    // Map wire numbers to highest and lowest in the plane
    std::map<unsigned int, double> wlo, whi;
    // On first pass, initialize
    for(size_t j = 0; j < hits.size(); ++j){

      // check that we are on the correct plane and TPC
      if(hits[j]->WireID().Plane    != plane        || 
	 hits[j]->WireID().TPC      != rawOpt->fTPC || 
	 hits[j]->WireID().Cryostat != rawOpt->fCryostat) continue;

      wlo[hits[j]->WireID().Wire] = hits[j]->PeakTime();
      whi[hits[j]->WireID().Wire] = hits[j]->PeakTime();
    }

    double t = 0.;

    // Finalize on second pass
    for (size_t j = 0; j < hits.size(); ++j) {
      t = hits[j]->PeakTime();

      if (t < wlo[hits[j]->WireID().Wire]) wlo[hits[j]->WireID().Wire] = t;
      if (t > whi[hits[j]->WireID().Wire]) whi[hits[j]->WireID().Wire] = t;
    }
    
    // Loop over wires and low times to make lines along bottom
    // edge. Work from upstream edge to downstream edge
    std::map<unsigned int, double>::iterator itr   (wlo.begin());
    std::map<unsigned int, double>::iterator itrEnd(wlo.end());
    for (; itr != itrEnd; ++itr) {
      unsigned int w = itr->first;
      t = itr->second;
      
      wpts.push_back(1.*w-0.1); tpts.push_back(t-0.1);
      wpts.push_back(1.*w+0.1); tpts.push_back(t-0.1);
    }
    
    // Loop over planes and high cells to make lines along top
    // edge. Work from downstream edge toward upstream edge
    std::map<unsigned int, double>::reverse_iterator ritr   (whi.rbegin());
    std::map<unsigned int, double>::reverse_iterator ritrEnd(whi.rend());
    for (; ritr != ritrEnd; ++ritr) {
      unsigned int w = ritr->first;
      t = ritr->second;
      
      wpts.push_back(1.*w+0.1); tpts.push_back(t+0.1);
      wpts.push_back(1.*w-0.1); tpts.push_back(t+0.1);
    }
    
    // Add link to starting point to close the box
    wpts.push_back(wpts[0]); tpts.push_back(tpts[0]);

    return;
  }

  //......................................................................
  void RecoBaseDrawer::DrawProng2D(std::vector<const recob::Hit*>&     hits,
				   evdb::View2D*                       view, 
				   unsigned int                        plane,
				   TVector3                     const& startPos,
				   TVector3                     const& startDir,	
				   int                                 id)   
  {
    art::ServiceHandle<util::DetectorProperties> detprop;
    art::ServiceHandle<evd::RawDrawingOptions>   rawOpt;
    art::ServiceHandle<geo::Geometry>            geo;
    art::ServiceHandle<util::LArProperties>      larp;
    unsigned int c = rawOpt->fCryostat;
    unsigned int t = rawOpt->fTPC;

    // first draw the hits
    this->Hit2D(hits, evd::kColor[id%evd::kNCOLS], view);

    // prepare to draw prongs
    double local[3] = {0.};
    double world[3] = {0.};
    geo->Cryostat(c).TPC(t).Plane(plane).LocalToWorld(local, world);
    world[1] = startPos.Y();
    world[2] = startPos.Z();

    // convert the starting position and direction from 3D to 2D coordinates
    double tick = detprop->ConvertXToTicks(startPos.X(), plane, t, c);
    double wire = 0.;
    try{
      wire = 1.*geo->NearestWire(world, plane, t, c);
    }
    catch(cet::exception &e){
      wire = 1.*atoi(e.explain_self().substr(e.explain_self().find("#")+1,5).c_str());
    }

   // thetawire is the angle measured CW from +z axis to wire
    double thetawire = geo->TPC(t).Plane(plane).Wire(0).ThetaZ();
    double wirePitch = geo->WirePitch(hits[0]->View());
    double driftvelocity = larp->DriftVelocity(); // cm/us
    double timetick = detprop->SamplingRate()*1e-3;  // time sample in us
    //rotate coord system CCW around x-axis by pi-thetawire
    //   new yprime direction is perpendicular to the wire direction
    //   in the same plane as the wires and in the direction of
    //   increasing wire number
    //use yprime-component of dir cos in rotated coord sys to get
    //   dTdW (number of time ticks per unit of wire pitch)
    double rotang = TMath::Pi() - thetawire;
    double yprime = std::cos(rotang)*startDir[1]
      +std::sin(rotang)*startDir[2];
    double dTdW = startDir[0]*wirePitch/driftvelocity/timetick/yprime;
    
    this->Draw2DSlopeEndPoints(wire, tick, dTdW, evd::kColor[id%evd::kNCOLS], view);	

    return;
  }

  //......................................................................
  void RecoBaseDrawer::Prong2D(const art::Event& evt,
			       evdb::View2D*     view,
			       unsigned int      plane) 
  {    
    art::ServiceHandle<evd::RawDrawingOptions>   rawOpt;
    art::ServiceHandle<evd::RecoDrawingOptions>  recoOpt;
    art::ServiceHandle<geo::Geometry>            geo;
    art::ServiceHandle<util::DetectorProperties> detprop;

    if(rawOpt->fDrawRawDataOrCalibWires < 1) return;

    geo::View_t gview = geo->TPC(rawOpt->fTPC).Plane(plane).View();

    // annoying for now, but have to have multiple copies of basically the
    // same code to draw prongs, showers and tracks so that we can use
    // the art::Assns to get the hits and clusters.

    unsigned int cstat = rawOpt->fCryostat;
    unsigned int tpc = rawOpt->fTPC;
    int tid = 0;

    if(recoOpt->fDrawTracks != 0){
      for(size_t imod = 0; imod < recoOpt->fTrackLabels.size(); ++imod){
	std::string const which = recoOpt->fTrackLabels[imod];

	art::View<recob::Track> track;
	this->GetTracks(evt, which, track);

	if(track.vals().size() < 1) continue;

	art::FindMany<recob::Hit>     fmh(track, evt, which);

	// loop over the prongs and get the clusters and hits associated with
	// them.  only keep those that are in this view
	for(size_t t = 0; t < track.vals().size(); ++t){

          if(recoOpt->fDrawTracks > 1) {
            // BB: draw the track ID at the end of the track
            double x = track.vals().at(t)->End()(0);
            double y = track.vals().at(t)->End()(1);
            double z = track.vals().at(t)->End()(2);
            double tick = 30 + detprop->ConvertXToTicks(x, plane, tpc, cstat);
            double wire = geo->WireCoordinate(y, z, plane, tpc, cstat);
            tid = track.vals().at(t)->ID();
            std::string s = std::to_string(tid);
            char const* txt = s.c_str();
	    TText& trkID = view->AddText(wire, tick, txt);
            trkID.SetTextColor(evd::kColor[tid%evd::kNCOLS]);
          }
	  
	  std::vector<const recob::Hit*> hits = fmh.at(t);

	  // only get the hits for the current view
	  std::vector<const recob::Hit*>::iterator itr = hits.begin(); 
	  while(itr < hits.end()){
	    if((*itr)->View() != gview) hits.erase(itr);
	    else itr++;
	  }

	  this->DrawProng2D(hits, view, plane, 
			    track.vals().at(t)->Vertex(), 
			    track.vals().at(t)->VertexDirection(), 
			    track.vals().at(t)->ID());
	}// end loop over prongs
      }// end loop over labels
    }// end draw tracks

    if(recoOpt->fDrawShowers != 0){
      for(size_t imod = 0; imod < recoOpt->fShowerLabels.size(); ++imod){
	std::string const which = recoOpt->fShowerLabels[imod];

	art::View<recob::Shower> shower;
	this->GetShowers(evt, which, shower);

	if(shower.vals().size() < 1) continue;

	art::FindMany<recob::Hit>     fmh(shower, evt, which);

	// loop over the prongs and get the clusters and hits associated with
	// them.  only keep those that are in this view
	for(size_t s = 0; s < shower.vals().size(); ++s){
	  
	  std::vector<const recob::Hit*> hits = fmh.at(s);
	  
	  // only get the hits for the current view
	  std::vector<const recob::Hit*>::iterator itr = hits.begin(); 
	  while(itr < hits.end()){
	    if((*itr)->View() != gview) hits.erase(itr);
	    else itr++;
	  }

	  // sort the hits
	  std::sort(hits.begin(), hits.end());

	  // use the first hit in the collection until
	  // we start filling shower vertex points
	  // get the center of the first hit for now
	  // to get the y and z position
	  double wireXYZ[3] = {0.};
	  unsigned int cstat = hits.front()->WireID().Cryostat;
	  unsigned int tpc   = hits.front()->WireID().TPC;
	  unsigned int plane = hits.front()->WireID().Plane;
	  unsigned int wire  = hits.front()->WireID().Wire;
	  geo->Cryostat(cstat).TPC(tpc).Plane(plane).Wire(wire).GetCenter(wireXYZ);
	  TVector3 startPos(detprop->ConvertTicksToX(hits.front()->PeakTime(), 
						     plane, rawOpt->fTPC, rawOpt->fCryostat),
			    wireXYZ[1], wireXYZ[2]);

	  this->DrawProng2D(hits, view, plane, 
			    startPos,
			    shower.vals().at(s)->Direction(),
			    shower.vals().at(s)->ID());
	}// end loop over prongs
      }// end loop over labels
    }// end draw showers

    return;
  }

  //......................................................................
  void RecoBaseDrawer::Vertex2D(const art::Event& evt,
				evdb::View2D*     view,
				unsigned int      plane) 
  {    
    art::ServiceHandle<evd::RawDrawingOptions>   rawOpt;
    art::ServiceHandle<evd::RecoDrawingOptions>  recoOpt;
    art::ServiceHandle<geo::Geometry>            geo;
    art::ServiceHandle<util::DetectorProperties> detprop;

    if(rawOpt->fDrawRawDataOrCalibWires < 1) return;

    if(recoOpt->fDrawVertices == 0) return;
    
    for(size_t imod = 0; imod < recoOpt->fVertexLabels.size(); ++imod) {
      std::string const which = recoOpt->fVertexLabels[imod];

      art::PtrVector<recob::Vertex> vertex;
      this->GetVertices(evt, which, vertex);

      if(vertex.size() < 1) continue;

      for(size_t v = 0; v < vertex.size(); ++v){
        // BB: draw polymarker at the vertex position in this plane
        double xyz[3];
        vertex[v]->XYZ(xyz);
        double wire = geo->WireCoordinate(xyz[1], xyz[2], plane, rawOpt->fTPC, rawOpt->fCryostat);
        double time = detprop->ConvertXToTicks(xyz[0], plane, rawOpt->fTPC, rawOpt->fCryostat);
        int color  = evd::kColor[vertex[v]->ID()%evd::kNCOLS];
        TMarker& strt = view->AddMarker(wire, time, color, 24, 1.0);
        strt.SetMarkerColor(color);
      } // v
    } // end loop over vertex module lables
    
    return;
  }

  //......................................................................
  void RecoBaseDrawer::Event2D(const art::Event& evt,
			       evdb::View2D*     view,
			       unsigned int      plane) 
  {    
    art::ServiceHandle<evd::RawDrawingOptions>   rawOpt;
    art::ServiceHandle<evd::RecoDrawingOptions>  recoOpt;
    art::ServiceHandle<geo::Geometry>            geo;

    if(rawOpt->fDrawRawDataOrCalibWires < 1) return;
   
    if(recoOpt->fDrawEvents != 0){
      
      geo::View_t gview = geo->TPC(rawOpt->fTPC).Plane(plane).View();

      for (unsigned int imod = 0; imod < recoOpt->fEventLabels.size(); ++imod) {
	std::string const which = recoOpt->fEventLabels[imod];

	art::PtrVector<recob::Event> event;
	this->GetEvents(evt, which, event);

	if(event.size() < 1) continue;

	art::FindMany<recob::Hit> fmh(event, evt, which);

	for(size_t e = 0; e < event.size(); ++e){

	  std::vector<const recob::Hit*> hits;

	  hits = fmh.at(e);
	  
	  // only get the hits for the current view
	  std::vector<const recob::Hit*>::iterator itr = hits.begin(); 
	  while(itr < hits.end()){
	    if((*itr)->View() != gview) hits.erase(itr);
	    else itr++;
	  }
	  
	  this->Hit2D(hits, evd::kColor[event[e]->ID()%evd::kNCOLS], view);
	} // end loop over events 
      } // end loop over event module lables
    } // end if we are drawing events

    return;
  }

  //......................................................................
  void RecoBaseDrawer::Seed3D(const art::Event&   evt,
			      evdb::View3D*       view)
  {
    art::ServiceHandle<evd::RawDrawingOptions>   rawOpt;
    art::ServiceHandle<evd::RecoDrawingOptions>  recoOpt;
    
    std::vector<std::string> labels;
    if(recoOpt->fDrawSeeds != 0)
      for(size_t imod = 0; imod < recoOpt->fSeedLabels.size(); ++imod)
	labels.push_back(recoOpt->fSeedLabels[imod]);
    
    for(size_t imod = 0; imod < labels.size(); ++imod) {
      std::string const which = labels[imod];

      art::PtrVector<recob::Seed> seeds;
      this->GetSeeds(evt, which, seeds);
      
      int color=0;
      
      if(seeds.size() < 1) continue;

      TPolyMarker3D&  pmrk = view->AddPolyMarker3D(seeds.size(), color, 4, 1);
      
      for(size_t iseed = 0; iseed != seeds.size(); ++iseed){
	double pt[3], pterr[3], dir[3], direrr[3];
	seeds.at(iseed)->GetPoint(pt, pterr);
	seeds.at(iseed)->GetDirection(dir, direrr);
	
	double end1[3], end2[3];
	for(int i = 0; i != 3; ++i){
	  end1[i] = pt[i] + dir[i] ;
	  end2[i] = pt[i] - dir[i] ;
	}
	
	TPolyLine3D& pline = view->AddPolyLine3D(2, color, 2, 0);
	
	pmrk.SetPoint(iseed, pt[0], pt[1], pt[2]);	  
	pline.SetPoint(0, end1[0], end1[1], end1[2]);
	pline.SetPoint(1, end2[0], end2[1], end2[2]);
      }// end loop over seeds
    }// end loop over module labels
    
    return;
  }

  //......................................................................
  void RecoBaseDrawer::SeedOrtho(const art::Event&   evt,
				 evd::OrthoProj_t    proj,
				 evdb::View2D*       view)
  {
    art::ServiceHandle<evd::RawDrawingOptions>   rawOpt;
    art::ServiceHandle<evd::RecoDrawingOptions>  recoOpt;
    
    std::vector<std::string> labels;
    if(recoOpt->fDrawSeeds != 0)
      for(size_t imod = 0; imod < recoOpt->fSeedLabels.size(); ++imod)
	labels.push_back(recoOpt->fSeedLabels[imod]);
    
    for(size_t imod = 0; imod < labels.size(); ++imod) {
      std::string const which = labels[imod];

      art::PtrVector<recob::Seed> seeds;
      this->GetSeeds(evt, which, seeds);
      
      int color=0;
      
      for(size_t iseed = 0; iseed != seeds.size(); ++iseed){
	double pt[3], pterr[3], dir[3], direrr[3];
	seeds.at(iseed)->GetPoint(pt, pterr);
	seeds.at(iseed)->GetDirection(dir, direrr);
	
	double end1[3], end2[3];
	for(int i = 0; i != 3; ++i){
	  end1[i] = pt[i] + dir[i] ;
	  end2[i] = pt[i] - dir[i] ;
	}
	
	if(proj == evd::kXY){
	  TMarker& strt = view->AddMarker(pt[1], pt[0], color, 4, 1.5);
	  TLine&   line = view->AddLine(end1[1], end1[0], end2[1], end2[0]);
	  strt.SetMarkerColor(evd::kColor[color]);
	  line.SetLineColor(evd::kColor[color]);
	  line.SetLineWidth(2.0);
	}
	else if(proj == evd::kXZ){
	  TMarker& strt = view->AddMarker(pt[2], pt[0], color, 4, 1.5);
	  TLine& line = view->AddLine(end1[2], end1[0], end2[2], end2[0]);
	  strt.SetMarkerColor(evd::kColor[color]);
	  line.SetLineColor(evd::kColor[color]);
	  line.SetLineWidth(2.0);
	}
	else{
	  if(proj != evd::kYZ)
	    throw cet::exception("RecoBaseDrawer:SeedOrtho") << "projection is not YZ as expected\n";

	  TMarker& strt = view->AddMarker(pt[2], pt[1], color, 4, 1.5);
	  TLine& line = view->AddLine(end1[2], end1[1], end2[2], end2[1]);
	  strt.SetMarkerColor(evd::kColor[color]);
	  line.SetLineColor(evd::kColor[color]);
	  line.SetLineWidth(2.0);
	}

      }
    }
  }

  //......................................................................
  void RecoBaseDrawer::SpacePoint3D(const art::Event& evt,
				    evdb::View3D*     view)
  {
    art::ServiceHandle<evd::RawDrawingOptions>   rawOpt;
    art::ServiceHandle<evd::RecoDrawingOptions>  recoOpt;

    if(rawOpt->fDrawRawDataOrCalibWires < 1) return;

    std::vector<std::string> labels;
    if(recoOpt->fDrawSpacePoints != 0){
      for(size_t imod = 0; imod < recoOpt->fSpacePointLabels.size(); ++imod) 
	labels.push_back(recoOpt->fSpacePointLabels[imod]);
    }

    for(size_t imod = 0; imod < labels.size(); ++imod) {
	std::string const which = labels[imod];
      
	std::vector<const recob::SpacePoint*> spts;
	this->GetSpacePoints(evt, which, spts);
	int color = imod;
	this->DrawSpacePoint3D(spts, color, view);
    }

    return;
  }

  //......................................................................
  void RecoBaseDrawer::Prong3D(const art::Event& evt,
			       evdb::View3D*     view)
  {
    art::ServiceHandle<evd::RawDrawingOptions>   rawOpt;
    art::ServiceHandle<evd::RecoDrawingOptions>  recoOpt;

    if(rawOpt->fDrawRawDataOrCalibWires < 1) return;

    // annoying for now, but have to have multiple copies of basically the
    // same code to draw prongs, showers and tracks so that we can use
    // the art::Assns to get the hits and clusters.

    // Tracks.

    if(recoOpt->fDrawTracks != 0){
      for(size_t imod = 0; imod < recoOpt->fTrackLabels.size(); ++imod) {
	std::string which = recoOpt->fTrackLabels[imod];
	art::View<recob::Track> track;
	this->GetTracks(evt, which, track);

	for(size_t t = 0; t < track.vals().size(); ++t) {
	  const recob::Track* ptrack = track.vals().at(t);
	  int color = ptrack->ID();

	  // Draw track using only embedded information.

	  DrawTrack3D(*ptrack, color, view);
	}
      }
    }

    // Showers.

    if(recoOpt->fDrawShowers != 0){
      for(size_t imod = 0; imod < recoOpt->fShowerLabels.size(); ++imod) {
	std::string which = recoOpt->fShowerLabels[imod];
	art::View<recob::Shower> shower;
	this->GetShowers(evt, which, shower);

	for(size_t s = 0; s < shower.vals().size(); ++s) {
	  const recob::Shower* pshower = shower.vals().at(s);
	  int color = pshower->ID();
	  DrawShower3D(*pshower, color, view);
	}
      }
    }

    return;
  }

  //......................................................................
  void RecoBaseDrawer::DrawSpacePoint3D(const std::vector<const recob::SpacePoint*>& spts,
					int                 color, 
					evdb::View3D*       view)
  {
    // Get services.

    art::ServiceHandle<evd::RecoDrawingOptions>  recoOpt;

    // Organize space points into separate collections according to the color 
    // we want them to be.
    // If option If option fColorSpacePointsByChisq is false, this means
    // having a single collection with color inherited from the prong
    // (specified by the argument color).

    std::map<int, std::vector<const recob::SpacePoint*> > spmap;   // Indexed by color.

    for(auto i : spts) {
      const recob::SpacePoint* pspt = i;
      if(pspt == 0) throw cet::exception("RecoBaseDrawer:DrawSpacePoint3D") << "space point is null\n";

      // By default use event display palette.

      int spcolor = evd::kColor[color%evd::kNCOLS];

      // For rainbow effect, choose root colors in range [51,100].
      // We are using 100=best (red), 51=worst (blue).

      if(recoOpt->fColorSpacePointsByChisq) {
	spcolor = 100 - 2.5 * pspt->Chisq();
	if(spcolor < 51)
	  spcolor = 51;
	if(spcolor > 100)
	  spcolor = 100;
      }
      spmap[spcolor].push_back(pspt);
    }

    // Loop over colors.
    // Note that larger (=better) space points are plotted on
    // top for optimal visibility.

    for(auto const icolor : spmap) {
      int spcolor = icolor.first;
      const std::vector<const recob::SpacePoint*>& psps = icolor.second;

      // Make and fill a polymarker.

      TPolyMarker3D& pm = view->AddPolyMarker3D(psps.size(), spcolor, 1, 3);

      for(size_t s = 0; s < psps.size(); ++s){
	const recob::SpacePoint& spt = *psps[s];

	const double *xyz = spt.XYZ();
	pm.SetPoint(s, xyz[0], xyz[1], xyz[2]);
      }
    }
    
    return;
  }

  //......................................................................
  void RecoBaseDrawer::DrawTrack3D(const recob::Track& track,
				   int                 color, 
				   evdb::View3D*       view)
  {
    // Get options.

    art::ServiceHandle<evd::RecoDrawingOptions> recoOpt;

    if(recoOpt->fDrawTrackSpacePoints) {

      // Use brute force to find the module label and index of this
      // track, so that we can find associated space points and draw
      // them.

      const art::Event *evt = evdb::EventHolder::Instance()->GetEvent();
      std::vector<art::Handle<std::vector<recob::Track> > > handles;
      evt->getManyByType(handles);
      for(auto ih : handles) {
	const art::Handle<std::vector<recob::Track> > handle = ih;
	
	if(handle.isValid()) {
	  const std::string& which = handle.provenance()->moduleLabel();
	  art::FindMany<recob::SpacePoint> fmsp(handle, *evt, which);

	  int n = handle->size();
	  for(int i=0; i<n; ++i) {
	    art::Ptr<recob::Track> p(handle, i);
	    if(&*p == &track) {
	      std::vector<const recob::SpacePoint*> spts = fmsp.at(i);
	      DrawSpacePoint3D(spts, color, view);
	    }
	  }
	}
      }
      
    }
    if(recoOpt->fDrawTrackTrajectoryPoints) {
    
      // Draw trajectory points as a polyline

      int np = track.NumberTrajectoryPoints();

      TPolyLine3D& pl = view->AddPolyLine3D(np, evd::kColor[color%evd::kNCOLS], 2, 0);
      
      for(int p = 0; p < np; ++p) {
	const TVector3& pos = track.LocationAtPoint(p);
	pl.SetPoint(p, pos.X(), pos.Y(), pos.Z());
      } // p
    } // recoOpt->fDrawTrackTrajectoryPoints

    return;
  }

  //......................................................................
  void RecoBaseDrawer::DrawShower3D(const recob::Shower& shower,
				    int                 color, 
				    evdb::View3D*       view)
  {
    // Use brute force to find the module label and index of this
    // shower, so that we can find associated space points and draw
    // them.

    const art::Event *evt = evdb::EventHolder::Instance()->GetEvent();
    std::vector<art::Handle<std::vector<recob::Shower> > > handles;
    evt->getManyByType(handles);
    for(auto ih : handles) {
      const art::Handle<std::vector<recob::Shower> > handle = ih;

      if(handle.isValid()) {

	const std::string& which = handle.provenance()->moduleLabel();
	art::FindMany<recob::SpacePoint> fmsp(handle, *evt, which);

	int n = handle->size();
	for(int i=0; i<n; ++i) {
	  art::Ptr<recob::Shower> p(handle, i);
	  if(&*p == &shower) {
	    std::vector<const recob::SpacePoint*> spts = fmsp.at(i);
	    DrawSpacePoint3D(spts, color, view);
	  }
	}
      }
    }
          
    return;
  }

  //......................................................................
  void RecoBaseDrawer::Vertex3D(const art::Event& evt,
				evdb::View3D*     view)
  {
    art::ServiceHandle<evd::RawDrawingOptions>   rawOpt;
    art::ServiceHandle<evd::RecoDrawingOptions>  recoOpt;
   
    if(rawOpt->fDrawRawDataOrCalibWires < 1) return;

    if(recoOpt->fDrawVertices != 0){

      for (size_t imod = 0; imod < recoOpt->fVertexLabels.size(); ++imod) {
	std::string const which = recoOpt->fVertexLabels[imod];
	
	art::PtrVector<recob::Vertex> vertex;
	this->GetVertices(evt, which, vertex);

	art::FindManyP<recob::Track>  fmt(vertex, evt, which);
	art::FindManyP<recob::Shower> fms(vertex, evt, which);
	
	for(size_t v = 0; v < vertex.size(); ++v){

	  std::vector< art::Ptr<recob::Track> >  tracks  = fmt.at(v);
	  std::vector< art::Ptr<recob::Shower> > showers = fms.at(v);
	  
	  // grab the Prongs from the vertex and draw those
	  for(size_t t = 0; t < tracks.size(); ++t)
	    this->DrawTrack3D(*(tracks[t]), vertex[v]->ID(), view);
	  
	  for(size_t s = 0; s < showers.size(); ++s)
	    this->DrawShower3D(*(showers[s]), vertex[v]->ID(), view);
	  
	  double xyz[3] = {0.};
	  vertex[v]->XYZ(xyz);
	  TPolyMarker3D& pm = view->AddPolyMarker3D(1, evd::kColor[vertex[v]->ID()%evd::kNCOLS], 29, 6);
	  pm.SetPoint(0, xyz[0], xyz[1], xyz[2]);
	  
	} // end loop over vertices to draw from this label
      } // end loop over vertex module lables
    } // end if we are drawing vertices

    return;
  }

  //......................................................................
  void RecoBaseDrawer::Event3D(const art::Event& evt,
			       evdb::View3D*     view) 
  {    
    art::ServiceHandle<evd::RawDrawingOptions>   rawOpt;
    art::ServiceHandle<evd::RecoDrawingOptions>  recoOpt;
    
    if(rawOpt->fDrawRawDataOrCalibWires < 1) return;
    if(recoOpt->fDrawEvents != 0){

      for (size_t imod = 0; imod < recoOpt->fEventLabels.size(); ++imod) {
	std::string const which = recoOpt->fEventLabels[imod];

	art::PtrVector<recob::Event> event;
	this->GetEvents(evt, which, event);

	if(event.size() < 1) continue;

	art::FindManyP<recob::Vertex> fmvp(event, evt, which);
	art::FindMany<recob::Vertex>  fmv(event, evt, which);

	for(size_t e = 0; e < event.size(); ++e){

	  // grab the vertices for this event
	  std::vector< art::Ptr<recob::Vertex> > vertex = fmvp.at(e);

	  if(vertex.size() < 1) continue;

	  art::FindManyP<recob::Track>  fmt(vertex, evt, recoOpt->fVertexLabels[0]);
	  art::FindManyP<recob::Shower> fms(vertex, evt, recoOpt->fVertexLabels[0]);

	  for(size_t v = 0; v < vertex.size(); ++v){
	    
	    /// \todo need a better way to grab the vertex module labels, 
	    // right now assume there is only 1 in the list
	    std::vector< art::Ptr<recob::Track>  > tracks  = fmt.at(v);
	    std::vector< art::Ptr<recob::Shower> > showers = fms.at(v);

	    // grab the Prongs from the vertex and draw those
	    for(size_t t = 0; t < tracks.size(); ++t)
	      this->DrawTrack3D(*(tracks[t]), event[e]->ID(), view);

	    for(size_t s = 0; s < showers.size(); ++s)
	      this->DrawShower3D(*(showers[s]), event[e]->ID(), view);

	  } // end loop over vertices from this event

	  double xyz[3] = {0.};
	  std::vector<const recob::Vertex*> vts = fmv.at(e);

	  event[e]->PrimaryVertex(vts)->XYZ(xyz);
	  TPolyMarker3D& pm = view->AddPolyMarker3D(1, evd::kColor[event[e]->ID()%evd::kNCOLS], 29, 6);
	  pm.SetPoint(0, xyz[0], xyz[1], xyz[2]);

	} // end loop over events 
      } // end loop over event module lables
    } // end if we are drawing events

    return;
  }

  //......................................................................
  void RecoBaseDrawer::SpacePointOrtho(const art::Event& evt,
				       evd::OrthoProj_t  proj,
				       double            msize,
				       evdb::View2D*     view)
  {
    art::ServiceHandle<evd::RawDrawingOptions>   rawOpt;
    art::ServiceHandle<evd::RecoDrawingOptions>  recoOpt;

    if(rawOpt->fDrawRawDataOrCalibWires < 1) return;

    std::vector<std::string> labels;
    if(recoOpt->fDrawSpacePoints != 0){
      for(size_t imod = 0; imod < recoOpt->fSpacePointLabels.size(); ++imod) 
	labels.push_back(recoOpt->fSpacePointLabels[imod]);
    }

    for(size_t imod = 0; imod < labels.size(); ++imod) {
	std::string const which = labels[imod];
      
	std::vector<const recob::SpacePoint*> spts;
	this->GetSpacePoints(evt, which, spts);
	int color = imod;
	this->DrawSpacePointOrtho(spts, color, proj, msize, view);
    }
      
    return;
  }

  //......................................................................
  void RecoBaseDrawer::ProngOrtho(const art::Event& evt,
				  evd::OrthoProj_t  proj,
				  double            msize,
				  evdb::View2D*     view)
  {
    art::ServiceHandle<evd::RawDrawingOptions>   rawOpt;
    art::ServiceHandle<evd::RecoDrawingOptions>  recoOpt;

    if(rawOpt->fDrawRawDataOrCalibWires < 1) return;

    // annoying for now, but have to have multiple copies of basically the
    // same code to draw prongs, showers and tracks so that we can use
    // the art::Assns to get the hits and clusters.

    // Tracks.

    if(recoOpt->fDrawTracks != 0){
      for(size_t imod = 0; imod < recoOpt->fTrackLabels.size(); ++imod) {
	std::string which = recoOpt->fTrackLabels[imod];
	art::View<recob::Track> track;
	this->GetTracks(evt, which, track);

	for(size_t t = 0; t < track.vals().size(); ++t) {
	  const recob::Track* ptrack = track.vals().at(t);
	  int color = ptrack->ID();

	  // Draw track using only embedded information.

	  DrawTrackOrtho(*ptrack, color, proj, msize, view);
	}
      }
    }

    // Showers.

    if(recoOpt->fDrawShowers != 0){
      for(size_t imod = 0; imod < recoOpt->fShowerLabels.size(); ++imod) {
	std::string which = recoOpt->fShowerLabels[imod];
	art::View<recob::Shower> shower;
	this->GetShowers(evt, which, shower);

	for(size_t s = 0; s < shower.vals().size(); ++s) {
	  const recob::Shower* pshower = shower.vals().at(s);
	  int color = pshower->ID();
	  DrawShowerOrtho(*pshower, color, proj, msize, view);
	}
      }
    }


    return;
  }

  //......................................................................
  void RecoBaseDrawer::DrawSpacePointOrtho(const std::vector<const recob::SpacePoint*>& spts,
					   int                 color, 
					   evd::OrthoProj_t    proj,
					   double              msize,
					   evdb::View2D*       view)
  {
    // Get services.

    art::ServiceHandle<evd::RecoDrawingOptions>  recoOpt;

    // Organize space points into separate collections according to the color 
    // we want them to be.
    // If option If option fColorSpacePointsByChisq is false, this means
    // having a single collection with color inherited from the prong
    // (specified by the argument color).

    std::map<int, std::vector<const recob::SpacePoint*> > spmap;   // Indexed by color.

    for(auto i : spts){
      const recob::SpacePoint* pspt = i;
      if(pspt == 0) throw cet::exception("RecoBaseDrawer:DrawSpacePointOrtho") << "spacepoint is null\n";

      // By default use event display palette.

      int spcolor = evd::kColor[color%evd::kNCOLS];

      // For rainbow effect, choose root colors in range [51,100].
      // We are using 100=best (red), 51=worst (blue).

      if(recoOpt->fColorSpacePointsByChisq) {
	spcolor = 100 - 2.5 * pspt->Chisq();
	if(spcolor < 51)
	  spcolor = 51;
	if(spcolor > 100)
	  spcolor = 100;
      }
      spmap[spcolor].push_back(pspt);
    }

    // Loop over colors.
    // Note that larger (=better) space points are plotted on
    // top for optimal visibility.

    for(auto icolor : spmap) {
      int spcolor = icolor.first;
      const std::vector<const recob::SpacePoint*>& psps = icolor.second;

      // Make and fill a polymarker.

      TPolyMarker& pm = view->AddPolyMarker(psps.size(), spcolor,
					    kFullCircle, msize);
      for(size_t s = 0; s < psps.size(); ++s){
	const recob::SpacePoint& spt = *psps[s];
	const double *xyz = spt.XYZ();
	switch (proj) {
	  case evd::kXY:
	    pm.SetPoint(s, xyz[0], xyz[1]);
	    break;
	  case evd::kXZ:
	    pm.SetPoint(s, xyz[2], xyz[0]);
	    break;
	  case evd::kYZ:
	    pm.SetPoint(s, xyz[2], xyz[1]);
	    break;
	  default:
	    throw cet::exception("RecoBaseDrawer") << __func__
	      << ": unknown projection #" << ((int) proj) << "\n";
	} // switch
      }
    }
    
    return;
  }

  //......................................................................
  void RecoBaseDrawer::DrawTrackOrtho(const recob::Track& track,
				      int                 color, 
				      evd::OrthoProj_t    proj,
				      double              msize,
				      evdb::View2D*       view)
  {
    // Get options.

    art::ServiceHandle<evd::RecoDrawingOptions> recoOpt;

    if(recoOpt->fDrawTrackSpacePoints) {

      // Use brute force to find the module label and index of this
      // track, so that we can find associated space points and draw
      // them.

      const art::Event *evt = evdb::EventHolder::Instance()->GetEvent();
      std::vector<art::Handle<std::vector<recob::Track> > > handles;
      evt->getManyByType(handles);
      for(auto ih : handles) {
	const art::Handle<std::vector<recob::Track> > handle = ih;
	if(handle.isValid()) {
	  const std::string& which = handle.provenance()->moduleLabel();
	  art::FindMany<recob::SpacePoint> fmsp(handle, *evt, which);

	  int n = handle->size();
	  for(int i=0; i<n; ++i) {
	    art::Ptr<recob::Track> p(handle, i);
	    if(&*p == &track) {
	      std::vector<const recob::SpacePoint*> spts = fmsp.at(i);
	      DrawSpacePointOrtho(spts, color, proj, msize, view);
	    }
	  }
	}
      }
      
    }
    if(recoOpt->fDrawTrackTrajectoryPoints) {

      // Draw trajectory points.

      int np = track.NumberTrajectoryPoints();

      // Make and fill a polymarker.

      TPolyMarker& pm = view->AddPolyMarker(np, evd::kColor[color%evd::kNCOLS], kFullCircle, msize);
      TPolyLine& pl = view->AddPolyLine(np, evd::kColor[color%evd::kNCOLS], 2, 0);
      for(int p = 0; p < np; ++p){
	const TVector3& pos = track.LocationAtPoint(p);
	switch (proj) {
	  case evd::kXY:
	    pm.SetPoint(p, pos.X(), pos.Y());
	    pl.SetPoint(p, pos.X(), pos.Y());
	    break;
	  case evd::kXZ:
	    pm.SetPoint(p, pos.Z(), pos.X());
	    pl.SetPoint(p, pos.Z(), pos.X());
	    break;
	  case evd::kYZ:
	    pm.SetPoint(p, pos.Z(), pos.Y());
	    pl.SetPoint(p, pos.Z(), pos.Y());
	    break;
	  default:
	    throw cet::exception("RecoBaseDrawer") << __func__
	      << ": unknown projection #" << ((int) proj) << "\n";
	} // switch
      } // p
      // BB: draw the track ID at the end of the track
      if(recoOpt->fDrawTracks > 1) {
        int tid = track.ID();
        std::string s = std::to_string(tid);
        char const* txt = s.c_str();
        double x = track.End()(0);
        double y = track.End()(1);
        double z = track.End()(2);
        if(proj == evd::kXY) {
	  TText& trkID = view->AddText(x, y, txt);
          trkID.SetTextColor(evd::kColor[tid]);
          trkID.SetTextSize(0.03);
        } else if(proj == evd::kXZ) {
	  TText& trkID = view->AddText(z, x, txt);
          trkID.SetTextColor(evd::kColor[tid]);
          trkID.SetTextSize(0.03);
        } else if(proj == evd::kYZ) {
	  TText& trkID = view->AddText(z, y, txt);
          trkID.SetTextColor(evd::kColor[tid]);
          trkID.SetTextSize(0.03);
        } // proj
      } // recoOpt->fDrawTracks > 1
    }

    return;
  }

  //......................................................................
  void RecoBaseDrawer::DrawShowerOrtho(const recob::Shower& shower,
				       int                 color, 
				       evd::OrthoProj_t    proj,
				       double              msize,
				       evdb::View2D*       view)
  {
    // Use brute force to find the module label and index of this
    // shower, so that we can find associated space points and draw
    // them.

    const art::Event *evt = evdb::EventHolder::Instance()->GetEvent();
    std::vector<art::Handle<std::vector<recob::Shower> > > handles;
    evt->getManyByType(handles);
    for(auto ih : handles) {
      const art::Handle<std::vector<recob::Shower> > handle = ih;
      if(handle.isValid()) {
	const std::string& which = handle.provenance()->moduleLabel();
	art::FindMany<recob::SpacePoint> fmsp(handle, *evt, which);

	int n = handle->size();
	for(int i=0; i<n; ++i) {
	  art::Ptr<recob::Shower> p(handle, i);
	  if(&*p == &shower) {
	    std::vector<const recob::SpacePoint*> spts = fmsp.at(i);
	    DrawSpacePointOrtho(spts, color, proj, msize, view);
	  }
	}
      }
    }
      
    return;
  }

  //......................................................................
  int RecoBaseDrawer::GetWires(const art::Event&            evt,
			       const std::string&           which,
			       art::PtrVector<recob::Wire>& wires) 
  {
    wires.clear();

    art::Handle< std::vector<recob::Wire> > wcol;
    art::PtrVector<recob::Wire> temp;

    try{
      evt.getByLabel(which, wcol);
      
      for(unsigned int i = 0; i < wcol->size(); ++i){
	art::Ptr<recob::Wire> w(wcol, i);
	temp.push_back(w);
      }
      temp.swap(wires);
    }
    catch(cet::exception& e){
      writeErrMsg("GetWires", e);
    }

    return wires.size();
  }

  //......................................................................
  int RecoBaseDrawer::GetHits(const art::Event&               evt,
			      const std::string&              which,
			      std::vector<const recob::Hit*>& hits,
			      unsigned int                    plane) 
  {
    art::ServiceHandle<evd::RawDrawingOptions>   rawOpt;

    hits.clear();

    std::vector<const recob::Hit*> temp;

    try{
      evt.getView(which, temp);
      for(size_t t = 0; t < temp.size(); ++t){

	if( temp[t]->WireID().Plane    == plane        && 
	    temp[t]->WireID().TPC      == rawOpt->fTPC && 
	    temp[t]->WireID().Cryostat == rawOpt->fCryostat) hits.push_back(temp[t]);
      }
    }
    catch(cet::exception& e){
      writeErrMsg("GetHits", e);
    }

    return hits.size();
  }

  //......................................................................
  int RecoBaseDrawer::GetClusters(const art::Event&               evt, 
				  const std::string&              which,
				  art::PtrVector<recob::Cluster>& clust)
  {
    clust.clear();
    art::PtrVector<recob::Cluster> temp;

    art::Handle< std::vector<recob::Cluster> > clcol;

    try{
      evt.getByLabel(which, clcol);
      temp.reserve(clcol->size());
      for(unsigned int i = 0; i < clcol->size(); ++i){
        temp.push_back(art::Ptr<recob::Cluster>(clcol, i));
      }
      temp.swap(clust);
    }
    catch(cet::exception& e){
      writeErrMsg("GetClusters", e);
    }

    return clust.size();
  }

  //......................................................................
  int RecoBaseDrawer::GetEndPoint2D(const art::Event&                  evt, 
                                    const std::string&                 which,
                                    art::PtrVector<recob::EndPoint2D>& ep2d)
  {
    ep2d.clear();
    art::PtrVector<recob::EndPoint2D> temp;
   
    art::Handle< std::vector<recob::EndPoint2D> > epcol;
    
    try{
      evt.getByLabel(which, epcol);
      for(unsigned int i = 0; i < epcol->size(); ++i){
        art::Ptr<recob::EndPoint2D> ep(epcol, i);
        temp.push_back(ep);
      }
      temp.swap(ep2d);
    }
    catch(cet::exception& e){
      writeErrMsg("GetEndPoint2D", e);
    }
    
    return ep2d.size();
  }

  //......................................................................

  int RecoBaseDrawer::GetOpFlashes(const art::Event&                 evt,
				   const std::string&                 which,
				   art::PtrVector<recob::OpFlash>&    opflashes)
  {
    opflashes.clear();
    art::PtrVector<recob::OpFlash> temp;

    art::Handle< std::vector<recob::OpFlash> > opflashcol;

    try{
      evt.getByLabel(which, opflashcol);
      for(unsigned int i = 0; i < opflashcol->size(); ++i){
	art::Ptr<recob::OpFlash> opf(opflashcol, i);
	temp.push_back(opf);
      }
      temp.swap(opflashes);
    }
    catch(cet::exception& e){
      writeErrMsg("GetOpFlashes", e);
    }

    return opflashes.size();
  }

  //......................................................................

  int RecoBaseDrawer::GetSeeds(const art::Event&                  evt, 
			       const std::string&                 which,
			       art::PtrVector<recob::Seed>&       seeds)
  {
    seeds.clear();
    art::PtrVector<recob::Seed> temp;

    art::Handle< std::vector<recob::Seed> > seedcol;

    try{
      evt.getByLabel(which, seedcol);
      for(unsigned int i = 0; i < seedcol->size(); ++i){
	art::Ptr<recob::Seed> sd(seedcol, i);
	temp.push_back(sd);
      }
      temp.swap(seeds);
    }
    catch(cet::exception& e){
      writeErrMsg("GetSeeds", e);
    }

    return seeds.size();
  }


  //......................................................................
  int RecoBaseDrawer::GetBezierTracks(const art::Event&                  evt,
                                      const std::string&                 which,
                                      art::PtrVector<recob::Track>&  btbs)
  {
    btbs.clear();
    art::PtrVector<recob::Track> temp;

    art::Handle< std::vector<recob::Track> > btbcol;
    
    try{
      evt.getByLabel(which, btbcol);
      for(unsigned int i = 0; i < btbcol->size(); ++i){
	art::Ptr<recob::Track> btb(btbcol, i);
        temp.push_back(btb);
      }
      temp.swap(btbs);
    }
    catch(cet::exception& e){
      writeErrMsg("GetBezierTracks", e);
    }
    
    return btbs.size();
  }




  //......................................................................
  int RecoBaseDrawer::GetSpacePoints(const art::Event&               evt,
				     const std::string&              which,
				     std::vector<const recob::SpacePoint*>& spts)
  {
    spts.clear();
    std::vector<const recob::SpacePoint*> temp;
   
    art::Handle< std::vector<recob::SpacePoint> > spcol;
    
    try{
      evt.getByLabel(which, spcol);
      temp.reserve(spcol->size());
      for(unsigned int i = 0; i < spcol->size(); ++i){
        art::Ptr<recob::SpacePoint> spt(spcol, i);
        temp.push_back(&*spt);
      }
      temp.swap(spts);
    }
    catch(cet::exception& e){
      writeErrMsg("GetSpacePoints", e);
    }
    
    return spts.size();
  }

  //......................................................................
  int RecoBaseDrawer::GetTracks(const art::Event&        evt, 
				const std::string&       which,
				art::View<recob::Track>& track)
  {
    try{
      evt.getView(which,track);
    }
    catch(cet::exception& e){
      writeErrMsg("GetTracks", e);
    }

    return track.vals().size();
  }

  //......................................................................
  int RecoBaseDrawer::GetShowers(const art::Event&        evt, 
				const std::string&        which,
				art::View<recob::Shower>& shower)
  {
    try{
      evt.getView(which,shower);
    }
    catch(cet::exception& e){
      writeErrMsg("GetShowers", e);
    }

    return shower.vals().size();
  }

  //......................................................................
  int RecoBaseDrawer::GetVertices(const art::Event&              evt, 
				  const std::string&             which,
				  art::PtrVector<recob::Vertex>& vertex)
  {
    vertex.clear();
    art::PtrVector<recob::Vertex> temp;

    art::Handle< std::vector<recob::Vertex> > vcol;

    try{
      evt.getByLabel(which, vcol);
      for(size_t i = 0; i < vcol->size(); ++i){
	art::Ptr<recob::Vertex> v(vcol, i);
	temp.push_back(v);
      }
      temp.swap(vertex);
    }
    catch(cet::exception& e){
      writeErrMsg("GetVertices", e);
    }

    return vertex.size();
  }

  //......................................................................
  int RecoBaseDrawer::GetEvents(const art::Event&             evt, 
				const std::string&            which,
				art::PtrVector<recob::Event>& event)
  {
    event.clear();
    art::PtrVector<recob::Event> temp;

    art::Handle< std::vector<recob::Event> > ecol;

    try{
      evt.getByLabel(which, ecol);
      for(size_t i = 0; i < ecol->size(); ++i){
	art::Ptr<recob::Event> e(ecol, i);
	temp.push_back(e);
      }
      temp.swap(event);
    }
    catch(cet::exception& e){
      writeErrMsg("GetEvents", e);
    }

    return event.size();
  }

  //......................................................................
  void RecoBaseDrawer::FillTQHisto(const art::Event&    evt,
				   unsigned int         plane,
				   unsigned int         wire,
				   TH1F*                histo,
				   std::vector<double>& hstart,
				   std::vector<double>& hend,
				   std::vector<double>& hitamplitudes,
				   std::vector<double>& hpeaktimes)
  {
    art::ServiceHandle<evd::RawDrawingOptions>   rawOpt;
    art::ServiceHandle<evd::RecoDrawingOptions>  recoOpt;
    art::ServiceHandle<geo::Geometry>            geo;

    // Check if we're supposed to draw raw hits at all
    if(rawOpt->fDrawRawDataOrCalibWires==0) return;

    for (size_t imod = 0; imod < recoOpt->fWireLabels.size(); ++imod) {
      std::string const which = recoOpt->fWireLabels[imod];

      art::PtrVector<recob::Wire> wires;
      this->GetWires(evt, which, wires);

      for (size_t i = 0; i < wires.size(); ++i) {

	std::vector<geo::WireID> wireids = geo->ChannelToWire(wires[i]->Channel());

	bool goodWID = false;
	for( auto const& wid : wireids ){
	  // check for correct plane, wire and tpc
	  if(wid.Plane    == plane        &&
	     wid.Wire     == wire         &&
	     wid.TPC      == rawOpt->fTPC && 
	     wid.Cryostat == rawOpt->fCryostat) goodWID = true;
	}
	if(!goodWID) continue;

        std::vector<float> wirSig = wires[i]->Signal();
        for(unsigned int ii = 0; ii < wirSig.size(); ++ii) 
          histo->Fill(1.*ii, wirSig[ii]);
      }//end loop over wires
    }//end loop over wire modules

    for (size_t imod = 0; imod < recoOpt->fHitLabels.size(); ++imod) {
      std::string const which = recoOpt->fHitLabels[imod];

      std::vector<const recob::Hit*> hits;
      this->GetHits(evt, which, hits, plane);

      for (size_t i = 0; i < hits.size(); ++i){

	// check for correct wire, plane, cryostat and tpc were checked in GetHits
	if(hits[i]->WireID().Wire != wire) continue;

	hstart.push_back(hits[i]->PeakTimeMinusRMS());
	hend.push_back(hits[i]->PeakTimePlusRMS());
	hitamplitudes.push_back(hits[i]->PeakAmplitude());
	hpeaktimes.push_back(hits[i]->PeakTime());
	
      }//end loop over reco hits
    }//end loop over HitFinding modules

    return;
  }

  //......................................................................
  void RecoBaseDrawer::FillQHisto(const art::Event& evt,
				  unsigned int      plane,
				  TH1F*             histo)
  {
    art::ServiceHandle<evd::RawDrawingOptions>   rawOpt;
    art::ServiceHandle<evd::RecoDrawingOptions>  recoOpt;
    art::ServiceHandle<geo::Geometry>            geo;
 
    // Check if we're supposed to draw raw hits at all
    if(rawOpt->fDrawRawDataOrCalibWires==0) return;
  
    for (size_t imod = 0; imod < recoOpt->fWireLabels.size(); ++imod) {
      std::string const which = recoOpt->fWireLabels[imod];

      art::PtrVector<recob::Wire> wires;
      this->GetWires(evt, which, wires);

      for (unsigned int i=0; i<wires.size(); ++i) {

	std::vector<geo::WireID> wireids = geo->ChannelToWire(wires[i]->Channel());

	bool goodWID = false;
	for( auto const& wid : wireids ){
	  // check for correct plane, wire and tpc
	  if(wid.Plane    == plane        &&
	     wid.TPC      == rawOpt->fTPC && 
	     wid.Cryostat == rawOpt->fCryostat) goodWID = true;
	}
	if(!goodWID) continue;
	std::vector<float> wirSig = wires[i]->Signal();
        for(unsigned int ii = 0; ii < wirSig.size(); ++ii) 
          histo->Fill(wirSig[ii]);
/*
	for(size_t s = 0; s < wires[i]->NSignal(); ++s)
	  histo->Fill(wires[i]->Signal()[s]);
*/
      }//end loop over raw hits
    }//end loop over Wire modules

    return;
  }

}// namespace
////////////////////////////////////////////////////////////////////////
