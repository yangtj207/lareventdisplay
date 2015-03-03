/// \file    RawDataDrawer.cxx
/// \brief   Class to aid in the rendering of RawData objects
/// \author  messier@indiana.edu
/// \version $Id: RawDataDrawer.cxx,v 1.2 2010/11/11 18:11:22 p-novaart Exp $
#include <cmath>
#include <stdint.h>

#include "TH1F.h"
#include "TPolyLine3D.h"
#include "TBox.h"

#include "EventDisplay/RawDataDrawer.h"
#include "EventDisplayBase/View2D.h"
#include "EventDisplayBase/EventHolder.h"
#include "EventDisplay/ColorDrawingOptions.h"
#include "EventDisplay/RawDrawingOptions.h"
#include "RawData/raw.h"
#include "RawData/RawDigit.h"
#include "Filters/ChannelFilter.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Event.h"
#include "art/Persistency/Common/Ptr.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace {
  // Utility function to make uniform error messages.
  void writeErrMsg(const char* fcn,
		   cet::exception const& e)
  {
    mf::LogWarning("RawDataDrawer") << "RawDataDrawer::" << fcn
				    << " failed with message:\n"
				    << e;
  }
}


namespace evd {

  //......................................................................
  RawDataDrawer::RawDataDrawer() :
    fTicks(2048)
  { 
    art::ServiceHandle<geo::Geometry> geo;

    art::ServiceHandle<evd::RawDrawingOptions> rawopt;
    
    fTicks = rawopt->fTicks;

    // set the list of bad channels in this detector
    filter::ChannelFilter cf; 
    unsigned int nplanes=geo->Nplanes();
    fWireMin.resize(nplanes,-1);   
    fWireMax.resize(nplanes,-1);    
    fTimeMin.resize(nplanes,-1);    
    fTimeMax.resize(nplanes,-1);    
    fRawCharge.resize(nplanes,0);   
    fConvertedCharge.resize(nplanes,0);
      
    for(size_t p = 0; p < nplanes; ++p){
      for(size_t w = 0; w < geo->Plane(p).Nwires(); ++w){
	uint32_t channel = geo->PlaneWireToChannel(p, w);
	if( cf.BadChannel(channel) ) fBadChannels.push_back(channel);
      }
    }

  }

  //......................................................................
  RawDataDrawer::~RawDataDrawer()
  { 
  }


  //......................................................................
  void RawDataDrawer::RawDigit2D(const art::Event& evt,
				 evdb::View2D*     view,
				 unsigned int      plane) 
  {
    // Check if we're supposed to draw raw hits at all
    art::ServiceHandle<evd::RawDrawingOptions> drawopt;
    if (drawopt->fDrawRawDataOrCalibWires == 1) return;

    art::ServiceHandle<geo::Geometry> geo;
    art::ServiceHandle<evd::ColorDrawingOptions> cst;

    art::ServiceHandle<util::LArProperties> larp;
    art::ServiceHandle<util::DetectorProperties> detp;
    
    geo::PlaneID pid(drawopt->fCryostat, drawopt->fTPC, plane);

    unsigned int w  = 0;
    fRawCharge[plane]       = 0;
    fConvertedCharge[plane] = 0;
    
    int ticksPerPoint = drawopt->fTicksPerPoint;

    // to make det independent later:
    double mint = 5000., maxt = 0., minw = 5000., maxw = 0.;
        
    art::PtrVector<raw::RawDigit> rawhits;
    this->GetRawDigits(evt, rawhits);

    if(rawhits.size() < 1) return;

    for (auto const& hit : rawhits) {

      geo::SigType_t sigType = geo->SignalType(hit->Channel());
      geo::View_t    v       = geo->View(hit->Channel());
      std::vector<geo::WireID> wireids = geo->ChannelToWire(hit->Channel());
      bool skipchan = true;
      for(auto const& wid : wireids){
	// check that the plane and tpc are the correct ones to draw
	if(wid.planeID() == pid){
	  skipchan = false;
	}
      }
      if (skipchan) continue;
      
      std::vector<short> uncompressed(hit->Samples());
      raw::Uncompress(hit->ADCs(), uncompressed, hit->Compression());
      
      for(auto const& wid : wireids){
	// check that the plane and tpc are the correct ones to draw
	if(wid.planeID() == pid){
	  
	  double wire = 1.*wid.Wire;
	  double tick = 0;
	  // get an iterator over the adc values
	  std::vector<short>::iterator itr = uncompressed.begin();
	  while( itr != uncompressed.end() ){
	    int ticksUsed = 0;
	    double tdcsum = 0.;
	    double adcsum = 0.;
	    while(ticksUsed < ticksPerPoint && itr != uncompressed.end()){
	      tdcsum  += tick;
	      adcsum  += (1.*(*itr)) - hit->GetPedestal();
	      ++ticksUsed;
	      tick += 1.;
	      itr++; // this advance of the iterator is sufficient for the external loop too
	    }
	    double adc = adcsum/ticksPerPoint;
	    double tdc = tdcsum/ticksPerPoint;
	    
	    if(std::abs(adc) < drawopt->fMinSignal) continue;
	
	    fRawCharge[plane] += std::abs(adc);
	    double dQdX = std::abs(adc)/geo->WirePitch(v)/detp->ElectronsToADC();
	    fConvertedCharge[plane] += larp->BirksCorrection(dQdX);
	
	    int    co = 0;
	    double sf = 1.;
	    double q0 = 1000.0;
	    
	    co = cst->RawQ(sigType).GetColor(adc);
	    if (drawopt->fScaleDigitsByCharge) {
	      sf = std::sqrt(adc/q0);
	      if (sf>1.0) sf = 1.0;
	    }
	    if(wire < minw)
	      minw = wire;
	    if(wire > maxw)
	      maxw = wire;
	    if(tdc < mint)
	      mint = tdc;
	    if(tdc > maxt)
	      maxt = tdc;
	
	    // don't draw boxes for tdc values that don't exist
	    if(tdc > fTicks || tdc < 0) continue;
    
	    if(drawopt->fAxisOrientation < 1){
	      TBox& b1 = view->AddBox(wire-sf*0.5,
				      tdc-sf*0.5*ticksPerPoint,
				      wire+sf*0.5,
				      tdc+sf*0.5*ticksPerPoint);
	      b1.SetFillStyle(1001);
	      b1.SetFillColor(co);    
	      b1.SetBit(kCannotPick);
	    }
	    else{
	      TBox &b1 = view->AddBox(tdc-sf*0.5*ticksPerPoint,
				      wire-sf*0.5,
				      tdc+sf*0.5*ticksPerPoint,
				      wire+sf*0.5);
	      b1.SetFillStyle(1001);
	      b1.SetFillColor(co);    
	      b1.SetBit(kCannotPick);
	    }	  
	    
	    // 	  TBox& b2 = view->AddBox(wire-0.1,tdc-0.1,wire+0.1,tdc+0.1);
	    // 	  b2.SetFillStyle(0);
	    // 	  b2.SetLineColor(15);
	    // 	  b2.SetBit(kCannotPick);
	 
	  }// end loop over samples 
	}// end if in the right plane
      }// end loopo over wireids
    }//end loop over raw hits

    fWireMin[plane] = minw;   
    fWireMax[plane] = maxw;    
    fTimeMin[plane] = mint;    
    fTimeMax[plane] = maxt; 
    
    // now loop over all the bad channels and set them to 0 adc
    for(size_t bc = 0; bc < fBadChannels.size(); ++bc){
      
      geo::SigType_t sigType = geo->SignalType(fBadChannels[bc]);
	
      std::vector<geo::WireID> wireids = geo->ChannelToWire(fBadChannels[bc]);
      
      // check this is the correct plane and tpc
      for( auto const& wid : wireids){
	if(wid.planeID() == pid){
	
	  if(drawopt->fMinSignal > 0) continue;
	
	  int      co = cst->RawQ(sigType).GetColor(0);
	  double wire = 1.*w;
	  
	  for(int i = 0; i < fTicks; i += ticksPerPoint){
	    double tdc = i + 0.5*ticksPerPoint;
	    
	    if(drawopt->fAxisOrientation < 1){
	      TBox& b1 = view->AddBox(wire-0.5,tdc-0.5*ticksPerPoint,wire+0.5,tdc+0.5*ticksPerPoint);
	      b1.SetFillStyle(1001);
	      b1.SetFillColor(co);    
	      b1.SetBit(kCannotPick);
	    }
	    else{
	      TBox &b1 = view->AddBox(tdc-0.5*ticksPerPoint,wire-0.5,tdc+0.5*ticksPerPoint,wire+0.5);
	      b1.SetFillStyle(1001);
	      b1.SetFillColor(co);    
	      b1.SetBit(kCannotPick);
	    }	  
	  }
	}// end if in the right plane
      }// end loop over wireids
    }// end loop over bad channels    
	
  }

  //........................................................................
  int RawDataDrawer::GetRegionOfInterest(int plane,int& minw,int& maxw,int& mint,int& maxt)
  {
    art::ServiceHandle<geo::Geometry> geo;
 
    if((unsigned int)plane>fWireMin.size())
      {mf::LogWarning  ("RawDataDrawer") << " Requested plane " << plane <<" is larger than those available " << std::endl;
	return -1;
      }
  
    minw=fWireMin[plane];
    maxw=fWireMax[plane];
    mint=fTimeMin[plane];
    maxt=fTimeMax[plane];
  
    //make values a bit larger, but make sure they don't go out of bounds 
    minw= (minw-30<0) ? 0 : minw-30;
    mint= (mint-10<0) ? 0 : mint-10;

    maxw= (maxw+10>(int)geo->Nwires(plane)) ? geo->Nwires(plane) : maxw+10;
    maxt= (maxt+10>TotalClockTicks()) ? TotalClockTicks() : maxt+10;
    
    return 0;
  }

  //......................................................................
  void RawDataDrawer::GetChargeSum(int plane,double& charge,double& convcharge)
  {
    charge=fRawCharge[plane]; 
    convcharge=fConvertedCharge[plane];    
    
  }

  //......................................................................
  void RawDataDrawer::FillQHisto(const art::Event& evt,
				 unsigned int      plane,
				 TH1F*             histo)
  {

    // Check if we're supposed to draw raw hits at all
    art::ServiceHandle<evd::RawDrawingOptions> drawopt;
    if (drawopt->fDrawRawDataOrCalibWires==1) return;

    art::ServiceHandle<geo::Geometry> geo;
  
    art::PtrVector<raw::RawDigit> rawhits;
    this->GetRawDigits(evt, rawhits);
    
    geo::PlaneID pid(drawopt->fCryostat, drawopt->fTPC, plane);

    for (auto const& hit : rawhits) {
      
      std::vector<geo::WireID> wireids = geo->ChannelToWire(hit->Channel());
      for(auto const& wid : wireids){
	// check that the plane and tpc are the correct ones to draw
	if(wid.planeID() == pid){

	  std::vector<short> uncompressed(hit->Samples());
	  raw::Uncompress(hit->ADCs(), uncompressed, hit->Compression());
      
	  for(unsigned int j = 0; j < uncompressed.size(); ++j)
	    histo->Fill(1.*uncompressed[j] - hit->GetPedestal());

	  // this channel is on the correct plane, don't double count the raw signal
	  // if there are more than one wids for the channel
	  break;
	}
      }// end loop over wids
    }//end loop over raw hits

  }

  //......................................................................
  void RawDataDrawer::FillTQHisto(const art::Event& evt,
				  unsigned int      plane,
				  unsigned int      wire,
				  TH1F*             histo)
  {

    // Check if we're supposed to draw raw hits at all
    art::ServiceHandle<evd::RawDrawingOptions> drawopt;
    if (drawopt->fDrawRawDataOrCalibWires==1) return;

    art::ServiceHandle<geo::Geometry> geo;

    art::PtrVector<raw::RawDigit> rawhits;
    this->GetRawDigits(evt, rawhits);

    if(rawhits.size() < 1) return;

    geo::PlaneID pid(drawopt->fCryostat, drawopt->fTPC, plane);

    for (auto const& hit : rawhits) {
      
      std::vector<geo::WireID> wireids = geo->ChannelToWire(hit->Channel());
      for(auto const& wid : wireids){
	// check that the plane and tpc are the correct ones to draw
	if(wid.planeID() == pid  && 
	   wid.Wire      == wire){
      
	  std::vector<short> uncompressed(hit->Samples());
	  raw::Uncompress(hit->ADCs(), uncompressed, hit->Compression());
      
	  for(unsigned int j = 0; j < uncompressed.size(); ++j)
	    histo->Fill(1.*j, 1.*uncompressed[j] - hit->GetPedestal());

	  // this channel is on the correct plane, don't double count the raw signal
	  // if there are more than one wids for the channel
	  break;
	}// end if channel is in the right location
      }// end loop over wire ids
    }//end loop over raw hits

  }

  //......................................................................

  //   void RawDataDrawer::RawDigit3D(const art::Event& evt,
  // 				 evdb::View3D*     view)
  //   {
  //     // Check if we're supposed to draw raw hits at all
  //     art::ServiceHandle<evd::RawDrawingOptions> drawopt;
  //     if (drawopt->fDrawRawOrCalibHits!=0) return;


  //     art::ServiceHandle<geo::Geometry> geom;

  //     HitTower tower;
  //     tower.fQscale = 0.01;

  //     for (unsigned int imod=0; imod<drawopt->fRawDigitModules.size(); ++imod) {
  //       const char* which = drawopt->fRawDigitModules[imod].c_str();

  //       art::PtrVector<raw::RawDigit> rawhits;
  //       this->GetRawDigits(evt, which, rawhits);

  //       for (unsigned int i=0; i<rawhits.size(); ++i) {
  // 	double t = 0;
  // 	double q = 0;
  // 	t = rawhits[i]->fTDC[0];
  // 	for (unsigned int j=0; j<rawhits[i]->NADC(); ++j) {
  // 	  q += rawhits[i]->ADC(j);
  // 	}
  // 	// Hack for now...
  // 	if (q<=0.0) q = 1+i%10;
      
  // 	// Get the cell geometry for the hit
  // 	int         iplane = cmap->GetPlane(rawhits[i].get());
  // 	int         icell  = cmap->GetCell(rawhits[i].get());
  // 	double      xyz[3];
  // 	double      dpos[3];
  // 	geo::View_t v;
  // 	geom->CellInfo(iplane, icell, &v, xyz, dpos);
      
  // 	switch (drawopt->fRawDigit3DStyle) {
  // 	case 1:
  // 	  //
  // 	  // Render digits as towers
  // 	  //
  // 	  if (v==geo::kX) {
  // 	    tower.AddHit(v, iplane, icell, xyz[0], xyz[2], q,  t);
  // 	  }
  // 	  else if (v==geo::kY) {
  // 	    tower.AddHit(v, iplane, icell, xyz[1], xyz[2], q, t);
  // 	  }
  // 	  else abort();
  // 	  break;
  // 	default:
  // 	  //
  // 	  // Render Digits as boxes 
  // 	  //
  // 	  TPolyLine3D& p = view->AddPolyLine3D(5,kGreen+1,1,2);
  // 	  double sf = std::sqrt(0.01*q);
  // 	  if (v==geo::kX) {
  // 	    double x1 = xyz[0] - sf*dpos[0];
  // 	    double x2 = xyz[0] + sf*dpos[0];
  // 	    double z1 = xyz[2] - sf*dpos[2];
  // 	    double z2 = xyz[2] + sf*dpos[2];
  // 	    p.SetPoint(0, x1, geom->DetHalfHeight(), z1);
  // 	    p.SetPoint(1, x2, geom->DetHalfHeight(), z1);
  // 	    p.SetPoint(2, x2, geom->DetHalfHeight(), z2);
  // 	    p.SetPoint(3, x1, geom->DetHalfHeight(), z2);
  // 	    p.SetPoint(4, x1, geom->DetHalfHeight(), z1);
  // 	  }
  // 	  else if (v==geo::kY) {
  // 	    double y1 = xyz[1] - sf*dpos[1];
  // 	    double y2 = xyz[1] + sf*dpos[1];
  // 	    double z1 = xyz[2] - sf*dpos[2];
  // 	    double z2 = xyz[2] + sf*dpos[2];
  // 	    p.SetPoint(0, geom->DetHalfWidth(), y1, z1);
  // 	    p.SetPoint(1, geom->DetHalfWidth(), y2, z1);
  // 	    p.SetPoint(2, geom->DetHalfWidth(), y2, z2);
  // 	    p.SetPoint(3, geom->DetHalfWidth(), y1, z2);
  // 	    p.SetPoint(4, geom->DetHalfWidth(), y1, z1);
  // 	  }
  // 	  else abort();
  // 	  break;
  // 	} // switch fRawDigit3DStyle    
  //       }//end loop over raw digits
  //     }// end loop over RawDigit modules
  
  //     // Render the towers for that style choice
  //     if (drawopt->fRawDigit3DStyle==1) tower.Draw(view);
  //   }

  //......................................................................    

  int RawDataDrawer::GetRawDigits(const art::Event&              evt,
				  art::PtrVector<raw::RawDigit>& rawhits)
  {
    rawhits.clear();
    art::PtrVector<raw::RawDigit> temp;

    art::ServiceHandle<evd::RawDrawingOptions> rawopt;
    art::Handle< std::vector<raw::RawDigit> >  rdcol;
    try{
      evt.getByLabel(rawopt->fRawDataLabel, rdcol);
      for(unsigned int i = 0; i < rdcol->size(); ++i){
	art::Ptr<raw::RawDigit> rd(rdcol, i);
	temp.push_back(rd);
      }
      temp.swap(rawhits);
    }
    catch(cet::exception& e){
      writeErrMsg("GetRawDigits", e);
    }
  
    return rawhits.size();
  }

}// namespace
////////////////////////////////////////////////////////////////////////
