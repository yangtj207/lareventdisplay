///
/// \file    Display3DPad.cxx
/// \brief   Drawing pad showing a 3D rendering of the detector
/// \author  messier@indiana.edu
///
#include <iostream>
#include "TPad.h"
#include "TView3D.h"
#include "TGLViewer.h"
#include "TPolyLine3D.h"

#include "lareventdisplay/EventDisplay/Display3DPad.h"
#include "nutools/EventDisplayBase/View3D.h"
#include "nutools/EventDisplayBase/EventHolder.h"
#include "larcore/Geometry/Geometry.h"
#include "lareventdisplay/EventDisplay/SimulationDrawer.h"
#include "lareventdisplay/EventDisplay/RecoBaseDrawer.h"
#include "lareventdisplay/EventDisplay/EvdLayoutOptions.h"
#include "lareventdisplay/EventDisplay/SimulationDrawingOptions.h"
#include "lareventdisplay/EventDisplay/HitSelector.h"
#include "lareventdisplay/EventDisplay/SimDrawers/ISim3DDrawer.h"
#include "lareventdisplay/EventDisplay/ExptDrawers/IExperimentDrawer.h"
#include "lardataobj/RecoBase/Seed.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Utilities/make_tool.h"

namespace evd{

///
/// Create a pad to show a 3D rendering of the detector and events
/// @param nm : Name of the pad
/// @param ti : Title of the pad
/// @param x1 : Location of left  edge of pad (0-1)
/// @param x2 : Location of right edge of pad (0-1)
/// @param y1 : Location of bottom edge of pad (0-1)
/// @param y2 : Location of top    edge of pad (0-1)
/// @param opt: Options. Currently just a place holder
///
Display3DPad::Display3DPad(const char* nm, const char* ti,
                           double      x1, double      y1,
                           double      x2, double      y2,
                           const char* /*opt*/) :
    DrawingPad(nm, ti, x1, y1, x2, y2)
{
    this->Pad()->SetFillColor(kBlack);
    this->Pad()->Draw();
    this->Pad()->cd();
    fView = new evdb::View3D();
    
    // Set up the 3D drawing tools
    art::ServiceHandle<evd::SimulationDrawingOptions> simDrawOpt;

    // Implement the tools for handling the responses
    const fhicl::ParameterSet& draw3DTools = simDrawOpt->f3DDrawerParams;
    
    for(const std::string& draw3DTool : draw3DTools.get_pset_names())
    {
        const fhicl::ParameterSet& draw3DToolParamSet = draw3DTools.get<fhicl::ParameterSet>(draw3DTool);
        
        fSim3DDrawerVec.push_back(art::make_tool<evdb_tool::ISim3DDrawer>(draw3DToolParamSet));
    }
    
    return;
}

//......................................................................

Display3DPad::~Display3DPad()
{
    if (fView) { delete fView; fView = 0; }
}

//......................................................................

void Display3DPad::Draw()
{
    fView->Clear();

    art::ServiceHandle<geo::Geometry> geo;

    // grab the event from the singleton
    const art::Event *evt = evdb::EventHolder::Instance()->GetEvent();

    if(evt){
        this->GeometryDraw()->DetOutline3D(fView);
//        this->SimulationDraw()->MCTruth3D    (*evt, fView);
        this->RecoBaseDraw()->  PFParticle3D (*evt, fView);
        this->RecoBaseDraw()->  Edge3D       (*evt, fView);
        this->RecoBaseDraw()->  SpacePoint3D (*evt, fView);
        this->RecoBaseDraw()->  Prong3D      (*evt, fView);
        this->RecoBaseDraw()->  Seed3D       (*evt, fView);
        this->RecoBaseDraw()->  Vertex3D     (*evt, fView);
        this->RecoBaseDraw()->  Event3D      (*evt, fView);
        
        for(auto& draw3D : fSim3DDrawerVec) draw3D->Draw(*evt, fView);
        
        art::ServiceHandle<evd::EvdLayoutOptions> evdlayoutoptions;
        if(evdlayoutoptions->fMakeSeeds) UpdateSeedCurve();
    }

    this->Pad()->Clear();
    this->Pad()->cd();
    if (fPad->GetView()==0) {
        int irep;
        double rmin[]={-2.1*geo->DetHalfWidth(),-2.1*geo->DetHalfHeight(),-0.5*geo->DetLength()};
        double rmax[]={ 2.1*geo->DetHalfWidth(), 2.1*geo->DetHalfHeight(), 0.5*geo->DetLength()};
        TView3D* v = new TView3D(1,rmin,rmax);
        v->SetPerspective();
        v->SetView(0.0,260.0,270.0,irep);
        fPad->SetView(v); // ROOT takes ownership of object *v
    }
    fView->Draw();
    fPad->Update();
}


////////////////////////////////////////////////////////////////////////


void Display3DPad::UpdateSeedCurve()
{
    std::vector<recob::Seed> SeedVec = this->HitSelectorGet()->SeedVector();
    std::cout<<"Seeds available to Display3DPad : " << SeedVec.size()<<std::endl;
    trkf::BezierTrack BTrack(SeedVec);
    
    int N=100;
    TPolyLine3D& pl = fView->AddPolyLine3D(N,kOrange+5,2,0);
    fView->Draw();
    
    
    double  xyzpt[3] ;
    
    for(int i=0; i!=N; i++)
    {
        BTrack.GetTrackPoint(float(i)/N,xyzpt );
        double x = xyzpt[0];
        double y = xyzpt[1];
        double z = xyzpt[2];
        
        pl.SetPoint(i,x,y,z);
        
    }
    
    fView->Draw();
  //  UpdatePad();
}

}//namespace
