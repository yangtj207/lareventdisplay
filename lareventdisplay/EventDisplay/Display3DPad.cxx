///
/// \file    Display3DPad.cxx
/// \brief   Drawing pad showing a 3D rendering of the detector
/// \author  messier@indiana.edu
///
#include "TPad.h"
#include "TView3D.h"

#include "lareventdisplay/EventDisplay/Display3DPad.h"
#include "nutools/EventDisplayBase/View3D.h"
#include "nutools/EventDisplayBase/EventHolder.h"
#include "larcore/Geometry/Geometry.h"
#include "lareventdisplay/EventDisplay/RecoBaseDrawer.h"
#include "lareventdisplay/EventDisplay/SimulationDrawingOptions.h"
#include "lareventdisplay/EventDisplay/RecoDrawingOptions.h"
#include "lareventdisplay/EventDisplay/SimDrawers/ISim3DDrawer.h"
#include "lareventdisplay/EventDisplay/ExptDrawers/IExperimentDrawer.h"
#include "lareventdisplay/EventDisplay/3DDrawers/I3DDrawer.h"

#include "art/Framework/Principal/fwd.h"
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
    
    // Set up the 3D drawing tools for the simulation
    art::ServiceHandle<evd::SimulationDrawingOptions> simDrawOpt;

    // Implement the tools for handling the 3D simulation drawing tools
    const fhicl::ParameterSet& drawSim3DTools = simDrawOpt->f3DDrawerParams;
    
    for(const std::string& draw3DTool : drawSim3DTools.get_pset_names())
    {
        const fhicl::ParameterSet& draw3DToolParamSet = drawSim3DTools.get<fhicl::ParameterSet>(draw3DTool);
        
        fSim3DDrawerVec.push_back(art::make_tool<evdb_tool::ISim3DDrawer>(draw3DToolParamSet));
    }
    
    // Set up the 3D drawing tools for the reconstruction
    art::ServiceHandle<evd::RecoDrawingOptions> recoDrawOpt;

    // Implement the tools for handling the 3D reco drawing tools
    const fhicl::ParameterSet& drawReco3DTools = recoDrawOpt->f3DDrawerParams;
    
    for(const std::string& draw3DTool : drawReco3DTools.get_pset_names())
    {
        const fhicl::ParameterSet& draw3DToolParamSet = drawReco3DTools.get<fhicl::ParameterSet>(draw3DTool);
        
        fReco3DDrawerVec.push_back(art::make_tool<evdb_tool::I3DDrawer>(draw3DToolParamSet));
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
        this->RecoBaseDraw()->  Slice3D      (*evt, fView);

        // Call the 3D simulation drawing tools
        for(auto& draw3D : fSim3DDrawerVec) draw3D->Draw(*evt, fView);
        
        // Call the 3D reco drawing tools
        for(auto& draw3D : fReco3DDrawerVec) draw3D->Draw(*evt, fView);
        
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


}//namespace
