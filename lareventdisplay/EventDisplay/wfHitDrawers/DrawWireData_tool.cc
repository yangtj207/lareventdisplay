////////////////////////////////////////////////////////////////////////
/// \file   DrawWireData_tool.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "lareventdisplay/EventDisplay/wfHitDrawers/IWFWireDrawer.h"
#include "art/Utilities/ToolMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "nutools/EventDisplayBase/EventHolder.h"
#include "lareventdisplay/EventDisplay/RecoDrawingOptions.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Wire.h"

#include "TF1.h"
#include "TPolyLine.h"

#include <fstream>

namespace evdb_tool
{

class DrawWireData : public IWFWireDrawer
{
public:
    explicit DrawWireData(const fhicl::ParameterSet& pset);
    
    ~DrawWireData();
    
    void configure(const fhicl::ParameterSet& pset)       override;
    void Draw(evdb::View2D&, raw::ChannelID_t&)     const override;
    
private:
    
    std::vector<int> fColorMap;
};
    
//----------------------------------------------------------------------
// Constructor.
DrawWireData::DrawWireData(const fhicl::ParameterSet& pset)
{
    configure(pset);
}
    
DrawWireData::~DrawWireData()
{
}
    
void DrawWireData::configure(const fhicl::ParameterSet& pset)
{
    fColorMap.push_back(kBlue);
    fColorMap.push_back(kMagenta);
    fColorMap.push_back(kBlack);
    fColorMap.push_back(kRed);
    
    return;
}

    
void DrawWireData::Draw(evdb::View2D&     view2D,
                        raw::ChannelID_t& channel) const
{
    art::ServiceHandle<evd::RecoDrawingOptions> recoOpt;
    
    //grab the singleton with the event
    const art::Event* event = evdb::EventHolder::Instance()->GetEvent();
    if(!event) return;
    
    for (size_t imod = 0; imod < recoOpt->fWireLabels.size(); ++imod)
    {
        // Step one is to recover the hits for this label that match the input channel
        art::InputTag const which = recoOpt->fWireLabels[imod];
        
        art::Handle< std::vector<recob::Wire> > wireVecHandle;
        event->getByLabel(which, wireVecHandle);
        
        for(size_t wireIdx = 0; wireIdx < wireVecHandle->size(); wireIdx++)
        {
            art::Ptr<recob::Wire> wire(wireVecHandle, wireIdx);
            
            if (wire->Channel() != channel) continue;

            // Recover a full wire version of the deconvolved wire data
            // (the ROIs don't tend to display well)
            std::vector<float> signal = wire->Signal();

            TPolyLine& wireWaveform = view2D.AddPolyLine(signal.size(), fColorMap[imod % fColorMap.size()], 2, 1);
                
            for(size_t idx = 0; idx < signal.size(); idx++) wireWaveform.SetPoint(idx,float(idx)+0.5,signal[idx]);
                
            wireWaveform.Draw("same");
        }
    }//end loop over HitFinding modules

    return;
}
    
DEFINE_ART_CLASS_TOOL(DrawWireData)
}
