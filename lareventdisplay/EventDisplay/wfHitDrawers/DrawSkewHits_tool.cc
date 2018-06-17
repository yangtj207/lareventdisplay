////////////////////////////////////////////////////////////////////////
/// \file   DrawSkewHits_tool.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "lareventdisplay/EventDisplay/wfHitDrawers/IWFHitDrawer.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

#include <fstream>

namespace evdb_tool
{

class DrawSkewHits : public IWFHitDrawer
{
public:
    explicit DrawSkewHits(const fhicl::ParameterSet& pset);
    
    ~DrawSkewHits();
    
    void configure(const fhicl::ParameterSet& pset)                                       override;
    void Draw(evdb::View2D&, std::vector<std::unique_ptr<TF1>>&, raw::ChannelID_t&) const override;
    
private:
};
    
//----------------------------------------------------------------------
// Constructor.
DrawSkewHits::DrawSkewHits(const fhicl::ParameterSet& pset)
{
    configure(pset);
}
    
DrawSkewHits::~DrawSkewHits()
{
}
    
void DrawSkewHits::configure(const fhicl::ParameterSet& pset)
{
    return;
}

    
    void DrawSkewHits::Draw(evdb::View2D&                      view2D,
                            std::vector<std::unique_ptr<TF1>>& hitFuncVec,
                            raw::ChannelID_t&                  channel) const
{
    return;
}
    
DEFINE_ART_CLASS_TOOL(DrawSkewHits)
}
