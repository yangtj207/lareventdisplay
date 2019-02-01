///////////////////////////////////////////////////////////////////////
///
/// \file   IWFWireDrawer.h
///
/// \brief  This provides an interface for tools which are tasked with
///         drawing the "wire" data (deconvolved waveforms)
///
/// \author T. Usher
///
////////////////////////////////////////////////////////////////////////

#ifndef IWFWireDrawer_H
#define IWFWireDrawer_H

#include "fhiclcpp/ParameterSet.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "nutools/EventDisplayBase/View2D.h"
#include "TF1.h"

namespace evdb_tool
{
    class IWFWireDrawer
    {
    public:
        virtual ~IWFWireDrawer() noexcept = default;
        
        virtual void configure(const fhicl::ParameterSet& pset)       = 0;
        virtual void Draw(evdb::View2D&, raw::ChannelID_t&)     const = 0;
    };
}

#endif
