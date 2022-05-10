///////////////////////////////////////////////////////////////////////
///
/// \file   IWFHitDrawers.h
///
/// \brief  This provides an interface for tools which are tasked with
///         drawing the reconstructed hits on waveforms
///
/// \author T. Usher
///
////////////////////////////////////////////////////////////////////////

#ifndef IWFHitDrawer_H
#define IWFHitDrawer_H

namespace fhicl { class ParameterSet; }
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
namespace evdb { class View2D; }

namespace evdb_tool
{
    class IWFHitDrawer
    {
    public:
        virtual ~IWFHitDrawer() noexcept = default;

        virtual void configure(const fhicl::ParameterSet& pset)       = 0;
        virtual void Draw(evdb::View2D&,
                          raw::ChannelID_t&)                    const = 0;
    };
}

#endif
