///////////////////////////////////////////////////////////////////////
///
/// \file   IWaveformDrawer.h
///
/// \brief  This provides an interface for tools which are tasked with
///         drawing the "wire" data (deconvolved waveforms)
///
/// \author T. Usher
///
////////////////////////////////////////////////////////////////////////

#ifndef IWaveformDrawer_H
#define IWaveformDrawer_H

#include "fhiclcpp/ParameterSet.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "nuevdb/EventDisplayBase/View2D.h"
#include "TF1.h"

namespace evdb_tool
{
    class IWaveformDrawer
    {
    public:
        virtual ~IWaveformDrawer() noexcept = default;

        virtual void configure(const fhicl::ParameterSet& pset)           = 0;
        virtual void Fill(evdb::View2D&, raw::ChannelID_t&, float, float) = 0;
        virtual void Draw(const std::string&,float,float)                 = 0;

        virtual float getMaximum() const                                  = 0;
        virtual float getMinimum() const                                  = 0;
    };
}

#endif
