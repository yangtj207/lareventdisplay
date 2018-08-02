///////////////////////////////////////////////////////////////////////
///
/// \file   I3DDrawers.h
///
/// \brief  This provides an interface for tools which are tasked with
///         drawing the simulated 3D objects
///
/// \author T. Usher
///
////////////////////////////////////////////////////////////////////////

#ifndef I3DDrawer_H
#define I3DDrawer_H

#include "fhiclcpp/ParameterSet.h"
#include "nutools/EventDisplayBase/View3D.h"
#include "art/Framework/Principal/Event.h"

namespace evdb_tool
{
    class I3DDrawer
    {
    public:
        virtual ~I3DDrawer() noexcept = default;
        
        virtual void Draw(const art::Event&, evdb::View3D*) const = 0;
    };
}

#endif
