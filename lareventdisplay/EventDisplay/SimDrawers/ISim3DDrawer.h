///////////////////////////////////////////////////////////////////////
///
/// \file   ISim3DDrawers.h
///
/// \brief  This provides an interface for tools which are tasked with
///         drawing the simulated 3D objects
///
/// \author T. Usher
///
////////////////////////////////////////////////////////////////////////

#ifndef ISim3DDrawer_H
#define ISim3DDrawer_H

namespace art { class Event; }
namespace evdb { class View3D; }

namespace evdb_tool
{
    class ISim3DDrawer
    {
    public:
        virtual ~ISim3DDrawer() noexcept = default;

        virtual void Draw(const art::Event&, evdb::View3D*) const = 0;
    };
}

#endif
