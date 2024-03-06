///////////////////////////////////////////////////////////////////////
///
/// \file   ISpacePoint3D.h
///
/// \brief  This provides an interface for tools which are tasked with
///         drawing the simulated 3D objects
///
/// \author T. Usher
///
////////////////////////////////////////////////////////////////////////

#ifndef ISpacePoints3D_H
#define ISpacePoints3D_H

#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"

namespace evdb {
  class View3D;
}

namespace recob {
  class Hit;
  class SpacePoint;
}

namespace evdb_tool {
  class ISpacePoints3D {
  public:
    virtual ~ISpacePoints3D() noexcept = default;

    virtual void Draw(const std::vector<art::Ptr<recob::SpacePoint>>&, // Space points
                      evdb::View3D*,                                   // 3D display
                      int = 1,                                         // Color
                      int = 1,                                         // Marker
                      float = 1.,                                      // Size
                      const art::FindManyP<recob::Hit>* = nullptr      // pointer
    ) const = 0;
  };
}

#endif
