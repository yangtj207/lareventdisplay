/// \file    GeometryDrawer.h
/// \brief   Class to aid in the rendering of Geometry objects
/// \author  messier@indiana.edu
#ifndef EVD_GEOMETRYDRAWER_H
#define EVD_GEOMETRYDRAWER_H
#include <vector>
class TH1F;
namespace evdb{ 
  class View2D;
  class View3D;
}

namespace evd {
  /// Aid in the rendering of Geometry objects
  class GeometryDrawer {
  public:
    GeometryDrawer();
    ~GeometryDrawer();
    void DetOutline3D(evdb::View3D* view);

  };
}

#endif
////////////////////////////////////////////////////////////////////////
