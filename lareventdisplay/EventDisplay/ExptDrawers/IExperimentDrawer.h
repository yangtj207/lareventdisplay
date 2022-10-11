///////////////////////////////////////////////////////////////////////
///
/// \file   IExperimentDrawer.h
///
/// \brief  This is the interface class for drawing 3D detector geometries
///
/// \author T. Usher
///
////////////////////////////////////////////////////////////////////////

#ifndef IExperimentDrawer_H
#define IExperimentDrawer_H

namespace evdb {
  class View3D;
}

namespace evd_tool {
  class IExperimentDrawer {
  public:
    virtual ~IExperimentDrawer() noexcept = default;

    virtual void DetOutline3D(evdb::View3D* view) = 0;
  };
}

#endif
