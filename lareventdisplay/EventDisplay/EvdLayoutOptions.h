////////////////////////////////////////////////////////////////////////
//
// Display parameters for the raw data
//
// \author brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef EVDLAYOUTOPTIONS_H
#define EVDLAYOUTOPTIONS_H
#ifndef __CINT__
#include <string>

#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "nuevdb/EventDisplayBase/Reconfigurable.h"

namespace evd {
  class EvdLayoutOptions : public evdb::Reconfigurable
  {
  public:
    explicit EvdLayoutOptions(fhicl::ParameterSet const& pset);

    void reconfigure(fhicl::ParameterSet const& pset) override;

    fhicl::ParameterSet const& fParameterSet;

      int         fShowSideBar;		           ///< 1 to show, 0 don't show
      int         fAutoZoomInterest;           ///< Set the automatic zoom to the interest region
      int         fPrintTotalCharge;           ///< Print out the total charge in an event
      int         fShowEndPointSection;        ///< Show section corresponding to EndPoint finding
      int         fShowEndPointMarkers;        ///< Draw EndPoint Markers if clicked.
      int 	      fShowClusterSection; 	       ///< Show section to make clusters
      int	      fMakeClusters;		       ///< Draw two lines to make clusters if clicked
      int 	      fChangeWire; 		           ///< 1 to click mouse and change wire, 0 don't
      int         fEnableMCTruthCheckBox;      ///< 1 to have the check box appear, 0 otherwise

      bool        fThreeWindow;                ///< true to draw rectangular box representing 3 windows
      bool        fDrawGrid;                   ///< true to draw backing grid
      bool        fDrawAxes;                   ///< true to draw coordinate axes
      bool        fDrawBadChannels;            ///< true to draw bad channels

      std::string fDisplayName;                ///< Name to apply to 2D display
  };
}//namespace
#endif // __CINT__
DECLARE_ART_SERVICE(evd::EvdLayoutOptions, LEGACY)
#endif
