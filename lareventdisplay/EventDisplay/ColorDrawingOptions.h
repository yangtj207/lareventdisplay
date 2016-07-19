/// \file  ColorDrawingOptions.h
/// \brief The color scales used by the event display
/// \author messier@indiana.edu
/// \version $Id: ColorScaleTable.h,v 1.1.1.1 2010/11/10 19:44:54 p-novaart Exp $
#ifndef EVD_COLORDRAWINGOPTIONS_H
#define EVD_COLORDRAWINGOPTIONS_H
#ifndef __CINT__
#include "nutools/EventDisplayBase/ColorScale.h"
#include "larcore/Geometry/Geometry.h"

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"

namespace evd {
  class ColorDrawingOptions {
  public:

    ColorDrawingOptions(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
    ~ColorDrawingOptions();

    void reconfigure(fhicl::ParameterSet const& pset);

    const evdb::ColorScale& RawQ(geo::SigType_t st) const;
    const evdb::ColorScale& CalQ(geo::SigType_t st) const;
    const evdb::ColorScale& RawT(geo::SigType_t st) const;
    const evdb::ColorScale& CalT(geo::SigType_t st) const;

    int                 fColorOrGray; ///< 0 = color, 1 = gray
    std::vector<int   > fRawDiv;      ///< number of divisions in raw
    std::vector<int   > fRecoDiv;     ///< number of divisions in raw
    std::vector<double> fRawQLow;     ///< low  edge of ADC values for drawing raw digits
    std::vector<double> fRawQHigh;    ///< high edge of ADC values for drawing raw digits
    std::vector<double> fRecoQLow;    ///< low  edge of ADC values for drawing raw digits
    std::vector<double> fRecoQHigh;   ///< high edge of ADC values for drawing raw digits

  private:

    void CheckInputVectorSizes();

    std::vector<evdb::ColorScale> fColorScaleRaw;
    std::vector<evdb::ColorScale> fColorScaleReco;
    std::vector<evdb::ColorScale> fGrayScaleRaw;
    std::vector<evdb::ColorScale> fGrayScaleReco;
  };
}
#endif // __CINT__
DECLARE_ART_SERVICE(evd::ColorDrawingOptions, LEGACY)
#endif
////////////////////////////////////////////////////////////////////////
