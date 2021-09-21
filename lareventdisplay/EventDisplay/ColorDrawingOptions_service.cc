////////////////////////////////////////////////////////////////////////
/// \file ColorDrawingOptions_service.cc
///
/// \author  brebel@fnal.gov

// Framework includes
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

/// LArSoft includes
#include "lareventdisplay/EventDisplay/ColorDrawingOptions.h"

namespace evd{

  //......................................................................
  ColorDrawingOptions::ColorDrawingOptions(fhicl::ParameterSet const& pset)
    : evdb::Reconfigurable{pset}
    , fColorOrGray(pset.get< int                 >("ColorOrGrayScale"))
    , fRawDiv     (pset.get< std::vector<int>    >("RawDiv")          )
    , fRecoDiv    (pset.get< std::vector<int>    >("RecoDiv")         )
    , fRawQLow    (pset.get< std::vector<double> >("RawQLow")         )
    , fRawQHigh   (pset.get< std::vector<double> >("RawQHigh")        )
    , fRecoQLow   (pset.get< std::vector<double> >("RecoQLow")        )
    , fRecoQHigh  (pset.get< std::vector<double> >("RecoQHigh")       )
  {
    this->CheckInputVectorSizes();

    for(size_t i = 0; i < fRawDiv.size(); ++i){
      fColorScaleRaw.push_back(evdb::ColorScale(fRawQLow[i], fRawQHigh[i],
						evdb::kBlueToRedII, evdb::kLinear,
						fRawDiv[i],
						285.0, 135.0, // angle in the color wheel
						0.65, 0.25)); // intensity from light to dark,
                                                              // starting with low color wheel value

      fGrayScaleRaw.push_back(evdb::ColorScale(fRawQLow[i], fRawQHigh[i],
					       evdb::kLinGray, evdb::kLinear,
					       fRawDiv[i],
					       270.0, 0.0, // angle in the color wheel
					       0.5, 0.5)); // intensity from light to dark,
                                                           // starting with low color wheel value
    }

    for(size_t i = 0; i < fRecoDiv.size(); ++i){
      fColorScaleReco.push_back(evdb::ColorScale(fRecoQLow[i], fRecoQHigh[i],
						 evdb::kBlueToRedII, evdb::kLinear,
						 fRecoDiv[i],
						 285.0, 135.0,
						 0.65, 0.25));
      fGrayScaleReco.push_back(evdb::ColorScale(fRecoQLow[i], fRecoQHigh[i],
						evdb::kLinGray, evdb::kLinear,
						fRecoDiv[i],
						270.0, 0.0,
						0.5, 0.5));
    }
  }

  //......................................................................
  void ColorDrawingOptions::CheckInputVectorSizes()
  {

    // compare all input vectors for reco and raw color scaling against
    // the number of possible signal types in the geometry
    if(fRawDiv.size() != geo::kMysteryType){
      if(fRawDiv.size() == 1 && fRawQLow.size() == 1 && fRawQHigh.size() == 1){
	// pad out the vectors to all have the same entries
	fRawDiv  .resize(geo::kMysteryType, fRawDiv[0]);
	fRawQLow .resize(geo::kMysteryType, fRawQLow[0]);
	fRawQHigh.resize(geo::kMysteryType, fRawQHigh[0]);
	mf::LogWarning("ColorDrawingOptions") << "only 1 value given for raw color scale: "
					      << "number of divisions, low and high values.\n"
					      << "Pad out those values for the number of signal types.";
      }
      else
	throw cet::exception("ColorDrawingOptionsUnclear") << "You have specified an incorrect number of "
							   << "values for the raw color/gray scale "
							   << "than there are types of signal planes in "
							   << "the detector.  It is unclear what your   "
							   << "intention is, so bail.\n";
    }// end check on the raw vector sizes

    if(fRecoDiv.size() != geo::kMysteryType){
      if(fRecoDiv.size() == 1 && fRecoQLow.size() == 1 && fRecoQHigh.size() == 1){
	// pad out the vectors to all have the same entries
	fRecoDiv  .resize(geo::kMysteryType, fRecoDiv[0]);
	fRecoQLow .resize(geo::kMysteryType, fRecoQLow[0]);
	fRecoQHigh.resize(geo::kMysteryType, fRecoQHigh[0]);
	mf::LogWarning("ColorDrawingOptions") << "only 1 value given for reco color scale: "
					      << "number of divisions, low and high values.\n"
					      << "Pad out those values for the number of signal types.";
      }
      else
	throw cet::exception("ColorDrawingOptionsUnclear") << "You have specified an incorrect number of "
							   << "values for the reco color/gray scale "
							   << "than there are types of signal planes in "
							   << "the detector.  It is unclear what your   "
							   << "intention is, so bail.\n";
    }// end check on the reco vector sizes

    return;
  }

  //......................................................................
  void ColorDrawingOptions::reconfigure(fhicl::ParameterSet const& pset)
  {
    fColorOrGray = pset.get< int    >("ColorOrGrayScale");
    fRawDiv      = pset.get< std::vector<int>    >("RawDiv");
    fRecoDiv     = pset.get< std::vector<int>    >("RecoDiv");
    fRawQLow     = pset.get< std::vector<double> >("RawQLow");
    fRawQHigh    = pset.get< std::vector<double> >("RawQHigh");
    fRecoQLow    = pset.get< std::vector<double> >("RecoQLow");
    fRecoQHigh   = pset.get< std::vector<double> >("RecoQHigh");

    this->CheckInputVectorSizes();

    for(size_t i = 0; i < fRawDiv.size(); ++i){
      fColorScaleRaw[i].SetBounds(fRawQLow[i],  fRawQHigh[i]);
      fGrayScaleRaw[i] .SetBounds(fRawQLow[i],  fRawQHigh[i]);
    }

    for(size_t i = 0; i < fRecoDiv.size(); ++i){
      fColorScaleReco[i].SetBounds(fRecoQLow[i], fRecoQHigh[i]);
      fGrayScaleReco[i] .SetBounds(fRecoQLow[i], fRecoQHigh[i]);
    }

    return;
  }

  //......................................................................
  const evdb::ColorScale& ColorDrawingOptions::RawQ(geo::SigType_t st) const
  {
    size_t pos = (size_t)st;

    if(st == geo::kMysteryType)
      throw cet::exception("ColorDrawingOptions") << "asked for RawQ with geo::kMysteryType, "
						  << "bad things will happen, so bail\n";

    if(fColorOrGray > 0) return fGrayScaleRaw[pos];

    return fColorScaleRaw[pos];
  }

  //......................................................................
  const evdb::ColorScale& ColorDrawingOptions::CalQ(geo::SigType_t st) const
  {
    size_t pos = (size_t)st;

    if(st == geo::kMysteryType)
      throw cet::exception("ColorDrawingOptions") << "asked for CalQ with geo::kMysteryType, "
						  << "bad things will happen, so bail\n";

    if(fColorOrGray > 0) return fGrayScaleReco[pos];

    return fColorScaleReco[pos];
  }

  //......................................................................
  const evdb::ColorScale& ColorDrawingOptions::RawT(geo::SigType_t st) const
  {
    size_t pos = (size_t)st;

    if(st == geo::kMysteryType)
      throw cet::exception("ColorDrawingOptions") << "asked for RawT with geo::kMysteryType, "
						  << "bad things will happen, so bail\n";


    if(fColorOrGray > 0) return fGrayScaleRaw[pos];

    return fColorScaleRaw[pos];
  }

  //......................................................................
  const evdb::ColorScale& ColorDrawingOptions::CalT(geo::SigType_t st) const
  {
    size_t pos = (size_t)st;

    if(st == geo::kMysteryType)
      throw cet::exception("ColorDrawingOptions") << "asked for CalT with geo::kMysteryType, "
						  << "bad things will happen, so bail\n";

    if(fColorOrGray > 0) return fGrayScaleReco[pos];

    return fColorScaleReco[pos];
  }
}// namespace

namespace evd {

  DEFINE_ART_SERVICE(ColorDrawingOptions)

} // namespace evd
////////////////////////////////////////////////////////////////////////
