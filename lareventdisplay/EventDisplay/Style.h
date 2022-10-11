////////////////////////////////////////////////////////////////////////
/// \file Style.h
//
/// \author messier@indiana.edu
////////////////////////////////////////////////////////////////////////
#ifndef EVD_LINESTYLE_H
#define EVD_LINESTYLE_H
class TLine;

namespace evd {
  /// Parameters for drawing options. Allow a consistent style for
  /// drawing particle tracks
  class Style {
  public:
    static const char* LatexName(int pdgcode);
    static void FromPDG(TLine& line, int pdgcode);
    static int ColorFromPDG(int pdgcode);
    static int LineStyleFromPDG(int pdgcode);
    static int LineWidthFromPDG(int pdgcode);
  };
}
#endif
