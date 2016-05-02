#ifndef ComparisonHists_H
#define ComparisonHists_H

#include "SFrameTools/include/BaseHists.h"
#include "TH1.h"
/**
 *   Substructure class for booking and filling histograms
 *
 *   This class books and fills a collection of histograms.
 *   It should have a unique name, such that the histograms
 *   of multiple instances of this class are ordered in the
 *   output file. 
 *   Always sort your histograms and used methods topically.
 *   This example collection can be used for data and reconstructed
 *   MC events.
 *   
 *   @version $Revision: 1.2 $
 */

class ComparisonHists : public BaseHists {

public:
   /// Named constructor
   ComparisonHists(const char* name);

   /// Default destructor
   ~ComparisonHists();

   void Init();

   void Fill();

   void Finish();

   TH1F *m_scale_factors;

private:

}; // class ComparisonHists


#endif // ComparisonHists_H
