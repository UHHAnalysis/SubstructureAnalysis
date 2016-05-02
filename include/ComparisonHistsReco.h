#ifndef ComparisonHistsReco_H
#define ComparisonHistsReco_H

#include "SFrameTools/include/BaseHists.h"

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

class ComparisonHistsReco : public BaseHists {

public:
   /// Named constructor
   ComparisonHistsReco(const char* name);

   /// Default destructor
   ~ComparisonHistsReco();

   void Init();

   void Fill();

   void Finish();

private:

}; // class ComparisonHistsReco


#endif // ComparisonHistsReco_H
