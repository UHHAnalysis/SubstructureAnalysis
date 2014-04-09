#ifndef SubstructureHistsNew_H
#define SubstructureHistsNew_H

#include "SFrameTools/include/AnalysisModule.h"

/**
 *   Substructure class for booking and filling histograms, the new version using AnalysisModule mechanisms.
 */

class SubstructureHistsNew: public Hists {
public:
    // use the same constructor arguments as Hists for forwarding:
    SubstructureHistsNew(Context & ctx, const string & dirname);

    virtual void fill(EventCalc & ev);
};


#endif
