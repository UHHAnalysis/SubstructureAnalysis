#ifndef GroomingHists_H
#define GroomingHists_H

#include "SFrameTools/include/BaseHists.h"

class GroomingHists : public BaseHists {

public:
   /// Named constructor
   GroomingHists(const char* name);

   /// Default destructor
   ~GroomingHists();

   void Init();

   void Fill();

   void Finish();

private:

}; 


#endif 
