#ifndef CutsCountMSSS_H
#define CutsCountMSSS_H

#include "EventBase.h"

class CutsCountMSSS  {

public:
  CutsCountMSSS(EventBase* eventBase);
  ~CutsCountMSSS();
  bool CountMSSSCutsOK();


private:
  
  EventBase* m_event;

};

#endif
