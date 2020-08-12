#include "CutsCountMSSS.h"
#include "ConfigSvc.h"

CutsCountMSSS::CutsCountMSSS(EventBase* eventBase) {
  m_event = eventBase;
}

CutsCountMSSS::~CutsCountMSSS() {
  
}

// Function that lists all of the common cuts for this Analysis
bool CutsCountMSSS::CountMSSSCutsOK() {
  // List of common cuts for this analysis into one cut
  return true;
}

