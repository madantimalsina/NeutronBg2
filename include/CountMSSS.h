#ifndef CountMSSS_H
#define CountMSSS_H

#include "Analysis.h"

#include "EventBase.h"
#include "CutsCountMSSS.h"

#include <TTreeReader.h>
#include <string>

class CountMSSS: public Analysis {

public:
  CountMSSS(); 
  ~CountMSSS();

  void Initialize();
  void Execute();
  void Finalize();
   int EventNumber;
  std::string ParentParticle;


protected:
  CutsCountMSSS* m_cutsCountMSSS;
  ofstream outfile;
  ConfigSvc* m_conf;
};

#endif
