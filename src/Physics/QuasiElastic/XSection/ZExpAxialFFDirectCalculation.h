//___________________________________________________________________________
/* !
\brief    Directly calculates axial form factor using the model-independent z-expansion technique
          Relies on a0...a_kmax coefficients specified in the ModelConfiguration
          Adapted from genie::ZExpAxialFormFactorModel

\author   Sarah Hawkins 

          Based off ZExpAxialFormFactorModel by
          Aaron Meyer <asmeyer2012 \at uchicago.edu>


\created  January 2026

\cpright  Copyright (c) 2003-2025, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org

*/
//____________________________________________________________________________

#ifndef _Z_EXPANSION_AXIAL_FF_DIRECT_CALCULATION_H_
#define _Z_EXPANSION_AXIAL_FF_DIRECT_CALCULATION_H_

#include "Physics/QuasiElastic/XSection/AxialFormFactorModelI.h"

namespace genie {

class ZExpAxialFFDirectCalculation : public AxialFormFactorModelI {

public:
  ZExpAxialFFDirectCalculation();
  ZExpAxialFFDirectCalculation(string config);
  virtual ~ZExpAxialFFDirectCalculation();

  // implement the AxialFormFactorModelI interface
  double FA (const Interaction * interaction) const;

  // overload Algorithm's Configure()
  void   Configure  (const Registry & config);
  void   Configure  (string param_set);

private:

  // calculate z parameter used in expansion
  double CalculateFA0(void) const;
  double CalculateZ(double q2) const;
  void LoadConfig(void);

  int    fKmax;
  double fT0;
  double fTcut;
  double* fZ_An;

};

}         // genie namespace

#endif    // _Z_EXPANSION_AXIAL_FF_DIRECT_CALCULATION_H_
