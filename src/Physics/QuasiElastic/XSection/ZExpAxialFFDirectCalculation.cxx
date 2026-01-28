//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 

\author   Sarah Hawkins
	  based off ZExpAxialFormFactorModel by Aaron Meyer <asmeyer2012 \at uchicago.edu>

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include <TMath.h>
#include <sstream>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Interaction/Interaction.h"
#include "Physics/QuasiElastic/XSection/ZExpAxialFFDirectCalculation.h"
#include "Framework/Messenger/Messenger.h"

using std::ostringstream;

using namespace genie;

//____________________________________________________________________________
ZExpAxialFFDirectCalculation::ZExpAxialFFDirectCalculation() :
AxialFormFactorModelI("genie::ZExpAxialFFDirectCalculation")
{

}
//____________________________________________________________________________
ZExpAxialFFDirectCalculation::ZExpAxialFFDirectCalculation(string config) :
AxialFormFactorModelI("genie::ZExpAxialFFDirectCalculation", config)
{

}
//____________________________________________________________________________
ZExpAxialFFDirectCalculation::~ZExpAxialFFDirectCalculation()
{
  delete[] fZ_An;
}
//____________________________________________________________________________
double ZExpAxialFFDirectCalculation::FA(const Interaction * interaction) const
{
  // calculate and return FA
  double q2 = interaction->KinePtr()->q2();
  double zparam = this->CalculateZ(q2);
  if (zparam != zparam) // checks for nan
  {
    LOG("ZExpAxialFFDirectCalculation",pWARN) << "Undefined expansion parameter";
    return 0.;
  }
  double fa = 0.;
  for (int ki=0;ki<=fKmax;ki++)
  {
    fa = fa + TMath::Power(zparam,ki) * fZ_An[ki];
  
  }

  return fa;
}
//___________________________________________________________________________
double ZExpAxialFFDirectCalculation::CalculateFA0(void) const
{

  // calculate and return FA
  double zparam = this->CalculateZ(0.);
  if (zparam != zparam) // checks for nan
  {
    LOG("ZExpAxialFFDirectCalculation",pWARN) << "Undefined expansion parameter";
    return 0.;
  }

  double fa0 = 0.;
  for (int ki=0;ki<=fKmax;ki++)
  {
    fa0 = fa0 + TMath::Power(zparam,ki) * fZ_An[ki];
  }
 
  return fa0;

}
//___________________________________________________________________________
double ZExpAxialFFDirectCalculation::CalculateZ(double q2) const
{

  // calculate z expansion parameter
  double znum  = TMath::Sqrt(fTcut - q2) - TMath::Sqrt(fTcut - fT0);
  double zden  = TMath::Sqrt(fTcut - q2) + TMath::Sqrt(fTcut - fT0);

  return znum/zden;
}
//____________________________________________________________________________
void ZExpAxialFFDirectCalculation::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void ZExpAxialFFDirectCalculation::Configure(string param_set)
{
  Algorithm::Configure(param_set);
  this->LoadConfig();
}
//____________________________________________________________________________
void ZExpAxialFFDirectCalculation::LoadConfig(void)
{
// get config options from the configuration registry or set defaults 
// from the global parameter list
  GetParam( "QEL-Kmax", fKmax ) ;
  GetParam( "QEL-T0", fT0 ) ;
  GetParam( "QEL-Tcut", fTcut ) ;

  assert(fKmax > 0);
  
  fZ_An = new double [fKmax]; 

  // load the user-defined coefficient values
  // all a_n for n<fKmax are explicit inputs
  for (int ip=0;ip<fKmax+1;ip++) {
    ostringstream alg_key;
    alg_key << "QEL-Z_A" << ip;
    GetParam( alg_key.str(), fZ_An[ip] ) ;
    
    // Add print statement
    std::cout << " fZ_An " << ip << " = " << fZ_An[ip] << std::endl;
  }

  double FA0 = this->CalculateFA0();

  std::cout << " FA0 = " << FA0 << std::endl;

}
//____________________________________________________________________________

