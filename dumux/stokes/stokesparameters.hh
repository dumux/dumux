// $Id$

#ifndef DUNE_STOKESPARAMETERS_HH
#define DUNE_STOKESPARAMETERS_HH
#include<iomanip>

#include<dune/common/fvector.hh>
#include<dune/common/fmatrix.hh>
#include<dune/common/exceptions.hh>
#include<dune/grid/common/grid.hh>
#include<dune/grid/common/referenceelements.hh>


/*
  Stokes Equation

  -(\mu \Delta u)+ \nabla p = f in \Omega

  \nabla \cdot = 0 in \Omega

  u= g on \partial \Omega

*/

class DGStokesParameters
{
public:
    DGStokesParameters() : mu(1.0), sigma(0.0), epsilon(1),dg_method("NIPG")
    {}
    //coeff of viscosity
    double mu;
    // sigma and epsilon are paremetrs for defining the different versions of DG
    /*
      if sigma > 0 and epsilon = +1 --> NIPG
      if sigma > 0 and epsilon = -1 --> SIPG
      if sigma = 0 and epsilon = +1 --> OBB
    */
    double sigma;
    int epsilon;
    std::string dg_method;
};



#endif
