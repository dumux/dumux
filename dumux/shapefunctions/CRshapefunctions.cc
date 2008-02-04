#include "config.h"
#include "CRshapefunctions.hh"

namespace Dune {

  template<typename C, typename T, int d>
  CRCubeShapeFunctionSetContainer<C,T,d> CRShapeFunctions<C,T,d>::cube;

  template<typename C, typename T, int d>
  CRShapeFunctionSetContainer<C,T,d> CRShapeFunctions<C,T,d>::general;

  template<typename C, typename T, int d>
  CRSimplexShapeFunctionSetContainer<C,T,d> CRShapeFunctions<C,T,d>::simplex;

  namespace {

    template <class C, class T, int d>
    struct InitCRShapefunctions
    {
      CRCubeShapeFunctionSetContainer<C,T,d> & f1;
      CRSimplexShapeFunctionSetContainer<C,T,d> & f2;
      CRShapeFunctionSetContainer<C,T,d> & f3;
      InitCRShapefunctions() :
        f1(CRShapeFunctions<C,T,d>::cube),
        f2(CRShapeFunctions<C,T,d>::simplex),
        f3(CRShapeFunctions<C,T,d>::general)
        {
          InitCRShapefunctions<C,T,d-1> i;
        };
    };
  
    template <class C, class T>
    struct InitCRShapefunctions<T,C,0>
    {
      enum { d=0 };
      InitCRShapefunctions()
        {
        };
    };
  
    // force creation of symbols and code ...
    void init_CRshapefunctions()
    {
      InitCRShapefunctions<double, double, 3> i1;
      InitCRShapefunctions<float, double, 3> i2;
      InitCRShapefunctions<double, float, 3> i3;
      InitCRShapefunctions<float, float, 3> i4;
    }
  }

}
