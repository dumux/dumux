// $Id:$
/*****************************************************************************
 *   Copyright (C) <YEARS> by <ADD_AUTHOR_HERE>                              *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/

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

/** \todo Please doc me! */

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

/** \todo Please doc me! */

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
