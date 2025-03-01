// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Discretization
 * \brief Update helper for grid variables.
 */

 #ifndef DUMUX_VARIABLES_UPDATE_HELPER_HH
 #define DUMUX_VARIABLES_UPDATE_HELPER_HH

 #include <type_traits>
 #include <dune/common/reservedvector.hh>

 #ifndef DOXYGEN

 namespace Dumux {

 template<class GridLocalVariables>
 class VariablesUpdateHelper
 {
     using ElementLocalVariables = typename GridLocalVariables::LocalView;
     using LocalVariables = typename ElementLocalVariables::LocalVariables;

 public:
     VariablesUpdateHelper(GridLocalVariables& gridLocalVariables,
                           ElementLocalVariables& elementLocalVariables)
     : gridLocalVariables_(gridLocalVariables)
     , elementLocalVariables_(elementLocalVariables)
     { }

     template<class ElementSolution, class Problem, class FVElementGeometry, class LocalDof>
     void update(const ElementSolution& elemSol,
                 const Problem& problem,
                 const FVElementGeometry& fvGeometry,
                 const LocalDof& localDof)
     { accessor_(localDof).update(elemSol, problem, fvGeometry, localDof); }

 private:
     template<class LocalDof>
     void accessor_(const LocalDof& localDof)
     {
         if constexpr (GridLocalVariables::cachingEnabled)
             return gridLocalVariables_.localVars(localDof);
         else
             return elementLocalVariables_[localDof];
     }

     GridLocalVariables& gridLocalVariables_;
     ElementLocalVariables& elementLocalVariables_;
 };

 } // end namespace Dumux

 #endif // DOXYGEN

 #endif
