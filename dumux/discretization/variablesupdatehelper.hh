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

 template<class GridVariablesCache>
 class VariablesUpdateHelper
 {
     using ElementVariables = typename GridVariablesCache::LocalView;
     using Variables = typename ElementVariables::Variables;

 public:
     VariablesUpdateHelper(GridVariablesCache& gridVariablesCache,
                           ElementVariables& elementVariables)
     : gridVariablesCache_(gridVariablesCache)
     , elementVariables_(elementVariables)
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
         if constexpr (GridVariablesCache::cachingEnabled)
             return gridVariablesCache_.localVars(localDof);
         else
             return elementVariables_[localDof];
     }

     GridVariablesCache& gridVariablesCache_;
     ElementVariables& elementVariables_;
 };

 } // end namespace Dumux

 #endif // DOXYGEN

 #endif
