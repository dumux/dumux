// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CVFEDiscretization
 */
#ifndef DUMUX_CVFE_LOCAL_VARIABLES_INTERFACE_HH
#define DUMUX_CVFE_LOCAL_VARIABLES_INTERFACE_HH

#include <dune/common/std/type_traits.hh>

namespace Dumux::Detail::CVFE {

//! helper struct detecting if volumeVariables class has new update function
template<class Imp, class ES, class P, class FVG, class DOF>
using UpdateFunctionDetector = decltype(
    std::declval<Imp>().update(std::declval<ES>(), std::declval<P>(), std::declval<FVG>(), std::declval<DOF>())
);

template<class Imp, class ES, class P, class FVG, class DOF>
constexpr inline bool hasUpdateFunction()
{ return Dune::Std::is_detected<UpdateFunctionDetector, Imp, ES, P, FVG, DOF>::value; }

/*!
 * \ingroup CVFEDiscretization
 * \brief A class for providing the new update interface of variables
 */
template <class VolumeVariables>
class VariablesInterface : public VolumeVariables
{
public:
    //! export the type used for the primary variables
    using PrimaryVariables = typename VolumeVariables::PrimaryVariables;

    //! export the indices type
    using Indices = typename VolumeVariables::Indices;

    /*!
     * \brief Update all quantities for a given control volume
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to
     *                be simulated
     * \param fvGeometry The local finite volume geometry
     * \param localDof The local dof
     */
    template<class ElementSolution, class Problem, class FVElementGeometry, class LocalDof>
    void update(const ElementSolution& elemSol,
                const Problem& problem,
                const FVElementGeometry& fvGeometry,
                const LocalDof& localDof)
    {
        if constexpr (hasUpdateFunction<VolumeVariables, ElementSolution, Problem, FVElementGeometry, LocalDof>())
            return VolumeVariables::update(elemSol, problem, fvGeometry, localDof);
        else
            return VolumeVariables::update(elemSol, problem, fvGeometry.element(), fvGeometry.scv(localDof.index()));
    };
};

} // end namespace Dumux

#endif
