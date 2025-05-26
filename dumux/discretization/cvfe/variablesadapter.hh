// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CVFEDiscretization
 */
#ifndef DUMUX_CVFE_LOCAL_VARIABLES_ADAPTER_HH
#define DUMUX_CVFE_LOCAL_VARIABLES_ADAPTER_HH

#include <dune/common/std/type_traits.hh>

namespace Dumux::Detail::CVFE {

//! helper struct detecting if volumeVariables class has update function for scvs (old interface)
template<class Imp, class ES, class P, class E, class SCV>
using UpdateFunctionDetector = decltype(
    std::declval<Imp>().update(std::declval<ES>(), std::declval<P>(), std::declval<E>(), std::declval<SCV>())
);

//! Whenever the old interface is supported, we update related to scvs
template<class Imp, class ES, class P, class E, class SCV>
constexpr inline bool hasUpdateFunctionForScvs()
{ return Dune::Std::is_detected<UpdateFunctionDetector, Imp, ES, P, E, SCV>::value; }

/*!
 * \ingroup CVFEDiscretization
 * \brief A class for providing the new update interface of variables. This allows to still use the VolumesVariables.
 */
template <class VolumeVariables>
class VariablesAdapter : public VolumeVariables
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
        // As default we assume that for each localDof there is a corresponding scv
        // such that the update interface of VolumeVariables can still be used.
        if constexpr (hasUpdateFunctionForScvs<VolumeVariables, ElementSolution, Problem,
                      typename FVElementGeometry::Element, typename FVElementGeometry::SubControlVolume>())
            VolumeVariables::update(elemSol, problem, fvGeometry.element(), fvGeometry.scv(localDof.index()));
        else
            VolumeVariables::update(elemSol, problem, fvGeometry, localDof);
    };
};

} // end namespace Dumux

#endif
