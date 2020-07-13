// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \copydoc Dumux::PNMNewtonConsistencyChecks
 */
#ifndef DUMUX_PNM_NEWTON_CONSISTENCY_CHECKS
#define DUMUX_PNM_NEWTON_CONSISTENCY_CHECKS

#include <vector>
#include <iostream>
#include <dune/common/exceptions.hh>
#include <dumux/discretization/localview.hh>

namespace Dumux {

namespace Detail {

// Primary template
template <typename T, typename U = int>
struct SaturationIndex {};

// Specialization for 2p
template <typename T>
struct SaturationIndex <T, decltype((void) T::saturationIdx, 0)>
{
    static constexpr auto value = T::saturationIdx;
};

// Specialization for 2pnc
template <typename T>
struct SaturationIndex <T, decltype((void) T::switchIdx, 0)>
{
    static constexpr auto value = T::switchIdx;
};

// Specialization for MPNC
template <typename T>
struct SaturationIndex <T, decltype((void) T::s0Idx, 0)>
{
    static constexpr auto value = T::s0Idx;
};
} // end namespace Detail

/*!
 * \ingroup PNMTwoP
 * \brief Consistency checks for the PNM two-phase model
 */
template<class GridVariables, class SolutionVector>
class PNMNewtonConsistencyChecks
{
public:

    /*!
     * \brief Perform all checks.
     */
    void performChecks(const GridVariables& gridVariables, const SolutionVector& uCurrentIter, const SolutionVector& prevSol) const
    {
        checkRelativeSaturationShift(gridVariables, uCurrentIter, prevSol);
        checkIfValuesArePhysical(gridVariables, uCurrentIter);
        checkIfCapillaryPressureIsCloseToEntryPressure(gridVariables, uCurrentIter);
    }

    /*!
     * \brief Checks if the relative shift of saturation  between to consecutive time steps is below a given treshold.
     */
    void checkRelativeSaturationShift(const GridVariables& gridVariables, const SolutionVector& uCurrentIter, const SolutionVector& prevSol) const
    {
        using Scalar = typename SolutionVector::block_type::value_type;
        const auto& problem = gridVariables.curGridVolVars().problem();

        static const Scalar allowedSaturationChange = getParamFromGroup<Scalar>(problem.paramGroup(), "Newton.AllowedSaturationChange", -1.0);
        if (allowedSaturationChange < 0.0)
            return;

        const auto& curGridVolVars = gridVariables.curGridVolVars();
        const auto& prevGridVolVars = gridVariables.prevGridVolVars();

        auto fvGeometry = localView(problem.gridGeometry());
        auto curElemVolVars = localView(curGridVolVars);
        auto prevElemVolVars = localView(prevGridVolVars);

        std::vector<bool> dofVisited(uCurrentIter.size(), false);

        for (const auto& element : elements(problem.gridGeometry().gridView()))
        {
            fvGeometry.bindElement(element);
            curElemVolVars.bindElement(element, fvGeometry, uCurrentIter);
            prevElemVolVars.bindElement(element, fvGeometry, prevSol);

            for (const auto& scv : scvs(fvGeometry))
            {
                if (dofVisited[scv.dofIndex()])
                    continue;
                dofVisited[scv.dofIndex()] = true;

                const Scalar satNew = curElemVolVars[scv].saturation(0);
                const Scalar satOld = prevElemVolVars[scv].saturation(0);
                using std::abs;
                const Scalar deltaS = abs(satNew - satOld);
                static const bool considerRelativeShift = getParamFromGroup<Scalar>(problem.paramGroup(), "Newton.SaturationChangeIsRelative", false);

                if (considerRelativeShift)
                {
                    // satOld mgiht be (close to) zero, so have to take care of this
                    if (satOld > 1e-3 && deltaS / satOld > allowedSaturationChange)
                        DUNE_THROW(NumericalProblem, "Saturation change too high at dof " << scv.dofIndex() << ", old sat. " << satOld << ", new sat. " << satNew << std::endl);
                }
                else if (deltaS > allowedSaturationChange)
                    DUNE_THROW(NumericalProblem, "Saturation change too high at dof " << scv.dofIndex() << ", old sat. " << satOld << ", new sat. " << satNew << std::endl);
            }
        }
    }

    /*!
     * \brief Checks if the saturation is between zero and one.
     */
    void checkIfValuesArePhysical(const GridVariables& gridVariables, const SolutionVector& uCurrentIter) const
    {
        const auto& problem = gridVariables.curGridVolVars().problem();

        static const bool doPlausibilityCheck = getParamFromGroup<bool>(problem.paramGroup(), "Newton.PlausibilityCheck", false);
        if (!doPlausibilityCheck)
            return;

        for (std::size_t i = 0; i < uCurrentIter.size(); ++i)
        {
            const auto& priVars = uCurrentIter[i];
            using Indices = typename GridVariables::VolumeVariables::Indices;
            static constexpr auto saturationOrMoleFractionIndex = Detail::SaturationIndex<Indices>::value;
            if (priVars[saturationOrMoleFractionIndex] < 0.0 || priVars[saturationOrMoleFractionIndex] > 1.0)
            {
                std::cout << "at dof " << i << ": saturation " << priVars[saturationOrMoleFractionIndex] << std::endl;
                DUNE_THROW(NumericalProblem, "Saturation is below 0 or above 1");
            }
        }
    }

    /*!
     * \brief Checks if the capillary pressure at pore from which a throat was invaded is sufficiently close to the throat's entry capillary pressure.
     */
    void checkIfCapillaryPressureIsCloseToEntryPressure(const GridVariables& gridVariables, const SolutionVector& uCurrentIter) const
    {
        // this check is implemented in the invasion state itself
        const auto& invasionState = gridVariables.gridFluxVarsCache().invasionState();
        invasionState.checkIfCapillaryPressureIsCloseToEntryPressure(uCurrentIter, gridVariables.curGridVolVars(), gridVariables.gridFluxVarsCache());
    }
};
} // end namespace Dumux
#endif
