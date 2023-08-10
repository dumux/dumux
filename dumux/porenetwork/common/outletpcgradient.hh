// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoreNetworkModels
 * \brief This file contains functions useful for all types of pore-network models,
 *        e.g. for applying zero capillary pressure gradient at the boundaries.
 */
#ifndef DUMUX_PNM_PCGradient_HH
#define DUMUX_PNM_PCGradient_HH

#include <cmath>
#include <vector>
#include <dumux/common/parameters.hh>

namespace Dumux::PoreNetwork {

/*!
 * \ingroup PoreNetworkModels
 * \brief A helper to apply capillary pressure gradient at the outlet
 */
template<class GridVariables, class SolutionVector>
class OutletCapPressureGradient
{
    using Scalar = typename GridVariables::VolumeVariables::PrimaryVariables::value_type;
    using VolumeVariables = typename GridVariables::VolumeVariables;

public:
    OutletCapPressureGradient(const GridVariables& gridVariables,
                                 const SolutionVector& sol)
        : gridVariables_(gridVariables)
        , sol_(sol)
    {}

    template <class Element, class SubcontrolVol>
    Scalar zeroPcGradientSw(const Element& element, const SubcontrolVol& scv)
    {
        const auto& fluidMatrixInteraction = problem_().spatialParams().fluidMatrixInteraction(element, scv, 0);
        Scalar pc = neighboringPorePc_(element, scv);
        return fluidMatrixInteraction.sw(pc);
    }

private:
    template <class Element, class SubcontrolVol>
    Scalar neighboringPorePc_(const Element& element, const SubcontrolVol& scv)
    {

        auto fvGeometry = localView(problem_().gridGeometry());
        auto elemVolVars = localView(gridVariables_.curGridVolVars());

        fvGeometry.bind(element);
        elemVolVars.bind(element, fvGeometry, sol_);
        const auto scvIdx = scv.localDofIndex();

        if (scvIdx == 1)
            return elemVolVars[scvIdx - 1].capillaryPressure();
        else
            return elemVolVars[scvIdx + 1].capillaryPressure();
    }

    const auto &problem_()
    {
        return gridVariables_.curGridVolVars().problem();
    }

    const GridVariables &gridVariables_;
    const SolutionVector &sol_;
};

} // end Dumux::PoreNetwork

#endif
