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
 *        e.g. for the calculation of fluxes at the boundary.
 */
#ifndef DUMUX_PNM_UTILITES_HH
#define DUMUX_PNM_UTILITES_HH

#include <cmath>
#include <vector>
#include <dumux/common/parameters.hh>

namespace Dumux::PoreNetwork {

/*!
 * \ingroup PoreNetworkModels
 * \brief Calculates averaged values of the network solution.
 */
template<class GridVariables, class SolutionVector>
class AveragedValues
{
    using Scalar = typename GridVariables::VolumeVariables::PrimaryVariables::value_type;
    using VolumeVariables = typename GridVariables::VolumeVariables;

    struct AveragedQuantityInfo
    {
        std::function<Scalar(const VolumeVariables&)> quantity;
        std::function<Scalar(const VolumeVariables&)> weight;
        std::string name;
    };

public:

    AveragedValues(const GridVariables& gridVariables,
                   const SolutionVector& sol)
    : gridVariables_(gridVariables)
    , sol_(sol)
    {}

    void addAveragedQuantity(std::function<Scalar(const VolumeVariables&)>&& q,
                             std::function<Scalar(const VolumeVariables&)>&& w,
                             const std::string& name)
    {
        averagedQuantityInfo_.push_back(AveragedQuantityInfo{q, w, name});
        averagedQuantity_.push_back(Scalar());
    }

    void eval(const std::vector<std::size_t>& dofsToNeglect = std::vector<std::size_t>())
    {
        for (auto& q : averagedQuantity_)
            q = 0.0;

        std::vector<bool> poreVisited(problem_().gridGeometry().numDofs(), false);
        std::vector<Scalar> weights(averagedQuantityInfo_.size(), 0.0);

        auto fvGeometry = localView(problem_().gridGeometry());
        auto elemVolVars = localView(gridVariables_.curGridVolVars());

        for (const auto& element : elements(problem_().gridGeometry().gridView()))
        {
            fvGeometry.bind(element);
            elemVolVars.bind(element, fvGeometry, sol_);

            for (int scvIdx = 0; scvIdx < fvGeometry.numScv(); ++scvIdx)
            {
                static constexpr auto dofCodim = std::decay_t<decltype(problem_().gridGeometry().gridView())>::dimension;
                const int dofIdxGlobal = problem_().gridGeometry().vertexMapper().subIndex(element, scvIdx, dofCodim);

                if (poreVisited[dofIdxGlobal])
                    continue;
                else if (!dofsToNeglect.empty() && std::any_of(dofsToNeglect.begin(), dofsToNeglect.end(), [&](int dofIdx){ return dofIdx == dofIdxGlobal; }))
                    continue;
                else
                {
                    const auto& volVars = elemVolVars[scvIdx];
                    for (int i = 0; i < averagedQuantityInfo_.size(); ++i)
                    {
                        const Scalar weight = averagedQuantityInfo_[i].weight(volVars);
                        averagedQuantity_[i] += averagedQuantityInfo_[i].quantity(volVars) * weight;
                        weights[i] += weight;
                    }
                    poreVisited[dofIdxGlobal] = true;
                }
            }
        }

        for (int i = 0; i < averagedQuantityInfo_.size(); ++i)
            averagedQuantity_[i] /= weights[i];
    }

    const Scalar& operator[](const std::string& name) const
    {
        for (int i = 0; i < averagedQuantityInfo_.size(); ++i)
            if (averagedQuantityInfo_[i].name == name)
                return averagedQuantity_[i];

        DUNE_THROW(Dune::InvalidStateException, name + " not found");
    }

private:

    std::vector<AveragedQuantityInfo> averagedQuantityInfo_;
    std::vector<Scalar> averagedQuantity_;

    const auto& problem_()
    { return gridVariables_.curGridVolVars().problem(); }

    const GridVariables& gridVariables_;
    const SolutionVector& sol_;
};

} // end Dumux::PoreNetwork

#endif
