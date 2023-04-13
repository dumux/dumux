// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//

/*!
 * \file
 * \ingroup PNMTwoPModel
 * \brief A (quasi-) static two-phase pore-network model for drainage processes.
 */
#ifndef DUMUX_PNM_TWOP_STATIC_DRAINAGE_HH
#define DUMUX_PNM_TWOP_STATIC_DRAINAGE_HH

#include <stack>
#include <vector>

namespace Dumux::PoreNetwork {

/*!
 * \ingroup PNMTwoPModel
 *
 * \brief A (quasi-) static two-phase pore-network model for drainage processes.
 *        This assumes that there are no pressure gradients within the phases and thus, no flow.
 */
template<class GridGeometry, class Scalar>
class TwoPStaticDrainage
{
    using GridView = typename GridGeometry::GridView;

public:

    TwoPStaticDrainage(const GridGeometry& gridGeometry,
                       const std::vector<Scalar>& pcEntry,
                       const std::vector<int>& throatLabel,
                       const int inletPoreLabel,
                       const int outletPoreLabel,
                       const bool allowDraingeOfOutlet = false)
    : gridView_(gridGeometry.gridView())
    , pcEntry_(pcEntry)
    , throatLabel_(throatLabel)
    , inletThroatLabel_(inletPoreLabel)
    , outletThroatLabel_(outletPoreLabel)
    , allowDraingeOfOutlet_(allowDraingeOfOutlet)
    {}

    /*!
     * \brief Updates the invasion state of the network for the given global capillary pressure.
     *
     * \param elementIsInvaded A vector storing the invasion state of the network.
     * \param pcGlobal The global capillary pressure to be applied.
     */
    void updateInvasionState(std::vector<bool>& elementIsInvaded, const Scalar pcGlobal)
    {
        // iterate over all elements (throats)
        for (const auto& element : elements(gridView_))
        {
            const auto eIdx = gridView_.indexSet().index(element);

            // if the throat is already invaded, do nothing and continue
            if (elementIsInvaded[eIdx])
                continue;

            // try to find a seed from which to start the search process
            bool isSeed = false;

            // start at the inlet boundary
            if (throatLabel_[eIdx] == inletThroatLabel_  && pcGlobal >= pcEntry_[eIdx])
            {
                ++numThroatsInvaded_;
                isSeed = true;
            }

            // if no throat gets invaded at the inlet, search the interior domain
            if (!isSeed)
            {
                // find a throat which is not invaded yet but has an invaded neighbor (invasion-percolation)
                for (const auto& intersection : intersections(gridView_, element))
                {
                    if (intersection.neighbor())
                    {
                        const auto& neighborElement = intersection.outside();
                        const auto nIdx = gridView_.indexSet().index(neighborElement);
                        if (elementIsInvaded[nIdx] && pcGlobal >= pcEntry_[eIdx] && (allowDraingeOfOutlet_ || throatLabel_[eIdx] != outletThroatLabel_))
                        {
                            ++numThroatsInvaded_;
                            isSeed = true;
                            break;
                        }
                    }
                }
            }

            // if an invaded throat is found, continue the search for neighboring invaded throats
            if (isSeed)
            {
                elementIsInvaded[eIdx] = true;

                // use iteration instead of recursion here because the recursion can get too deep
                using Element = typename GridView::template Codim<0>::Entity;
                std::stack<Element> elementStack;
                elementStack.push(element);
                while (!elementStack.empty())
                {
                    auto e = elementStack.top();
                    elementStack.pop();

                    for (const auto& intersection : intersections(gridView_, e))
                    {
                        if (intersection.neighbor())
                        {
                            const auto& neighborElement = intersection.outside();
                            const auto nIdx = gridView_.indexSet().index(neighborElement);
                            if (!elementIsInvaded[nIdx] && pcGlobal >= pcEntry_[nIdx] && (allowDraingeOfOutlet_ || throatLabel_[nIdx] != outletThroatLabel_))
                            {
                                ++numThroatsInvaded_;
                                elementIsInvaded[nIdx] = true;
                                elementStack.push(neighborElement);
                            }
                        }
                    }
                }
            }
        }
    }

    /*!
     * \brief Returns the number of invaded throats.
     */
    std::size_t numThroatsInvaded() const
    { return numThroatsInvaded_; }

private:
    const GridView& gridView_;
    const std::vector<Scalar>& pcEntry_;
    const std::vector<int>& throatLabel_;
    const int inletThroatLabel_;
    const int outletThroatLabel_;
    const int allowDraingeOfOutlet_;
    std::size_t numThroatsInvaded_ = 0;
};

} // namespace Dumux::PoreNetwork

#endif
