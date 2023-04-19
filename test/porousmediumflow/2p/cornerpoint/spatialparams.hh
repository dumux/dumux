// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPTests
 * \brief The spatial params for the 2p cornerpoint test.
 */

#ifndef DUMUX_TWOP_CORNERPOINT_TEST_SPATIAL_PARAMS_HH
#define DUMUX_TWOP_CORNERPOINT_TEST_SPATIAL_PARAMS_HH

#if HAVE_OPM_GRID
#include <dune/common/version.hh>

#if DUNE_VERSION_GTE(OPM_GRID, 2022, 10)
#include <opm/input/eclipse/Deck/Deck.hpp>
#else
#include <opm/parser/eclipse/Deck/Deck.hpp>
#endif

#include <dumux/porousmediumflow/fvspatialparamsmp.hh>
#include <dumux/material/fluidmatrixinteractions/2p/vangenuchten.hh>

namespace Dumux {

/*!
 * \ingroup TwoPTests
 * \brief The spatial params for the 2p cornerpoint test.
 */
template<class GridGeometry, class Scalar>
class TwoPCornerPointTestSpatialParams
: public FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar, TwoPCornerPointTestSpatialParams<GridGeometry, Scalar>>
{
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using ThisType = TwoPCornerPointTestSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar, ThisType>;

    static constexpr int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using PcKrSwCurve = FluidMatrix::VanGenuchtenDefault<Scalar>;

    using DimWorldMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

public:
    using PermeabilityType = DimWorldMatrix;

    TwoPCornerPointTestSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry,
                                     std::shared_ptr<const Opm::Deck> deck)
    : ParentType(gridGeometry)
    , deck_(deck)
    , pcCurve_("SpatialParams")
    {
        homogeneous_ = getParam<bool>("Problem.Homogeneous");

        const auto getDeckData = [&] (const auto& name) {
#if DUNE_VERSION_GTE(OPM_GRID, 2022, 10)
            return (*deck_)[deck->index(name)[0]].getRawDoubleData();
#else
            return deck_->getKeyword(name).getRawDoubleData();
#endif
        };

        const std::vector<int>& globalCell = this->gridGeometry().gridView().grid().globalCell();

        if (deck_->hasKeyword("PORO")) {
            std::cout << "Found PORO..." << std::endl;
            std::vector<double> eclVector = getDeckData("PORO");
            porosity_.resize(globalCell.size());

            for (size_t i = 0; i < globalCell.size(); ++i) {
                if (homogeneous_)
                    porosity_[i] = 0.2;
                else
                    porosity_[i] = eclVector[globalCell[i]];
            }
        }

        if (deck_->hasKeyword("PERMX")) {
            std::cout << "Found PERMX..." << std::endl;
            std::vector<double> eclVector = getDeckData("PERMX");
            permX_.resize(globalCell.size());

            for (size_t i = 0; i < globalCell.size(); ++i) {
                // assume that values are given in mD = 9.86923e-16 m^2
                if (homogeneous_)
                    permX_[i] = 9.86923e-16;
                else
                    permX_[i] = 9.86923e-16*eclVector[globalCell[i]];
            }
        }

        if (deck_->hasKeyword("PERMZ")) {
            std::cout << "Found PERMZ..." << std::endl;
            std::vector<double> eclVector = getDeckData("PERMZ");
            permZ_.resize(globalCell.size());

            for (size_t i = 0; i < globalCell.size(); ++i) {
                // assume that values are given in mD = 9.86923e-16 m^2
                if (homogeneous_)
                    permZ_[i] = 9.86923e-16;
                else
                    permZ_[i] = 9.86923e-16*eclVector[globalCell[i]];
            }
        }
        else {
            std::cout << "PERMZ not found, set equal to PERMX..." << std::endl;
            permZ_ = permX_;
        }

    }

    /*!
     * \brief Function for defining the (intrinsic) permeability \f$[m^2]\f$.
     *
     * In this test, we use element-wise distributed permeabilities.
     *
     * \param element The current element
     * \param scv The sub-control volume inside the element.
     * \param elemSol The solution at the dofs connected to the element.
     * \return The permeability
     */
    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol) const
    {
        int eIdx = this->gridGeometry().gridView().indexSet().index(element);

        PermeabilityType K(0);
        K[0][0] = K[1][1] = permX_[eIdx];
        K[2][2] = permZ_[eIdx];

        return K;
    }

    /*!
     * \brief Returns the porosity \f$[-]\f$
     *
     * \param element The current element
     * \param scv The sub-control volume inside the element.
     * \param elemSol The solution at the dofs connected to the element.
     * \return The porosity
     */
    template<class ElementSolution>
    Scalar porosity(const Element& element,
                    const SubControlVolume& scv,
                    const ElementSolution& elemSol) const
    {
        int eIdx = this->gridGeometry().gridView().indexSet().index(element);
        return porosity_[eIdx];
    }


    /*!
     * \brief Returns the parameter object for the material law.
     *
     * \param element The current element
     * \param scv The sub-control volume inside the element.
     * \param elemSol The solution at the dofs connected to the element.
     * \return The material parameters object
     */
    template<class ElementSolution>
    auto fluidMatrixInteraction(const Element& element,
                                 const SubControlVolume& scv,
                                 const ElementSolution& elemSol) const
    {
        return makeFluidMatrixInteraction(pcCurve_);
    }

    /*!
     * \brief Function for defining which phase is to be considered as the wetting phase.
     *
     * \param globalPos The global position
     * \return The wetting phase index
     */
    template<class FluidSystem>
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    {
        return FluidSystem::phase0Idx;
    }

    Scalar permeabilityX(int eIdx)
    { return permX_[eIdx]; }

    Scalar permeabilityZ(int eIdx)
    { return permZ_[eIdx]; }

private:
    std::shared_ptr<const Opm::Deck> deck_; //!< the eclipse deck
    PcKrSwCurve pcCurve_;
    std::vector<Scalar> porosity_;
    std::vector<Scalar> permX_;
    std::vector<Scalar> permZ_;
    bool homogeneous_;
};

} // end namespace Dumux

#else
#warning "The opm-grid module is needed to use this class!"
#endif

#endif
