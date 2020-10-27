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
 * \ingroup PoromechanicsTests
 * \brief The spatial parameters class for the two-phase sub problem in the el2p test problem.
 */

#ifndef DUMUX_2P_TEST_SPATIALPARAMS_HH
#define DUMUX_2P_TEST_SPATIALPARAMS_HH

#include <dumux/discretization/elementsolution.hh>

#include <dumux/material/spatialparams/fv.hh>
#include <dumux/material/fluidmatrixinteractions/2p/brookscorey.hh>
#include <dumux/material/spatialparams/gstatrandomfield.hh>
#include <dumux/material/fluidmatrixinteractions/porositydeformation.hh>
#include <dumux/material/fluidmatrixinteractions/permeabilitykozenycarman.hh>

namespace Dumux {

/*!
 * \ingroup PoromechanicsTests
 * \brief The spatial parameters class for the two-phase sub problem in the el2p test problem.
 */
template<class GridGeometry, class Scalar, class CouplingManager>
class TwoPSpatialParams : public FVSpatialParams<GridGeometry, Scalar,
                                                 TwoPSpatialParams<GridGeometry, Scalar, CouplingManager>>
{
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using ThisType = TwoPSpatialParams<GridGeometry, Scalar, CouplingManager>;
    using ParentType = FVSpatialParams<GridGeometry, Scalar, ThisType>;

    using PcKrSwCurve = FluidMatrix::BrooksCoreyDefault<Scalar>;

public:
    // export permeability type
    using PermeabilityType = Scalar;

    TwoPSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry,
                      std::shared_ptr<CouplingManager> couplingManagerPtr)
    : ParentType(gridGeometry)
    , couplingManagerPtr_(couplingManagerPtr)
    , initPermeability_(getParam<Scalar>("SpatialParams.Permeability"))
    , initPorosity_(getParam<Scalar>("SpatialParams.InitialPorosity"))
    {
        // given Van Genuchten m
        const Scalar m = 0.457;
        // Brooks Corey lambda
        using std::pow;
        const Scalar brooksCoreyLambda = m / (1 - m) * (1 - pow(0.5, 1/m));

        auto baseParams = PcKrSwCurve::makeBasicParams("SpatialParams");
        baseParams.setLambda(brooksCoreyLambda);
        const auto effToAbsParams = PcKrSwCurve::makeEffToAbsParams("SpatialParams");

        pcKrSwCurve_ = std::make_unique<PcKrSwCurve>(baseParams, effToAbsParams);
    }

    //! Returns the porosity for a sub-control volume.
    template< class ElementSolution >
    Scalar porosity(const Element& element,
                    const SubControlVolume& scv,
                    const ElementSolution& elemSol) const
    {
        static constexpr auto poroMechId = CouplingManager::poroMechId;

        const auto& poroMechGridGeom = couplingManagerPtr_->problem(poroMechId).gridGeometry();
        const auto poroMechElemSol = elementSolution(element, couplingManagerPtr_->curSol()[poroMechId], poroMechGridGeom);

        // evaluate the deformation-dependent porosity at the scv center
        return PorosityDeformation<Scalar>::evaluatePorosity(poroMechGridGeom, element, scv.center(), poroMechElemSol, initPorosity_);
    }

    //! Functions for defining the (intrinsic) permeability \f$[m^2]\f$.
     template< class ElementSolution >
     PermeabilityType permeability(const Element& element,
                                   const SubControlVolume& scv,
                                   const ElementSolution& elemSol) const
     {
         PermeabilityKozenyCarman<PermeabilityType> permLaw;
         return permLaw.evaluatePermeability(initPermeability_, initPorosity_, porosity(element, scv, elemSol));
     }

     /*!
      * \brief Returns the parameters for the material law at a given location
      *
      * \param globalPos The global coordinates for the given location
      */
     auto fluidMatrixInteractionAtPos(const GlobalPosition& globalPos) const
     {
         return makeFluidMatrixInteraction(*pcKrSwCurve_);
     }

    /*!
     * \brief Function for defining which phase is to be considered as the wetting phase.
     *
     * \return the wetting phase index
     * \param globalPos The global position
     */
    template<class FluidSystem>
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    {
        return FluidSystem::phase0Idx;
    }

    //! Returns reference to the coupling manager.
    const CouplingManager& couplingManager() const
    { return *couplingManagerPtr_; }

private:
    std::shared_ptr<const CouplingManager> couplingManagerPtr_;
    Scalar initPermeability_;
    Scalar initPorosity_;
    std::unique_ptr<const PcKrSwCurve> pcKrSwCurve_;
};

} // end namespace Dumux

#endif
