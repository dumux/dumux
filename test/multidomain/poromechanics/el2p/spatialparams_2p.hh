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
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>
#include <dumux/material/spatialparams/gstatrandomfield.hh>
#include <dumux/material/fluidmatrixinteractions/porositydeformation.hh>
#include <dumux/material/fluidmatrixinteractions/permeabilitykozenycarman.hh>

namespace Dumux {

/*!
 * \ingroup PoromechanicsTests
 * \brief The spatial parameters class for the two-phase sub problem in the el2p test problem.
 */
template<class FVGridGeometry, class Scalar, class CouplingManager>
class TwoPSpatialParams : public FVSpatialParams<FVGridGeometry, Scalar,
                                                 TwoPSpatialParams<FVGridGeometry, Scalar, CouplingManager>>
{
    using SubControlVolume = typename FVGridGeometry::SubControlVolume;
    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using ThisType = TwoPSpatialParams<FVGridGeometry, Scalar, CouplingManager>;
    using ParentType = FVSpatialParams<FVGridGeometry, Scalar, ThisType>;

public:
    using EffectiveLaw = RegularizedBrooksCorey<Scalar>;
    using MaterialLaw = EffToAbsLaw<EffectiveLaw>;
    using MaterialLawParams = typename MaterialLaw::Params;
    // export permeability type
    using PermeabilityType = Scalar;

    TwoPSpatialParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry,
                      std::shared_ptr<CouplingManager> couplingManagerPtr)
    : ParentType(fvGridGeometry)
    , couplingManagerPtr_(couplingManagerPtr)
    , initPermeability_(getParam<Scalar>("SpatialParams.Permeability"))
    , initPorosity_(getParam<Scalar>("SpatialParams.InitialPorosity"))
    {
        // given Van Genuchten m
        Scalar m = 0.457;
        // Brooks Corey lambda
        using std::pow;
        Scalar brooksCoreyLambda = m / (1 - m) * (1 - pow(0.5, 1/m));

        // residual saturations
        myMaterialParams_.setSwr(0.3);
        myMaterialParams_.setSnr(0.05);

        // parameters for the Brooks Corey law
        myMaterialParams_.setPe(1.99e4);
        myMaterialParams_.setLambda(brooksCoreyLambda);
    }

    //! Returns the porosity for a sub-control volume.
    template< class ElementSolution >
    Scalar porosity(const Element& element,
                    const SubControlVolume& scv,
                    const ElementSolution& elemSol) const
    {
        static constexpr auto poroMechId = CouplingManager::poroMechId;

        const auto& poroMechGridGeom = couplingManagerPtr_->template problem<poroMechId>().fvGridGeometry();
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
     * \brief Returns the parameter object for the Brooks-Corey material law.
     *
     * In this test, we use element-wise distributed material parameters.
     *
     * \param element The current element
     * \param scv The sub-control volume inside the element.
     * \param elemSol The solution at the dofs connected to the element.
     * \return The material parameters object
     */
    template<class ElementSolution>
    const MaterialLawParams& materialLawParams(const Element& element,
                                               const SubControlVolume& scv,
                                               const ElementSolution& elemSol) const
    {
        // do not use different parameters in the test with inverted wettability
        return myMaterialParams_;
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

    MaterialLawParams myMaterialParams_;
};

} // end namespace Dumux

#endif
