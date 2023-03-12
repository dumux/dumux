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
 * \ingroup PNMTwoPModel
 * \ingroup SpatialParameters
 * \brief The two-phase spatial parameters for pore-network models.
 */
#ifndef DUMUX_PNM_2P_SPATIAL_PARAMS_HH
#define DUMUX_PNM_2P_SPATIAL_PARAMS_HH

#include <dumux/material/fluidmatrixinteractions/fluidmatrixinteraction.hh>
#include <dumux/material/fluidmatrixinteractions/porenetwork/throat/thresholdcapillarypressures.hh>

#include <dumux/porenetwork/common/spatialparams.hh>
#include <dumux/porenetwork/common/poreproperties.hh>
#include <dumux/porenetwork/common/throatproperties.hh>

#include <dumux/material/fluidmatrixinteractions/porenetwork/pore/2p/localrulesforplatonicbody.hh>
#include <dumux/io/plotpnmmateriallaw.hh>

namespace Dumux::PoreNetwork {

/*!
 * \ingroup PNMTwoPModel
 * \ingroup SpatialParameters
 * \brief The base class for spatial parameters for pore-network models.
 */
template<class GridGeometry, class Scalar, class LocalRules, class Implementation>
class TwoPSpatialParams
: public SpatialParams<GridGeometry, Scalar, Implementation>
{
    using ParentType = SpatialParams<GridGeometry, Scalar, Implementation>;
    using GridView = typename GridGeometry::GridView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    using MaterialLaw = Dumux::PoreNetwork::FluidMatrix::TwoPLocalRulesPlatonicBodyRegularization<Scalar, Dumux::PoreNetwork::FluidMatrix::TwoPLocalRulesPlatonicBody<Dumux::PoreNetwork::Pore::Shape::cube>>;

    TwoPSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        if (!gridGeometry->useSameGeometryForAllPores() && LocalRules::supportsMultipleGeometries())
            DUNE_THROW(Dune::InvalidStateException, "Your MaterialLaw does not support multiple pore body shapes.");

        if (this->gridGeometry().useSameShapeForAllThroats())
        {
            cornerHalfAngles_.resize(1);
            const auto& shape = this->gridGeometry().throatCrossSectionShape(/*eIdx*/0);
            cornerHalfAngles_[0] = Throat::cornerHalfAngles<Scalar>(shape);
        }
        else
        {
            cornerHalfAngles_.resize(this->gridGeometry().gridView().size(0));
            for (auto&& element : elements(this->gridGeometry().gridView()))
            {
                const auto eIdx = this->gridGeometry().elementMapper().index(element);
                const auto& shape = this->gridGeometry().throatCrossSectionShape(eIdx);
                cornerHalfAngles_[eIdx] = Throat::cornerHalfAngles<Scalar>(shape);
            }
        }
    }

    /*!
     * \brief The index of the wetting phase within a pore throat
     */
    template<class FS, class ElementVolumeVariables>
    int wettingPhase(const Element& element,
                     const ElementVolumeVariables& elemVolVars) const
    { return this->asImp_().template wettingPhaseAtPos<FS>(element.geometry().center()); }

    /*!
     * \brief The index of the wetting phase within a pore body
     */
    template<class FS, class ElementSolutionVector>
    int wettingPhase(const Element& element,
                     const SubControlVolume& scv,
                     const ElementSolutionVector& elemSol) const
    { return this->asImp_().template wettingPhaseAtPos<FS>(scv.center()); }

    /*!
     * \brief Function for defining which phase is to be considered as the wetting phase.
     *
     * \return the wetting phase index
     * \param globalPos The global position
     */
    template<class FluidSystem>
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    {
        DUNE_THROW(Dune::InvalidStateException,
                   "The spatial parameters do not provide "
                   "a wettingPhaseAtPos() method.");
    }

    /*!
     * \brief The contact angle within a pore throat \f$[rad]\f$.
     * \note Overload for solution-dependent values.
     *
     *  \param element The element
     *  \param elemVolVars The element  volume variables
     */
    template<class ElementVolumeVariables>
    Scalar contactAngle(const Element& element,
                        const ElementVolumeVariables& elemVolVars) const
    { return this->asImp_().contactAngleAtPos(element.geometry().center()); }

    /*!
     * \brief The contact angle within a pore body \f$[rad]\f$.
     * \note Overload for solution-dependent values.
     *
     *  \param element The element
     *  \param scv The sub-control volume
     *  \param elemSol The element solution
     */
    template<class ElementSolutionVector>
    Scalar contactAngle(const Element& element,
                        const SubControlVolume& scv,
                        const ElementSolutionVector& elemSol) const
    { return this->asImp_().contactAngleAtPos(scv.center()); }

    //! \brief Function for defining the Contact Angle
    int contactAngleAtPos(const GlobalPosition& globalPos) const
    {
        DUNE_THROW(Dune::InvalidStateException,
                   "The spatial parameters do not provide "
                   "a contactAngleAtPos() method.");
    }

    /*!
     * \brief Returns the surface tension \f$ [N/m] \f$
     *
     * \param element The current element
     * \param scv The sub-control volume inside the element.
     * \param elemSol The solution at the dofs connected to the element.
     */
    template<class ElementSolution>
    Scalar surfaceTension(const Element& element,
                          const SubControlVolume& scv,
                          const ElementSolution& elemSol) const
    { return this->asImp_().surfaceTensionAtPos(scv.center()); }

    //! \brief Function for defining the surface Tension
    Scalar surfaceTensionAtPos(const GlobalPosition& globalPos) const
    {
        DUNE_THROW(Dune::InvalidStateException,
                   "The spatial parameters do not provide "
                   "a surfaceTensionAtPos() method.");
    }

    /*!
     * \brief Return the element (throat) specific entry capillary pressure \f$ Pa\f$
     * \param element The current element
     * \param elemVolVars The element volume variables
     */
    template<class ElementVolumeVariables>
    const Scalar pcEntry(const Element& element,
                         const ElementVolumeVariables& elemVolVars) const
    {
        const auto eIdx = this->gridGeometry().elementMapper().index(element);
        // take the average of both adjacent pores TODO: is this correct?
        const Scalar surfaceTension = 0.5*(elemVolVars[0].surfaceTension() + elemVolVars[1].surfaceTension());
        return ThresholdCapillaryPressures::pcEntry(surfaceTension,
                                                    this->asImp_().contactAngle(element, elemVolVars),
                                                    this->asImp_().throatInscribedRadius(element, elemVolVars),
                                                    this->gridGeometry().throatShapeFactor(eIdx));
    }

    /*!
     * \brief Return the element (throat) specific snap-off capillary pressure \f$ Pa\f$
     * \param element The current element
     * \param elemVolVars The element volume variables
     */
    template<class ElementVolumeVariables>
    const Scalar pcSnapoff(const Element& element,
                           const ElementVolumeVariables& elemVolVars) const
    {
        // take the average of both adjacent pores TODO: is this correct?
        const Scalar surfaceTension = 0.5*(elemVolVars[0].surfaceTension() + elemVolVars[1].surfaceTension());
        return ThresholdCapillaryPressures::pcSnapoff(surfaceTension,
                                                      this->asImp_().contactAngle(element, elemVolVars),
                                                      this->asImp_().throatInscribedRadius(element, elemVolVars));
    }

    /*!
     * \brief Returns the parameter object for the pore-local pc-Sw law
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
        const auto params = typename LocalRules::BasicParams(*this, element, scv, elemSol);
        return makeFluidMatrixInteraction(LocalRules(params, "SpatialParams"));
    }

    const Dune::ReservedVector<Scalar, 4>& cornerHalfAngles(const Element& element) const
    {
        if (this->gridGeometry().useSameShapeForAllThroats())
            return cornerHalfAngles_[0];
        else
        {
            const auto eIdx = this->gridGeometry().gridView().indexSet().index(element);
            return cornerHalfAngles_[eIdx];
        }
    }

     /*!
     * \brief Plots the pore scale capillary pressure-saturation curve
     */
    template<class VolumeVariables, class DofIndices, class Problem>
    void plotPcSw(const DofIndices& dofIndices,
                  const Problem& problem,
                  const Scalar lowerSat = 0.0,
                  const Scalar upperSat = 1.0,
                  std::string plotName = "") const
    {
        if (dofIndices.empty())
            return;

        PlotLocalRules<MaterialLaw, Scalar> plot;
        auto fvGeometry = localView(this->gridGeometry());

        DofIndices handledDofs;
        handledDofs.reserve(dofIndices.size());

        auto dofFound = [&](const auto& dofs, const auto dofIdx)
        {
            return std::find(dofs.begin(), dofs.end(), dofIdx) != dofs.end();
        };

        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            // check if we already plotted all dofs
            if (handledDofs.size() == dofIndices.size())
                return;

            fvGeometry.bindElement(element);

            for (const auto& scv : scvs(fvGeometry))
            {
                // make sure to plot each dof only once
                if (!dofFound(handledDofs, scv.dofIndex()) && dofFound(dofIndices, scv.dofIndex()))
                {
                    handledDofs.push_back(scv.dofIndex());
                    plot.template plotPcOverSw<VolumeVariables>(problem, element, fvGeometry, scv, lowerSat, upperSat, plotName);
                }
            }
        }
    }

    //  /*!
    //  * \brief Plots the throat transmissibilities
    //  */
    // template<class VolumeVariables, class ElementIndices, class Problem, class ElementFluxVarsCache>
    // void plotTransmissibilities(const ElementIndices& elementIndices,
    //                             const Problem& problem,
    //                             const ElementFluxVarsCache& elementfluxvarscache) const
    // {
    //     if (elementIndices.empty())
    //         return;

    //     PlotLocalRules<MaterialLaw, Scalar> plot;

    //     ElementIndices handledElements;
    //     handledElements.reserve(elementIndices.size());

    //     auto elementFound = [&](const auto& elems, const auto eIdx)
    //     {
    //         return std::find(elems.begin(), elems.end(), eIdx) != elems.end();
    //     };

    //     auto fvGeometry = localView(this->gridGeometry());

    //     for (const auto& element : elements(this->gridGeometry().gridView()))
    //     {
    //         // check if we already plotted all throats
    //         if (handledElements.size() == elementIndices.size())
    //             return;

    //         const auto eIdx = this->gridGeometry().elementMapper().index(element);

    //         // make sure to plot each throat only once
    //         if (!elementFound(handledElements, eIdx) && elementFound(elementIndices, eIdx))
    //         {
    //                 handledElements.push_back(eIdx);
    //                 fvGeometry.bindElement(element);
    //                 plot.template plotTransmissibilities<VolumeVariables>(problem, elementfluxvarscache, element, fvGeometry);
    //         }
    //     }
    // }

private:
    std::vector<Dune::ReservedVector<Scalar, 4>> cornerHalfAngles_;
};

/*!
 * \ingroup PNMTwoPModel
 * \ingroup SpatialParameters
 * \brief The default class for spatial parameters for two-phase pore-network models.
 */
template<class GridGeometry, class Scalar, class LocalRules>
class TwoPDefaultSpatialParams
: public TwoPSpatialParams<GridGeometry, Scalar, LocalRules, TwoPDefaultSpatialParams<GridGeometry, Scalar, LocalRules>>
{
    using ParentType = TwoPSpatialParams<GridGeometry, Scalar, LocalRules,
                                         TwoPDefaultSpatialParams<GridGeometry, Scalar, LocalRules>>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
public:
    using ParentType::ParentType;

    TwoPDefaultSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    { }

    /*!
     * \brief Function for defining which phase is to be considered as the wetting phase.
     *
     * \return the wetting phase index
     * \param globalPos The global position
     */
    template<class FluidSystem>
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    { return FluidSystem::phase0Idx; }

    /*!
     * \brief Function for defining the Contact Angle
     *
     * \return the contact angle
     * \param globalPos The global position
     */
    int contactAngleAtPos(const GlobalPosition& globalPos) const
    {
        static const Scalar theta = getParam<Scalar>("SpatialParams.ContactAngle", 0.0);
        return theta;
    }

    /*!
     * \brief Function for defining the surface Tension
     *
     * \return the surface tension
     * \param globalPos The global position
     */
    Scalar surfaceTensionAtPos(const GlobalPosition& globalPos) const
    {
        static const Scalar gamma = getParam<Scalar>("SpatialParams.SurfaceTension", 0.0725); // default to surface tension of water/air
        return gamma;
    }

};

} // namespace Dumux::PoreNetwork

#endif
