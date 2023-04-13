// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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
                                                      this->asImp_().throatInscribedRadius(element, elemVolVars),
                                                      this->gridGeometry().throatCrossSectionShape(/*eIdx*/ 0));
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
