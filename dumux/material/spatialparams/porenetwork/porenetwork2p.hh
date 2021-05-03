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
 *
 * \ingroup SpatialParameters
 * \brief The two-phase spatial parameters for pore-network models.
 */
#ifndef DUMUX_PNM2P_SPATIAL_PARAMS_HH
#define DUMUX_PNM2P_SPATIAL_PARAMS_HH

#include <dumux/material/fluidmatrixinteractions/fluidmatrixinteraction.hh>
#include <dumux/material/fluidmatrixinteractions/porenetwork/throat/thresholdcapillarypressures.hh>

#include <dumux/porenetwork/common/poreproperties.hh>
#include <dumux/porenetwork/common/throatproperties.hh>

#include "porenetworkbase.hh"

namespace Dumux::PoreNetwork {

/*!
 * \ingroup SpatialParameters
 * \ingroup PNMTwoPModel
 */

/**
 * \brief The base class for spatial parameters for pore-network models.
 */
template<class GridGeometry, class Scalar, class LocalRules, class Implementation>
class TwoPBaseSpatialParams
: public BaseSpatialParams<GridGeometry, Scalar, TwoPBaseSpatialParams<GridGeometry, Scalar, LocalRules, Implementation>>
{
    using ParentType = BaseSpatialParams<GridGeometry, Scalar, TwoPBaseSpatialParams<GridGeometry, Scalar, LocalRules, Implementation>>;
    using GridView = typename GridGeometry::GridView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;

public:

    TwoPBaseSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
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
    int wettingPhase(const Element&, const ElementVolumeVariables& elemVolVars) const
    { return 0; }

    /*!
     * \brief The index of the wetting phase within a pore body
     */
    template<class FS, class ElementSolutionVector>
    int wettingPhase(const Element&, const SubControlVolume& scv, const ElementSolutionVector& elemSol) const
    { return 0; }

    /*!
     *\brief The contact angle within a pore throat \f$[rad]\f$.
     *\note Overload for solution-dependent values.
     *
     *  \param element The element
     *  \param elemVolVars The element  volume variables
     */
    template<class ElementVolumeVariables>
    Scalar contactAngle(const Element& element,
                        const ElementVolumeVariables& elemVolVars) const
    {
        static const Scalar theta = getParam<Scalar>("SpatialParams.ContactAngle", 0.0);
        return theta;
    }

    /*!
     *\brief The contact angle within a pore body \f$[rad]\f$.
     *\note Overload for solution-dependent values.
     *
     *  \param element The element
     *  \param scv The sub-control volume
     *  \param elemSol The element solution
     */
    template<class ElementSolutionVector>
    Scalar contactAngle(const Element& element,
                        const SubControlVolume& scv,
                        const ElementSolutionVector& elemSol) const
    {
        static const Scalar theta = getParam<Scalar>("SpatialParams.ContactAngle", 0.0);
        return theta;
    }

    /*!
     * \brief Return the element (throat) specific entry capillary pressure \f$ Pa\f$
     * \param element The current element
     * \param elemVolVars The element volume variables
     */
    template<class ElementVolumeVariables>
    const Scalar pcEntry(const Element& element, const ElementVolumeVariables& elemVolVars) const
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
    const Scalar pcSnapoff(const Element& element, const ElementVolumeVariables& elemVolVars) const
    {
        // take the average of both adjacent pores TODO: is this correct?
        const Scalar surfaceTension = 0.5*(elemVolVars[0].surfaceTension() + elemVolVars[1].surfaceTension());
        return ThresholdCapillaryPressures::pcSnapoff(surfaceTension,
                                                      this->asImp_().contactAngle(element, elemVolVars),
                                                      this->asImp_().throatInscribedRadius(element, elemVolVars));
    }

    /*!
     * \brief Returns the surface tension \f$ N/m\f$
     *
     * \param element The current element
     * \param scv The sub-control volume inside the element.
     * \param elemSol The solution at the dofs connected to the element.
     */
    template<class ElementSolution>
    Scalar surfaceTension(const Element& element,
                          const SubControlVolume& scv,
                          const ElementSolution& elemSol) const
    {
        static const Scalar gamma = getParam<Scalar>("SpatialParams.SurfaceTension", 0.0725); // default to surface tension of water/air
        return gamma;
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

// TODO docme
template<class GridGeometry, class Scalar, class MaterialLawT>
class TwoPDefaultSpatialParams : public TwoPBaseSpatialParams<GridGeometry, Scalar, MaterialLawT,
                                                              TwoPDefaultSpatialParams<GridGeometry, Scalar, MaterialLawT>>
{
    using ParentType = TwoPBaseSpatialParams<GridGeometry, Scalar, MaterialLawT,
                                             TwoPDefaultSpatialParams<GridGeometry, Scalar, MaterialLawT>>;
public:
    using ParentType::ParentType;
};

} // namespace Dumux::PoreNetwork

#endif
