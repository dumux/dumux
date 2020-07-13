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
 * \brief The two phase spatial parameters for pore network models.
 */
#ifndef DUMUX_PNM2P_SPATIAL_PARAMS_HH
#define DUMUX_PNM2P_SPATIAL_PARAMS_HH

#include <dumux/material/fluidmatrixinteractions/porenetwork/thresholdcapillarypressures.hh>
#include "porenetworkbase.hh"
#include <dumux/porenetworkflow/common/poreproperties.hh>
#include <dumux/porenetworkflow/common/throatproperties.hh>
#include <dumux/porenetworkflow/common/geometry.hh>

namespace Dumux
{

/*!
 * \ingroup SpatialParameters
 */

/**
 * \brief The base class for spatial parameters for pore network models.
 */
template<class GridGeometry, class Scalar, class MaterialLawT>
class PNMTwoPSpatialParams
: public PNMBaseSpatialParams<GridGeometry, Scalar, PNMTwoPSpatialParams<GridGeometry, Scalar, MaterialLawT>>
{
    using ParentType = PNMBaseSpatialParams<GridGeometry, Scalar, PNMTwoPSpatialParams<GridGeometry, Scalar, MaterialLawT>>;
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Element = typename GridView::template Codim<0>::Entity;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:

    using MaterialLaw = MaterialLawT;
    using MaterialLawParams = typename MaterialLaw::Params;

    PNMTwoPSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        setParams();
    }

    void setParams()
    {
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

    template<class FS, class ElementVolumeVariables>
    int wettingPhase(const Element&, const ElementVolumeVariables& elemVolVars) const
    { return 0; }

    template<class FS, class ElementSolutionVector>
    int wettingPhase(const Element&, const SubControlVolume& scv, const ElementSolutionVector& elemSol) const
    { return 0; }

    template<class ElementVolumeVariables>
    Scalar contactAngle(const Element& element,
                        const ElementVolumeVariables& elemVolVars) const
    {
        static const Scalar theta = getParam<Scalar>("SpatialParameters.ContactAngle", 0.0);
        return theta; // overload for different contact angles
    }

    template<class ElementSolutionVector>
    Scalar contactAngle(const Element& element,
                        const SubControlVolume& scv,
                        const ElementSolutionVector& elemSol) const
    {
        static const Scalar theta = getParam<Scalar>("SpatialParameters.ContactAngle", 0.0);
        return theta; // overload for different contact angles
    }

    /*!
     * \brief Return the element (throat) specific entry capillary pressure \f$ Pa\f$
     *
     * \param element The current element
     */
    template<class ElementVolumeVariables>
    const Scalar pcEntry(const Element& element, const ElementVolumeVariables& elemVolVars) const
    {
        const auto eIdx = this->gridGeometry().elementMapper().index(element);
        // take the average of both adjacent pores TODO: is this correct?
        const Scalar surfaceTension = 0.5*(elemVolVars[0].surfaceTension() + elemVolVars[1].surfaceTension());
        return ThresholdCapillaryPressures::pcEntry(surfaceTension,
                                                    this->asImp_().contactAngle(element, elemVolVars),
                                                    this->asImp_().throatRadius(element, elemVolVars),
                                                    this->gridGeometry().throatShapeFactor(eIdx));
    }

    /*!
     * \brief Return the element (throat) specific snap-off capillary pressure \f$ Pa\f$
     */
    template<class ElementVolumeVariables>
    const Scalar pcSnapoff(const Element& element, const ElementVolumeVariables& elemVolVars) const
    {
        // take the average of both adjacent pores TODO: is this correct?
        const Scalar surfaceTension = 0.5*(elemVolVars[0].surfaceTension() + elemVolVars[1].surfaceTension());
        return ThresholdCapillaryPressures::pcSnapoff(surfaceTension,
                                                      this->asImp_().contactAngle(element, elemVolVars),
                                                      this->asImp_().throatRadius(element, elemVolVars));
    }

    /*!
     * \brief Returns the parameter object for the PNM material law
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry of the element
     * \param scvIdx The local index of the sub-control volume
     */
    template<class ElementSolutionVector>
    MaterialLawParams materialLawParams(const Element& element,
                                        const SubControlVolume& scv,
                                        const ElementSolutionVector& elemSol) const
    {
        static const Scalar surfaceTension = getParam<Scalar>("SpatialParameters.SurfaceTension", 0.0725); // TODO
        const Scalar contactAngle = this->asImp_().contactAngle(element, scv, elemSol);
        const Scalar poreRadius = this->asImp_().poreRadius(element, scv, elemSol);
        return MaterialLawParams(surfaceTension, contactAngle, poreRadius);
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

} // namespace Dumux

#endif
