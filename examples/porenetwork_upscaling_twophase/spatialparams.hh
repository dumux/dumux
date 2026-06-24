// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoreNetworkModels
 * \ingroup SpatialParameters
 * \brief The default class for spatial parameters for single-phase pore-network models.
 */
#ifndef DUMUX_PNM_ONEP_PERMEABILITY_UPSCALING_SPATIAL_PARAMS_1P_HH
#define DUMUX_PNM_ONEP_PERMEABILITY_UPSCALING_SPATIAL_PARAMS_1P_HH

#include <dumux/porenetwork/common/spatialparams.hh>

namespace Dumux::PoreNetwork {

/*!
 * \ingroup PoreNetworkModels
 * \ingroup SpatialParameters
 * \brief Spatial parameters for the upscaling example
 */
template <class GridGeometry, class Scalar>
class UpscalingSpatialParams : public SpatialParams<GridGeometry, Scalar,
                                                    UpscalingSpatialParams<GridGeometry, Scalar>>
{
    using ParentType = SpatialParams<GridGeometry, Scalar,
                                     UpscalingSpatialParams<GridGeometry, Scalar>>;
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
public:

    UpscalingSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        std::vector<Scalar> totalAreaSharedWithAdjacentThroats(gridGeometry->numDofs());
        poreShapeFactor_.resize(gridGeometry->numDofs());

        for (const auto& element : elements(gridGeometry->gridView()))
        {
            const auto eIdx = gridGeometry->elementMapper().index(element);
            for (int i = 0; i < 2; ++i)
            {
                const auto vIdx = gridGeometry->gridView().indexSet().subIndex(element, i, GridView::dimension);
                totalAreaSharedWithAdjacentThroats[vIdx] += gridGeometry->throatCrossSectionalArea(eIdx);
                poreShapeFactor_[vIdx] += gridGeometry->throatShapeFactor(eIdx) * gridGeometry->throatCrossSectionalArea(eIdx);
            }
        }

        for (int i = 0; i < totalAreaSharedWithAdjacentThroats.size(); ++i)
            poreShapeFactor_[i] /= totalAreaSharedWithAdjacentThroats[i];
    }

    Scalar temperatureAtPos(const GlobalPosition& globalPos) const
    { return 283.15; }

    template<class ElementSolutionVector>
    Scalar poreLength(const Element& element,
                      const SubControlVolume& scv,
                      const ElementSolutionVector& elemSol) const
    {
        // we assume the pore length to be equal to the inscribed pore radius
        return this->asImp_().poreInscribedRadius(element, scv, elemSol);
    }

    template<class ElementSolutionVector>
    Scalar poreShapeFactor(const Element& element,
                           const SubControlVolume& scv,
                           const ElementSolutionVector& elemSol) const
    {
        // we assume the pore length to be equal to the inscribed pore radius
        return poreShapeFactor_[scv.dofIndex()];
    }

    template<class ElementSolutionVector>
    Scalar poreCrossSectionalArea(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolutionVector& elemSol) const
    {
        const Scalar r = this->asImp_().poreInscribedRadius(element, scv, elemSol);

        // we assume the pore cross sectional area to be equal to the area of
        // circle with the pore's inscribed radius
        return M_PI * r * r;
    }

    // dimensionless kinetic-energy coefficient which for non-creeping flow
    template<class ElementSolutionVector>
    Scalar kineticEnergyCoefficient(const Element& element,
                 const SubControlVolume& scv,
                 const ElementSolutionVector& elemSol) const
    { return 1.0; }

    // dimensionless momentum coefficient which for non-creeping flow
    template<class ElementSolutionVector>
    Scalar momentumCoefficient(const Element& element,
              const SubControlVolume& scv,
              const ElementSolutionVector& elemSol) const
    { return 1.0; }

private:
    std::vector<Scalar> poreShapeFactor_;
};

} // namespace Dumux::PoreNetwork

#endif
