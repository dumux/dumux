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
 * \brief Spatial parameters for a pore-network model with non-creeping flow.
 */
#ifndef DUMUX_PNM_NONCREEPING_SPATIAL_PARAMS_1P_HH
#define DUMUX_PNM_NONCREEPING_SPATIAL_PARAMS_1P_HH

#include <dumux/material/spatialparams/porenetwork/porenetworkbase.hh>
#include <dumux/porenetwork/common/poreproperties.hh>
#include <dumux/porenetwork/common/throatproperties.hh>

namespace Dumux::PoreNetwork {

template<class GridGeometry, class Scalar>
class NonCreepingSpatialParams : public BaseSpatialParams<GridGeometry, Scalar, NonCreepingSpatialParams<GridGeometry, Scalar>>
{
    using ParentType = BaseSpatialParams<GridGeometry, Scalar, NonCreepingSpatialParams<GridGeometry, Scalar>>;

    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;

public:
    using PermeabilityType = Scalar;
    using ParentType::ParentType;

    template<class ElementSolutionVector>
    auto poreShapeFactor(const Element& element,
                         const SubControlVolume& scv,
                         const ElementSolutionVector& elemSol) const
    {
        return Throat::shapeFactorCircle<Scalar>();
    }

    template<class ElementSolutionVector>
    Scalar poreCrossSectionalArea(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolutionVector& elemSol) const
    {
        const auto poreRadius = this->gridGeometry().poreInscribedRadius(scv.dofIndex());
        return poreRadius*poreRadius*M_PI;
    }

    template<class ElementSolutionVector>
    Scalar poreLength(const Element& element,
                      const SubControlVolume& scv,
                      const ElementSolutionVector& elemSol) const
    {
        return this->gridGeometry().poreInscribedRadius(scv.dofIndex());
    }

    // dimensionless kinetic-energy coeffiecient which for non-creeping flow is equal to 1.0
    template<class ElementSolutionVector>
    Scalar kineticEnergyCoefficient(const Element& element,
                 const SubControlVolume& scv,
                 const ElementSolutionVector& elemSol) const
    { return 1.0; }

    // dimensionless momentum coeffiecient which for non-creeping flow is equal to 1.0
    template<class ElementSolutionVector>
    Scalar momentumCoefficient(const Element& element,
              const SubControlVolume& scv,
              const ElementSolutionVector& elemSol) const
    { return 1.0; }

};
} // end namespace Dumux::PoreNetwork

#endif
