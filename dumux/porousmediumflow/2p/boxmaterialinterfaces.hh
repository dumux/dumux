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
 * \ingroup TwoPModel
 * \copydoc Dumux::BoxMaterialInterfaces
 */

#ifndef DUMUX_2P_BOX_MATERIAL_INTERFACES_HH
#define DUMUX_2P_BOX_MATERIAL_INTERFACES_HH

#include <dune/common/exceptions.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/box/elementsolution.hh>

namespace Dumux {

/*!
 * \ingroup TwoPModel
 * \brief Class that determines the material with the lowest capillary
 *        pressure (under fully water-saturated conditions) around the nodes
 *        of a grid.
 *
 * These parameters are then associated with the global degree
 * of freedom. On the basis of these parameters, the saturations in the
 * remaining sub-control volumes connected to the vertex can be reconstructed.
 */
template<class GridGeometry, class PcKrSw>
class BoxMaterialInterfaces
{
    using SubControlVolume = typename GridGeometry::SubControlVolume;

public:
    template<class SpatialParams, class SolutionVector>
    BoxMaterialInterfaces(const GridGeometry& gridGeometry,
                          const SpatialParams& spatialParams,
                          const SolutionVector& x)
    {
        update(gridGeometry, spatialParams, x);
    }

    /*!
     * \brief Updates the scv -> dofparameter map
     *
     * \param gridGeometry The finite volume grid geometry
     * \param spatialParams Class encapsulating the spatial parameters
     * \param x The current state of the solution vector
     */
    template<class SpatialParams, class SolutionVector>
    void update(const GridGeometry& gridGeometry,
                const SpatialParams& spatialParams,
                const SolutionVector& x)
    {
        // make sure this is only called for geometries of the box method!
        if (GridGeometry::discMethod != DiscretizationMethod::box)
            DUNE_THROW(Dune::InvalidStateException, "Determination of the interface material parameters with "
                                                    "this class only makes sense when using the box method!");

        isOnMaterialInterface_.resize(gridGeometry.numDofs(), false);
        pcSwAtDof_.resize(gridGeometry.numDofs(), nullptr);
        for (const auto& element : elements(gridGeometry.gridView()))
        {
            const auto elemSol = elementSolution(element, x, gridGeometry);

            auto fvGeometry = localView(gridGeometry);
            fvGeometry.bind(element);
            for (const auto& scv : scvs(fvGeometry))
            {
                const auto fluidMatrixInteraction = spatialParams.fluidMatrixInteraction(element, scv, elemSol);
                const auto& pcKrSw = fluidMatrixInteraction.pcSwCurve();

                // assert current preconditions on requiring the spatial params to store the pckrSw curve
                static_assert(std::is_lvalue_reference<typename std::decay_t<decltype(fluidMatrixInteraction)>::PcKrSwType>::value,
                              "In order to use the box-interface solver please provide access "
                              "to the material law parameters via returning (const) references");

                // if no parameters had been set, set them now
                if (!pcSwAtDof_[scv.dofIndex()])
                    pcSwAtDof_[scv.dofIndex()] = &pcKrSw;

                // otherwise only use the current ones if endPointPc (e.g. Brooks-Corey entry pressure) is lower
                else if (pcKrSw.endPointPc() < pcSwAtDof_[scv.dofIndex()]->endPointPc())
                {
                    pcSwAtDof_[scv.dofIndex()] = &pcKrSw;
                    isOnMaterialInterface_[scv.dofIndex()] = true;
                }

                // keep track of material interfaces in any case
                else if ( !(pcKrSw == *(pcSwAtDof_[scv.dofIndex()])) )
                    isOnMaterialInterface_[scv.dofIndex()] = true;
            }
        }
    }

    //! Returns if this scv is connected to a material interface
    bool isOnMaterialInterface(const SubControlVolume& scv) const
    { return isOnMaterialInterface_[scv.dofIndex()]; }

    //! Returns the material parameters associated with the dof
    const PcKrSw& pcSwAtDof(const SubControlVolume& scv) const
    { return *(pcSwAtDof_[scv.dofIndex()]); }

private:
    std::vector<bool> isOnMaterialInterface_;
    std::vector<const PcKrSw*> pcSwAtDof_;
};

} // end namespace Dumux

#endif
