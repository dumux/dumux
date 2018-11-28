// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \brief copydoc Dumux::BoxMaterialInterfaceParams
 */
#ifndef DUMUX_2P_BOX_MATERIAL_INTERFACE_PARAMS_HH
#define DUMUX_2P_BOX_MATERIAL_INTERFACE_PARAMS_HH

#include <dune/common/exceptions.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/box/elementsolution.hh>

namespace Dumux {

/*!
 * \ingroup TwoPModel
 * \brief Class that determines the material with the lowest capillary
 *        pressure (under fully water-saturated conditions) around the nodes
 *        of a grid. These parameters are then associated with the global degree
 *        of freedom. On the basis of these parameters, the saturations in the
 *        remaining sub-control volumes connected to the vertex can be reconstructed.
 */
template<class SpatialParams>
class BoxMaterialInterfaceParams
{
public:
    using MaterialLawParams = typename SpatialParams::MaterialLaw::Params;

    /*!
     * \brief Update the scv -> dofparameter map
     *
     * \param fvGridGeometry The finite volume grid geometry
     * \param spatialParams Class encapsulating the spatial parameters
     * \param x The current state of the solution vector
     */
    template<class FVGridGeometry, class SolutionVector>
    void update(const FVGridGeometry& fvGridGeometry,
                const SpatialParams& spatialParams,
                const SolutionVector& x)
    {
        using MaterialLaw = typename SpatialParams::MaterialLaw;

        // Make sure the spatial params return a const ref and no copy!
        using Elem = typename FVGridGeometry::GridView::template Codim<0>::Entity;
        using ElemSol = decltype( elementSolution(Elem(), x, fvGridGeometry) );
        using Scv = typename FVGridGeometry::SubControlVolume;
        using ReturnType = decltype(spatialParams.materialLawParams(Elem(), Scv(), ElemSol()));
        static_assert(std::is_lvalue_reference<ReturnType>::value,
                      "In order to use the box-interface solver please provide access "
                      "to the material law parameters via returning (const) references");

        // make sure this is only called for geometries of the box method!
        if (FVGridGeometry::discMethod != DiscretizationMethod::box)
            DUNE_THROW(Dune::InvalidStateException, "Determination of the interface material parameters with "
                                                    "this class only makes sense when using the box method!");

        isUpdated_ = true;
        isOnMaterialInterface_.resize(fvGridGeometry.numDofs(), false);
        dofParams_.resize(fvGridGeometry.numDofs(), nullptr);
        for (const auto& element : elements(fvGridGeometry.gridView()))
        {
            const auto elemSol = elementSolution(element, x, fvGridGeometry);

            auto fvGeometry = localView(fvGridGeometry);
            fvGeometry.bind(element);
            for (const auto& scv : scvs(fvGeometry))
            {
                const auto& params = spatialParams.materialLawParams(element, scv, elemSol);

                // if no parameters had been set, set them now
                if (dofParams_[scv.dofIndex()] == nullptr)
                    dofParams_[scv.dofIndex()] = &params;

                // otherwise only use the current ones if endPointPc (e.g. Brooks-Corey entry pressure) is lower
                else if (MaterialLaw::endPointPc( params ) < MaterialLaw::endPointPc( *(dofParams_[scv.dofIndex()]) ))
                {
                    dofParams_[scv.dofIndex()] = &params;
                    isOnMaterialInterface_[scv.dofIndex()] = true;
                }

                // keep track of material interfaces in any case
                else if ( !(params == *(dofParams_[scv.dofIndex()])) )
                    isOnMaterialInterface_[scv.dofIndex()] = true;
            }
        }
    }

    //! Return if this scv is connected to a material interface
    template<class Scv>
    bool isOnMaterialInterface(const Scv& scv) const
    { assert(isUpdated_); return isOnMaterialInterface_[scv.dofIndex()]; }

    //! Return the material parameters associated with the dof
    template<class Scv>
    const MaterialLawParams& getDofParams(const Scv& scv) const
    { assert(isUpdated_); return *(dofParams_[scv.dofIndex()]); }

private:
    bool isUpdated_{false};
    std::vector<bool> isOnMaterialInterface_;
    std::vector<const MaterialLawParams*> dofParams_;
};

}

#endif
