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
 * \brief This file contains the data which is required to calculate
 *        diffusive mass fluxes due to molecular diffusion with Fourier's law.
 */
#ifndef DUMUX_DISCRETIZATION_STAGGERED_FOURIERS_LAW_HH
#define DUMUX_DISCRETIZATION_STAGGERED_FOURIERS_LAW_HH

#include <dune/common/float_cmp.hh>

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>

#include <dumux/implicit/properties.hh>
#include <dumux/discretization/methods.hh>
#include <dumux/discretization/fluxvariablescaching.hh>

namespace Dumux
{

namespace Properties
{
// forward declaration of properties
NEW_PROP_TAG(CellCenterPrimaryVariables);
}

/*!
 * \ingroup StaggeredFouriersLaw
 * \brief Specialization of Fourier's Law for the staggered free flow method.
 */
template <class TypeTag>
class FouriersLawImplementation<TypeTag, DiscretizationMethods::Staggered >
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using Element = typename GridView::template Codim<0>::Entity;
    using CellCenterPrimaryVariables = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;

    enum {
        massBalanceIdx = Indices::massBalanceIdx,
        energyBalanceIdx = Indices::energyBalanceIdx,
        phaseIdx = Indices::phaseIdx
    };

public:
    // state the discretization method this implementation belongs to
    static const DiscretizationMethods myDiscretizationMethod = DiscretizationMethods::Staggered;

    //! state the type for the corresponding cache and its filler
    //! We don't cache anything for this law
    using Cache = FluxVariablesCaching::EmptyDiffusionCache;
    using CacheFiller = FluxVariablesCaching::EmptyCacheFiller<TypeTag>;

    static CellCenterPrimaryVariables diffusiveFluxForCellCenter(const Problem& problem,
                                                           const Element& element,
                                                           const FVElementGeometry& fvGeometry,
                                                           const ElementVolumeVariables& elemVolVars,
                                                           const SubControlVolumeFace &scvf)
    {
        CellCenterPrimaryVariables flux(0.0);

        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& outsideScv = fvGeometry.scv(scvf.outsideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];
        const auto& outsideVolVars = scvf.boundary() ?  insideVolVars : elemVolVars[scvf.outsideScvIdx()];

        // effective conductivity tensors
        auto insideLambda = insideVolVars.thermalConductivity();
        auto outsideLambda = outsideVolVars.thermalConductivity();

        // scale by extrusion factor
        insideLambda *= insideVolVars.extrusionFactor();
        outsideLambda *= outsideVolVars.extrusionFactor();

        // the resulting averaged conductivity tensor
        const auto lambda = harmonicMean(insideLambda, outsideLambda);

        const Scalar insideTemp = insideVolVars.temperature();
        const Scalar outsideTemp = scvf.boundary() ? problem.dirichletAtPos(scvf.center())[energyBalanceIdx]
                                                   : outsideVolVars.temperature();
        const Scalar distance = scvf.boundary() ? (insideScv.dofPosition() - scvf.ipGlobal()).two_norm()
                                                : (outsideScv.dofPosition() - scvf.ipGlobal()).two_norm();

        if(scvf.boundary())
        {
            const auto bcTypes = problem.boundaryTypesAtPos(scvf.center());
            if(bcTypes.isOutflow(energyBalanceIdx) || bcTypes.isNeumann(energyBalanceIdx))
                return flux;
        }

        flux[energyBalanceIdx] = -1.0 * (insideTemp - outsideTemp);
        flux[energyBalanceIdx] *= lambda / distance;
        flux *= scvf.area() * sign(scvf.outerNormalScalar());
        return flux;
    }


};
} // end namespace

#endif
