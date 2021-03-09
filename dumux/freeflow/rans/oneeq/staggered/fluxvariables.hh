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
 * \ingroup OneEqModel
 * \copydoc Dumux::OneEqFluxVariablesImpl
 */
#ifndef DUMUX_ONEEQ_STAGGERED_FLUXVARIABLES_HH
#define DUMUX_ONEEQ_STAGGERED_FLUXVARIABLES_HH

#include <numeric>
#include <dumux/common/properties.hh>
#include <dumux/flux/fluxvariablesbase.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/freeflow/navierstokes/fluxvariables.hh>
#include <dumux/freeflow/rans/oneeq/fluxvariables.hh>

namespace Dumux {

/*!
 * \ingroup OneEqModel
 * \brief The flux variables class for the one-equation model by Spalart-Allmaras
 *        using the staggered grid discretization.
 */

// forward declaration
template<class TypeTag, class BaseFluxVariables, DiscretizationMethod discMethod>
class OneEqFluxVariablesImpl;

template<class TypeTag, class BaseFluxVariables>
class OneEqFluxVariablesImpl<TypeTag, BaseFluxVariables, DiscretizationMethod::staggered>
: public BaseFluxVariables
{
    using ParentType = BaseFluxVariables;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;

    using GridVolumeVariables = typename GridVariables::GridVolumeVariables;
    using ElementVolumeVariables = typename GridVolumeVariables::LocalView;
    using VolumeVariables = typename GridVolumeVariables::VolumeVariables;

    using GridFluxVariablesCache = typename GridVariables::GridFluxVariablesCache;
    using FluxVariablesCache = typename GridFluxVariablesCache::FluxVariablesCache;

    using GridFaceVariables = typename GridVariables::GridFaceVariables;
    using ElementFaceVariables = typename GridFaceVariables::LocalView;
    using FaceVariables = typename GridFaceVariables::FaceVariables;

    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Extrusion = Extrusion_t<GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using Element = typename GridView::template Codim<0>::Entity;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using CellCenterPrimaryVariables = GetPropType<TypeTag, Properties::CellCenterPrimaryVariables>;

    static constexpr int viscosityTildeEqIdx = Indices::viscosityTildeEqIdx - ModelTraits::dim();

public:

    /*!
     * \brief Computes the flux for the cell center residual.
     */
    CellCenterPrimaryVariables computeMassFlux(const Problem& problem,
                                               const Element &element,
                                               const FVElementGeometry& fvGeometry,
                                               const ElementVolumeVariables& elemVolVars,
                                               const ElementFaceVariables& elemFaceVars,
                                               const SubControlVolumeFace &scvf,
                                               const FluxVariablesCache& fluxVarsCache)
    {
        CellCenterPrimaryVariables flux = ParentType::computeMassFlux(problem, element, fvGeometry,
                                                                      elemVolVars, elemFaceVars, scvf, fluxVarsCache);

        // calculate advective flux
        auto upwindTermK = [](const auto& volVars)
        {
            return volVars.viscosityTilde() * volVars.density();
        };

        flux[viscosityTildeEqIdx]
            = ParentType::advectiveFluxForCellCenter(problem, elemVolVars, elemFaceVars, scvf, upwindTermK);

        // calculate diffusive flux
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& outsideScv = fvGeometry.scv(scvf.outsideScvIdx());
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];

        // effective diffusion coefficients
        Scalar insideCoeff = (insideVolVars.viscosity() + insideVolVars.viscosityTilde() * insideVolVars.density())
                             / insideVolVars.sigma();
        Scalar outsideCoeff = (outsideVolVars.viscosity() + outsideVolVars.viscosityTilde() * outsideVolVars.density())
                              / outsideVolVars.sigma();

        // scale by extrusion factor
        insideCoeff *= insideVolVars.extrusionFactor();
        outsideCoeff *= outsideVolVars.extrusionFactor();


        Scalar coeff = 0.0;
        Scalar distance = 0.0;
        if (scvf.boundary())
        {
            coeff = insideCoeff;
            distance = (insideScv.dofPosition() - scvf.ipGlobal()).two_norm();
        }
        else
        {
            // average and distance
            coeff = arithmeticMean(insideCoeff, outsideCoeff,
                                  (outsideScv.dofPosition() - scvf.ipGlobal()).two_norm(),
                                  (insideScv.dofPosition() - scvf.ipGlobal()).two_norm());
            distance = (outsideScv.dofPosition() - insideScv.dofPosition()).two_norm();
        }

        const auto bcTypes = problem.boundaryTypes(element, scvf);
        if (!(scvf.boundary() && (bcTypes.isOutflow(Indices::viscosityTildeEqIdx)
                                  || bcTypes.isSymmetry())))
        {
            flux[viscosityTildeEqIdx]
                += coeff / distance
                   * (insideVolVars.viscosityTilde() - outsideVolVars.viscosityTilde())
                   * Extrusion::area(scvf);
        }
        return flux;
    }
};

} // end namespace

#endif
