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
 *        volume and mass fluxes of fluid phases over a face of a finite volume by means
 *        of the Darcy approximation. Specializations are provided for the different discretization methods.
 */
#ifndef DUMUX_DISCRETIZATION_STAGGERED_DARCYS_LAW_HH
#define DUMUX_DISCRETIZATION_STAGGERED_DARCYS_LAW_HH

#include <memory>

#include <dune/common/float_cmp.hh>

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>

#include <dumux/implicit/properties.hh>
#include <dumux/discretization/methods.hh>

namespace Dumux
{

namespace Properties
{
// forward declaration of properties
NEW_PROP_TAG(ProblemEnableGravity);
}

/*!
 * \ingroup DarcysLaw
 * \brief Specialization of Darcy's Law for the CCTpfa method. TODO: dummy, remove this class
 */
template <class TypeTag>
class DarcysLawImplementation<TypeTag, DiscretizationMethods::Staggered>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using FluxVarsCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);

    using Element = typename GridView::template Codim<0>::Entity;
    using IndexType = typename GridView::IndexSet::IndexType;
    using Stencil = std::vector<IndexType>;

    enum { dim = GridView::dimension} ;
    enum { dimWorld = GridView::dimensionworld} ;

    using DimWorldMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:

    static Scalar flux(const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolumeFace& scvFace,
                       int phaseIdx,
                       const FluxVarsCache& fluxVarsCache)
    {
//         const auto& tij = fluxVarsCache.tij();
//
//         // Get the inside and outside volume variables
//         const auto& insideScv = fvGeometry.scv(scvFace.insideScvIdx());
//         const auto& insideVolVars = elemVolVars[insideScv];
//         const auto& outsideVolVars = elemVolVars[scvFace.outsideScvIdx()];
//
//         auto hInside = insideVolVars.pressure(phaseIdx);
//         auto hOutside = outsideVolVars.pressure(phaseIdx);
//
//         if (GET_PARAM_FROM_GROUP(TypeTag, bool, Problem, EnableGravity))
//         {
//             // do averaging for the density
//             const auto rhoInside = insideVolVars.density(phaseIdx);
//             const auto rhoOutide = outsideVolVars.density(phaseIdx);
//             const auto rho = (rhoInside + rhoOutide)*0.5;
//
//             // ask for the gravitational acceleration in the inside neighbor
//             const auto xInside = insideScv.center();
//             const auto gInside = problem.gravityAtPos(xInside);
//
//             hInside -= rho*(gInside*xInside);
//
//             // and the outside neighbor
//             if (scvFace.boundary())
//             {
//                 const auto xOutside = scvFace.center();
//                 const auto gOutside = problem.gravityAtPos(xOutside);
//                 hOutside -= rho*(gOutside*xOutside);
//             }
//             else
//             {
//                 const auto outsideScvIdx = scvFace.outsideScvIdx();
//                 // as we assemble fluxes from the neighbor to our element the outside index
//                 // refers to the scv of our element, so we use the scv method
//                 const auto& outsideScv = fvGeometry.scv(outsideScvIdx);
//                 const auto xOutside = outsideScv.center();
//                 const auto gOutside = problem.gravityAtPos(xOutside);
//                 hOutside -= rho*(gOutside*xOutside);
//             }
//         }

            return 1.0;

//
    }

    // The flux variables cache has to be bound to an element prior to flux calculations
    // During the binding, the transmissibilities will be computed and stored using the method below.
    static Scalar calculateTransmissibilities(const Problem& problem,
                                              const Element& element,
                                              const FVElementGeometry& fvGeometry,
                                              const ElementVolumeVariables& elemVolVars,
                                              const SubControlVolumeFace& scvFace)
    {
        return 0;
    }


};

} // end namespace

#endif
