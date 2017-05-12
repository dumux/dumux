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
 * \brief Calculates the residual of models based on the box scheme element-wise.
 */
#ifndef DUMUX_CC_MPFA_LOCAL_RESIDUAL_HH
#define DUMUX_CC_MPFA_LOCAL_RESIDUAL_HH

#include <dune/istl/matrix.hh>

#include <dumux/common/valgrind.hh>
#include <dumux/implicit/localresidual.hh>
#include <dumux/implicit/cellcentered/localresidual.hh>

#include "properties.hh"

namespace Dumux
{
/*!
 * \ingroup CCMpfaModel
 * \ingroup CCMpfaLocalResidual
 * \brief Element-wise calculation of the residual for models
 *        based on the fully implicit cell-centered mpfa scheme.
 *
 * \todo Please doc me more!
 */
template<class TypeTag>
class CCMpfaLocalResidual : public CCLocalResidual<TypeTag>
{
    friend typename GET_PROP_TYPE(TypeTag, LocalJacobian);
    using ParentType = CCLocalResidual<TypeTag>;
    friend class CCLocalResidual<TypeTag>;
    friend class ImplicitLocalResidual<TypeTag>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);

    enum { numEq = GET_PROP_VALUE(TypeTag, NumEq) };

    using Element = typename GridView::template Codim<0>::Entity;
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using ElementBoundaryTypes = typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);

    static const bool useTpfaBoundary = GET_PROP_VALUE(TypeTag, UseTpfaBoundary);

public:
    // copying the local residual class is not a good idea
    CCMpfaLocalResidual(const CCMpfaLocalResidual &) = delete;

    CCMpfaLocalResidual() : ParentType() {}

protected:

    PrimaryVariables computeFlux_(const Element &element,
                                  const FVElementGeometry& fvGeometry,
                                  const ElementVolumeVariables& elemVolVars,
                                  const SubControlVolumeFace &scvf,
                                  const ElementFluxVariablesCache& elemFluxVarsCache)
    {
        if (!scvf.boundary() || !useTpfaBoundary)
            return this->asImp_().computeFlux(element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);
        else
            return PrimaryVariables(0.0);

    }

    PrimaryVariables evalFlux_(const Element &element,
                               const FVElementGeometry& fvGeometry,
                               const ElementVolumeVariables& elemVolVars,
                               const SubControlVolumeFace& scvf,
                               const ElementFluxVariablesCache& elemFluxVarsCache)
    {
        if (!scvf.boundary() || !useTpfaBoundary)
          return this->asImp_().computeFlux_(element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);
        else
        {
          auto bcTypes = this->problem().boundaryTypes(element, scvf);
          return this->asImp_().evalBoundaryFluxes_(element, fvGeometry, elemVolVars, scvf, bcTypes, elemFluxVarsCache);
        }
    }

    void evalBoundary_(const Element &element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const ElementBoundaryTypes& elemBcTypes,
                       const ElementFluxVariablesCache& elemFluxVarsCache)
    {
        if (!useTpfaBoundary)
            return;

        for (auto&& scvf : scvfs(fvGeometry))
        {
            if (scvf.boundary())
            {
                auto bcTypes = this->problem().boundaryTypes(element, scvf);
                this->residual_[0] += this->asImp_().evalBoundaryFluxes_(element, fvGeometry, elemVolVars, scvf, bcTypes, elemFluxVarsCache);
            }
        }
    }

    /*!
     * \brief Add all fluxes resulting from Neumann, outflow and pure Dirichlet
     *        boundary conditions to the local residual
     */
    PrimaryVariables evalBoundaryFluxes_(const Element &element,
                                         const FVElementGeometry& fvGeometry,
                                         const ElementVolumeVariables& elemVolVars,
                                         const SubControlVolumeFace &scvf,
                                         const BoundaryTypes& bcTypes,
                                         const ElementFluxVariablesCache& elemFluxVarsCache)
    {
        // return the boundary fluxes
        if (bcTypes.hasOnlyDirichlet())
            return this->asImp_().computeFlux(element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);
        else if (bcTypes.hasOnlyNeumann())
            return this->asImp_().evalNeumannSegment_(element, fvGeometry, elemVolVars, scvf, bcTypes);
        else if (bcTypes.hasOutflow())
            DUNE_THROW(Dune::NotImplemented, "Outflow BC for mpfa methods.");
        else
            DUNE_THROW(Dune::InvalidStateException, "Mixed BC can not be set when using an mpfa method");
    }
};

}

#endif   // DUMUX_CC_LOCAL_RESIDUAL_HH
