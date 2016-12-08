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
#ifndef DUMUX_BOX_LOCAL_RESIDUAL_HH
#define DUMUX_BOX_LOCAL_RESIDUAL_HH

#include <dune/geometry/type.hh>
#include <dune/istl/matrix.hh>

#include <dumux/common/valgrind.hh>
#include <dumux/implicit/localresidual.hh>

#include "properties.hh"

namespace Dumux
{
/*!
 * \ingroup BoxModel
 * \ingroup ImplicitLocalResidual
 * \brief Element-wise calculation of the residual for models
 *        based on the fully implicit box scheme.
 *
 * \todo Please doc me more!
 */
template<class TypeTag>
class BoxLocalResidual : public ImplicitLocalResidual<TypeTag>
{
    using ParentType = ImplicitLocalResidual<TypeTag>;
    friend class ImplicitLocalResidual<TypeTag>;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);

    enum {
        numEq = GET_PROP_VALUE(TypeTag, NumEq),
        dim = GridView::dimension
    };

    using Element = typename GridView::template Codim<0>::Entity;
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using ElementBoundaryTypes = typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);

public:
    // copying the local residual class is not a good idea
    BoxLocalResidual(const BoxLocalResidual &) = delete;

    BoxLocalResidual() : ParentType() {}

protected:

    void evalFluxes_(const Element& element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& elemVolVars,
                     const ElementBoundaryTypes& bcTypes,
                     const ElementFluxVariablesCache& elemFluxVarsCache)
    {
        // calculate the mass flux over the scv faces and subtract
        for (auto&& scvf : scvfs(fvGeometry))
        {
            if (!scvf.boundary())
            {
                const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
                const auto& outsideScv = fvGeometry.scv(scvf.outsideScvIdx());

                auto flux = this->asImp_().computeFlux(element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);

                this->residual_[insideScv.index()] += flux;
                this->residual_[outsideScv.index()] -= flux;
            }
        }
    }

    void evalBoundary_(const Element &element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const ElementBoundaryTypes& bcTypes,
                       const ElementFluxVariablesCache& elemFluxVarsCache)
    {
        if (bcTypes.hasNeumann() || bcTypes.hasOutflow())
        {
            for (auto&& scvf : scvfs(fvGeometry))
            {
                if (scvf.boundary())
                {
                    auto&& scv = fvGeometry.scv(scvf.insideScvIdx());
                    this->asImp_().evalBoundaryFluxes_(element, fvGeometry, elemVolVars, scvf, scv, bcTypes[scv.index()], elemFluxVarsCache);
                }
            }
        }

        // additionally treat mixed D/N conditions in a strong sense
        if (bcTypes.hasDirichlet())
        {
            for (auto&& scv : scvs(fvGeometry))
            {
                auto scvBcTypes = bcTypes[scv.index()];
                if (scvBcTypes.hasDirichlet())
                    this->asImp_().evalDirichlet_(element, fvGeometry, elemVolVars, scv, scvBcTypes);
            }
        }
    }

    /*!
     * \brief Set the values of the Dirichlet boundary control volumes
     *        of the current element.
     */
    void evalDirichlet_(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const SubControlVolume &scv,
                        const BoundaryTypes& bcTypes)
    {
        PrimaryVariables dirichletValues = this->problem().dirichlet(element, scv);

        // set the dirichlet conditions
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
        {
            if (bcTypes.isDirichlet(eqIdx))
            {
                int pvIdx = bcTypes.eqToDirichletIndex(eqIdx);
                assert(0 <= pvIdx && pvIdx < numEq);
                Valgrind::CheckDefined(dirichletValues[pvIdx]);

                // get the primary variables
                const auto& priVars = elemVolVars[scv].priVars();

                this->residual_[scv.index()][eqIdx] = priVars[pvIdx] - dirichletValues[pvIdx];
            }
        }
    }

    /*!
     * \brief Add all fluxes resulting from Neumann and outflow boundary conditions to the local residual.
     */
    void evalBoundaryFluxes_(const Element& element,
                             const FVElementGeometry& fvGeometry,
                             const ElementVolumeVariables& elemVolVars,
                             const SubControlVolumeFace &scvf,
                             const SubControlVolume& insideScv,
                             const BoundaryTypes& bcTypes,
                             const ElementFluxVariablesCache& elemFluxVarsCache)
    {

        // evaluate the Neumann conditions at the boundary face
        if (bcTypes.hasNeumann())
            this->residual_[insideScv.index()] += this->asImp_().evalNeumannSegment_(element, fvGeometry, elemVolVars, scvf, insideScv, bcTypes);

        // TODO: evaluate the outflow conditions at the boundary face
        //if (bcTypes.hasOutflow())
        //    flux += this->asImp_().evalOutflowSegment_(&intersection, bcTypes);
    }

    /*!
     * \brief Add Neumann boundary conditions for a single scv face
     */
    PrimaryVariables evalNeumannSegment_(const Element& element,
                                         const FVElementGeometry& fvGeometry,
                                         const ElementVolumeVariables& elemVolVars,
                                         const SubControlVolumeFace &scvf,
                                         const SubControlVolume& insideScv,
                                         const BoundaryTypes &bcTypes)
    {
        // temporary vector to store the neumann boundary fluxes
        PrimaryVariables flux(0);

        auto neumannFluxes = this->problem().neumann(element, fvGeometry, elemVolVars, scvf);

        // multiply neumann fluxes with the area and the extrusion factor
        neumannFluxes *= scvf.area()*elemVolVars[insideScv].extrusionFactor();

        // add fluxes to the temporary vector
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
            if (bcTypes.isNeumann(eqIdx))
                flux[eqIdx] += neumannFluxes[eqIdx];

        return flux;
    }

    /*!
    * \brief Add outflow boundary conditions for a single sub-control
    *        volume face to the local residual.
    *
    * \param isIt   The intersection iterator of current element
    * \param scvIdx The index of the considered face of the sub-control volume
    * \param boundaryFaceIdx The index of the considered boundary face of the sub control volume
    */
    // template <class IntersectionIterator>
    // void evalOutflowSegment_(const IntersectionIterator &isIt,
    //                         const int scvIdx,
    //                         const int boundaryFaceIdx)
    // {
    //     const BoundaryTypes &bcTypes = this->bcTypes_(scvIdx);
    //     // deal with outflow boundaries
    //     if (bcTypes.hasOutflow())
    //     {
    //         //calculate outflow fluxes
    //         PrimaryVariables values(0.0);
    //         this->asImp_().computeFlux(values, boundaryFaceIdx, true);
    //         Valgrind::CheckDefined(values);

    //         for (int equationIdx = 0; equationIdx < numEq; ++equationIdx)
    //         {
    //             if (!bcTypes.isOutflow(equationIdx) )
    //                 continue;
    //             // deduce outflow
    //             this->residual_[scvIdx][equationIdx] += values[equationIdx];
    //         }
    //     }
    // }
};

}

#endif
