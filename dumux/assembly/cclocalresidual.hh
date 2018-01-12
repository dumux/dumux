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
 * \ingroup Assembly
 * \ingroup CCDiscretization
 * \brief Calculates the element-wise residual for cell-centered discretization schemes
 */
#ifndef DUMUX_CC_LOCAL_RESIDUAL_HH
#define DUMUX_CC_LOCAL_RESIDUAL_HH

#include <dumux/common/reservedblockvector.hh>
#include <dumux/common/properties.hh>
#include <dumux/assembly/fvlocalresidual.hh>

namespace Dumux {

/*!
 * \ingroup Assembly
 * \ingroup CCDiscretization
 * \brief Calculates the element-wise residual for the cell-centered discretization schemes
 */
template<class TypeTag>
class CCLocalResidual : public FVLocalResidual<TypeTag>
{
    using ParentType = FVLocalResidual<TypeTag>;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Element = typename GET_PROP_TYPE(TypeTag, GridView)::template Codim<0>::Entity;
    using ResidualVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using ElementBoundaryTypes = typename GET_PROP_TYPE(TypeTag, ElementBoundaryTypes);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);

public:
    using ElementResidualVector = typename ParentType::ElementResidualVector;
    using ParentType::ParentType;

    //! evaluate the flux residual for a sub control volume face and add to residual
    void evalFlux(ElementResidualVector& residual,
                  const Problem& problem,
                  const Element& element,
                  const FVElementGeometry& fvGeometry,
                  const ElementVolumeVariables& elemVolVars,
                  const ElementBoundaryTypes& elemBcTypes,
                  const ElementFluxVariablesCache& elemFluxVarsCache,
                  const SubControlVolumeFace& scvf) const
    {
        const auto& scv = fvGeometry.scv(scvf.insideScvIdx());
        const auto localScvIdx = scv.indexInElement();
        residual[localScvIdx] += evalFlux(problem, element, fvGeometry, elemVolVars, elemFluxVarsCache, scvf);
    }

    //! evaluate the flux residual for a sub control volume face
    ResidualVector evalFlux(const Problem& problem,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const ElementFluxVariablesCache& elemFluxVarsCache,
                            const SubControlVolumeFace& scvf) const
    {
        ResidualVector flux(0.0);

        // inner faces
        if (!scvf.boundary())
        {
            flux += this->asImp().computeFlux(problem, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);
        }

        // boundary faces
        else
        {
            const auto& bcTypes = problem.boundaryTypes(element, scvf);

            // Dirichlet boundaries
            if (bcTypes.hasDirichlet() && !bcTypes.hasNeumann())
                flux += this->asImp().computeFlux(problem, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);

            // Neumann and Robin ("solution dependent Neumann") boundary conditions
            else if (bcTypes.hasNeumann() && !bcTypes.hasDirichlet())
            {
                auto neumannFluxes = problem.neumann(element, fvGeometry, elemVolVars, scvf);

                // multiply neumann fluxes with the area and the extrusion factor
                const auto& scv = fvGeometry.scv(scvf.insideScvIdx());
                neumannFluxes *= scvf.area()*elemVolVars[scv].extrusionFactor();

                flux += neumannFluxes;
            }

            else
                DUNE_THROW(Dune::NotImplemented, "Mixed boundary conditions. Use pure boundary conditions by converting Dirichlet BCs to Robin BCs");
        }

        return flux;
    }
};

} // end namespace Dumux

#endif
