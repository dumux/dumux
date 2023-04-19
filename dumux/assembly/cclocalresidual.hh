// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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
#include <dumux/common/numeqvector.hh>
#include <dumux/assembly/fvlocalresidual.hh>
#include <dumux/discretization/extrusion.hh>

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
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using Element = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView::template Codim<0>::Entity;
    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;
    using ElementBoundaryTypes = GetPropType<TypeTag, Properties::ElementBoundaryTypes>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using ElementFluxVariablesCache = typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using Extrusion = Extrusion_t<GridGeometry>;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

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
        const auto localScvIdx = scv.localDofIndex();
        residual[localScvIdx] += this->asImp().evalFlux(problem, element, fvGeometry, elemVolVars, elemFluxVarsCache, scvf);
    }

    //! evaluate the flux residual for a sub control volume face
    NumEqVector evalFlux(const Problem& problem,
                         const Element& element,
                         const FVElementGeometry& fvGeometry,
                         const ElementVolumeVariables& elemVolVars,
                         const ElementFluxVariablesCache& elemFluxVarsCache,
                         const SubControlVolumeFace& scvf) const
    {
        NumEqVector flux(0.0);

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
                auto neumannFluxes = problem.neumann(element, fvGeometry, elemVolVars, elemFluxVarsCache, scvf);

                // multiply neumann fluxes with the area and the extrusion factor
                const auto& scv = fvGeometry.scv(scvf.insideScvIdx());
                neumannFluxes *= Extrusion::area(fvGeometry, scvf)*elemVolVars[scv].extrusionFactor();

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
