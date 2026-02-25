// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Assembly
 * \ingroup CVFEDiscretization
 * \brief Calculates the element-wise residual for control-volume finite element schemes
 */
#ifndef DUMUX_CVFE_LOCAL_RESIDUAL__HH
#define DUMUX_CVFE_LOCAL_RESIDUAL__HH

#include <dune/common/std/type_traits.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/assembly/localresidual.hh>

namespace Dumux::Experimental::Detail {

template<class Imp>
using SCVFIsOverlappingDetector = decltype(
    std::declval<Imp>().isOverlapping()
);

template<class Imp>
constexpr inline bool hasScvfIsOverlapping()
{ return Dune::Std::is_detected<SCVFIsOverlappingDetector, Imp>::value; }

} // end namespace Dumux::Experimental::Detail


namespace Dumux::Experimental {

/*!
 * \ingroup Assembly
 * \ingroup CVFEDiscretization
 * \brief The element-wise residual for control-volume finite element schemes
 * \tparam TypeTag the TypeTag
 */
template<class TypeTag>
class CVFELocalResidual : public Dumux::LocalResidual<TypeTag>
{
    using ParentType = Dumux::LocalResidual<TypeTag>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementBoundaryTypes = GetPropType<TypeTag, Properties::ElementBoundaryTypes>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using GridVariablesCache = typename GridVariables::GridVariablesCache;
    using ElementVariables = typename GridVariablesCache::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;

public:
    using ElementResidualVector = typename ParentType::ElementResidualVector;
    using ParentType::ParentType;

    //! evaluate flux residuals for one sub control volume face and add to residual
    void evalFlux(ElementResidualVector& residual,
                  const Problem& problem,
                  const Element& element,
                  const FVElementGeometry& fvGeometry,
                  const ElementVariables& elemVolVars,
                  const ElementBoundaryTypes& elemBcTypes,
                  const SubControlVolumeFace& scvf) const
    {
        const auto flux = this->asImp().evalFlux(problem, element, fvGeometry, elemVolVars, elemBcTypes, scvf);
        if (!scvf.boundary())
        {
            const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
            const auto& outsideScv = fvGeometry.scv(scvf.outsideScvIdx());
            residual[insideScv.localDofIndex()] += flux;

            // for control-volume finite element schemes with overlapping control volumes
            if constexpr (Detail::hasScvfIsOverlapping<SubControlVolumeFace>())
            {
                if (!scvf.isOverlapping())
                    residual[outsideScv.localDofIndex()] -= flux;
            }
            else
                residual[outsideScv.localDofIndex()] -= flux;
        }
        else
        {
            const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
            residual[insideScv.localDofIndex()] += flux;
        }
    }

    //! evaluate flux residuals for one sub control volume face
    NumEqVector evalFlux(const Problem& problem,
                         const Element& element,
                         const FVElementGeometry& fvGeometry,
                         const ElementVariables& elemVolVars,
                         const ElementBoundaryTypes& elemBcTypes,
                         const SubControlVolumeFace& scvf) const
    {
        NumEqVector flux(0.0);

        if(!scvf.boundary())
            flux += this->asImp().fluxIntegral(fvGeometry, elemVolVars, scvf);
        else
        {
            const auto& bcTypes = elemBcTypes.get(fvGeometry, scvf);

            // Treat Neumann and Robin ("solution dependent Neumann") boundary conditions.
            // For Dirichlet there is no addition to the residual here but they
            // are enforced strongly by replacing the residual entry afterwards.
            if (bcTypes.hasNeumann())
            {
                NumEqVector boundaryFluxes = problem.boundaryFluxIntegral(fvGeometry, elemVolVars, scvf);

                // only add fluxes to equations for which Neumann is set
                for (int eqIdx = 0; eqIdx < NumEqVector::dimension; ++eqIdx)
                    if (bcTypes.isNeumann(eqIdx))
                        flux[eqIdx] += boundaryFluxes[eqIdx];
            }
        }

        return flux;
    }
};

} // end namespace Dumux

#endif
