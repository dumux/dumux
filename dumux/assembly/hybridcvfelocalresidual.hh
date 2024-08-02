// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Assembly
 * \ingroup CVFEDiscretization
 * \brief Calculates the element-wise residual for hybrid control-volume finite element schemes
 */
#ifndef DUMUX_HYBRID_CVFE_LOCAL_RESIDUAL_HH
#define DUMUX_HYBRID_CVFE_LOCAL_RESIDUAL_HH

#include <dumux/common/properties.hh>

#include "cvfelocalresidual.hh"

namespace Dumux {

/*!
 * \ingroup Assembly
 * \ingroup CVFEDiscretization
 * \brief The element-wise residual for control-volume finite element schemes
 * \tparam TypeTag the TypeTag
 */
template<class TypeTag>
class HybridCVFELocalResidual : public CVFELocalResidual<TypeTag>
{
    using ParentType = CVFELocalResidual<TypeTag>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementBoundaryTypes = GetPropType<TypeTag, Properties::ElementBoundaryTypes>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using GridVolumeVariables = GetPropType<TypeTag, Properties::GridVolumeVariables>;
    using ElementVolumeVariables = typename GridVolumeVariables::LocalView;
    using ElementFluxVariablesCache = typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView;

public:
    using ElementResidualVector = typename ParentType::ElementResidualVector;
    using ParentType::ParentType;

    /*!
     * \name Main interface
     * \note Methods used by the assembler to compute derivatives and residual
     */
    // \{

    /*!
     * \brief Compute the storage local residual, i.e. the deviation of the
     *        storage term from zero for instationary problems.
     *
     * \param element The DUNE Codim<0> entity for which the residual
     *                ought to be calculated
     * \param fvGeometry The finite-volume geometry of the element
     * \param prevElemVolVars The volume averaged variables for all
     *                        sub-control volumes of the element at the previous time level
     * \param curElemVolVars The volume averaged variables for all
     *                       sub-control volumes of the element at the current  time level
     */
    using ParentType::evalStorage;
    ElementResidualVector evalStorage(const Element& element,
                                      const FVElementGeometry& fvGeometry,
                                      const ElementVolumeVariables& prevElemVolVars,
                                      const ElementVolumeVariables& curElemVolVars) const
    {
        assert(!this->isStationary() && "no time loop set for storage term evaluation");

        // initialize the residual vector for all scvs in this element
        ElementResidualVector residual(fvGeometry.numScv());

        // evaluate the volume terms (storage + source terms)
        // forward to the local residual specialized for the discretization methods
        for (const auto& localDof : fvLocalDofs(fvGeometry))
            this->asImp().evalStorage(residual, this->problem(), element, fvGeometry, prevElemVolVars, curElemVolVars, localDof.scv());

        this->asImp().evalElementStorage(residual, this->problem(), element, fvGeometry, prevElemVolVars, curElemVolVars);

        return residual;
    }

    /*!
     * \brief Compute the flux and source
     *
     * \param element The DUNE Codim<0> entity for which the residual
     *                ought to be calculated
     * \param fvGeometry The finite-volume geometry of the element
     * \param curElemVolVars The volume averaged variables for all
     *                       sub-control volumes of the element at the current  time level
     * \param elemFluxVarsCache The element flux variables cache
     * \param bcTypes The element boundary types
     */
    ElementResidualVector evalFluxAndSource(const Element& element,
                                            const FVElementGeometry& fvGeometry,
                                            const ElementVolumeVariables& elemVolVars,
                                            const ElementFluxVariablesCache& elemFluxVarsCache,
                                            const ElementBoundaryTypes &bcTypes) const
    {
        // initialize the residual vector for all scvs in this element
        ElementResidualVector residual(fvGeometry.numScv());

        // evaluate the volume terms (storage + source terms)
        // forward to the local residual specialized for the discretization methods
        for (const auto& localDof : fvLocalDofs(fvGeometry))
            this->asImp().evalSource(residual, this->problem(), element, fvGeometry, elemVolVars, localDof.scv());

        // forward to the local residual specialized for the discretization methods
        for (auto&& scvf : scvfs(fvGeometry))
            this->asImp().evalFlux(residual, this->problem(), element, fvGeometry, elemVolVars, bcTypes, elemFluxVarsCache, scvf);

        // ToDo: Better name? E.g. evalElementRhsTerms
        this->asImp().evalElementFluxAndSource(residual, this->problem(), element, fvGeometry, elemVolVars, bcTypes);

        return residual;
    }

    // \}

    /*!
     * \name Model specific interface to account for hybrid dofs
     * \note The following methods are hybrid model specific implementations of the local residual
     */
    // \{

    void evalElementStorage(ElementResidualVector& residual,
                            const Problem& problem,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& prevElemVolVars,
                            const ElementVolumeVariables& curElemVolVars) const
    {
        DUNE_THROW(Dune::NotImplemented, "This hybrid model does not implement an evalElementStorage method!");
    }

    void evalElementFluxAndSource(ElementResidualVector& residual,
                                  const Problem& problem,
                                  const Element& element,
                                  const FVElementGeometry& fvGeometry,
                                  const ElementVolumeVariables& curElemVolVars,
                                  const ElementBoundaryTypes &bcTypes) const
    {
        DUNE_THROW(Dune::NotImplemented, "This hybrid model does not implement an evalElementFluxAndSource method!");
    }

    // \}
};

} // end namespace Dumux

#endif
