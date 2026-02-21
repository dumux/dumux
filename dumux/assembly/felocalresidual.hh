// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Assembly
 * \ingroup FEDiscretization
 * \brief Calculates the element-wise residual for control-volume finite element schemes
 */
#ifndef DUMUX_FE_LOCAL_RESIDUAL_HH
#define DUMUX_FE_LOCAL_RESIDUAL_HH

#include <dune/common/std/type_traits.hh>
#include <dune/geometry/type.hh>
#include <dune/istl/matrix.hh>

#include <dumux/common/typetraits/localdofs_.hh>
#include <dumux/common/typetraits/boundary_.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/timeloop.hh>
#include <dumux/common/reservedblockvector.hh>

namespace Dumux {

/*!
 * \ingroup Assembly
 * \ingroup FEDiscretization
 * \brief The element-wise residual for finite element schemes
 * \tparam TypeTag the TypeTag
 */
template<class TypeTag>
class FELocalResidual
{
    using Implementation = GetPropType<TypeTag, Properties::LocalResidual>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementBoundaryTypes = GetPropType<TypeTag, Properties::ElementBoundaryTypes>;
    using FEElementGeometry = typename GridGeometry::LocalView;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    // ToDo: Replace them once we have established the new variables concept
    using GridVariablesCache = typename GridVariables::GridVolumeVariables;
    using ElementVariables = typename GridVariablesCache::LocalView;
    using ElementDataCache = typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView;
    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;
    using TimeLoop = TimeLoopBase<Scalar>;

public:
    //! the container storing all element residuals
    using ElementResidualVector = ReservedBlockVector<NumEqVector, Detail::LocalDofs::maxNumLocalDofs<FEElementGeometry>()>;

    //! the constructor
    FELocalResidual(const Problem* problem,
                    const TimeLoop* timeLoop = nullptr)
    : problem_(problem)
    , timeLoop_(timeLoop)
    {}

    /*!
     * \brief Compute the storage local residual, i.e. the deviation of the
     *        storage term from zero for instationary problems.
     *
     * \param element The DUNE Codim<0> entity for which the residual
     *                ought to be calculated
     * \param fvGeometry The finite-volume geometry of the element
     * \param prevElemVars The variables related to the element at the previous time level
     * \param curElemVars The variables related to the element at the current time level
     */
    ElementResidualVector evalStorage(const Element& element,
                                      const FEElementGeometry& fvGeometry,
                                      const ElementVariables& prevElemVars,
                                      const ElementVariables& curElemVars) const
    {
        assert(!this->isStationary() && "no time loop set for storage term evaluation");

        // initialize the residual vector for all scvs in this element
        ElementResidualVector residual(Detail::LocalDofs::numLocalDofs(fvGeometry));

        // allow for additional contributions (e.g. hybrid CVFE schemes)
        this->asImp().addToElementStorageResidual(residual, this->problem(), element, fvGeometry, prevElemVars, curElemVars);

        return residual;
    }

    /*!
     * \brief Compute the flux and source
     *
     * \param element The DUNE Codim<0> entity for which the residual
     *                ought to be calculated
     * \param fvGeometry The finite-volume geometry of the element
     * \param elemVars The variables related to the element at the current time level
     * \param elemDataCache The element data cache
     * \param bcTypes The element boundary types
     */
    ElementResidualVector evalFluxAndSource(const Element& element,
                                            const FEElementGeometry& fvGeometry,
                                            const ElementVariables& elemVars,
                                            const ElementDataCache& elemDataCache,
                                            const ElementBoundaryTypes &bcTypes) const
    {
        // initialize the residual vector for all scvs in this element
        ElementResidualVector residual(Detail::LocalDofs::numLocalDofs(fvGeometry));

        // allow for additional contributions (e.g. hybrid CVFE schemes)
        this->asImp().addToElementFluxAndSourceResidual(residual, this->problem(), element, fvGeometry, elemVars, elemDataCache, bcTypes);

        return residual;
    }

    //! add additional storage contributions (e.g. hybrid CVFE schemes)
    void addToElementStorageResidual(ElementResidualVector& residual,
                                     const Problem& problem,
                                     const Element& element,
                                     const FEElementGeometry& fvGeometry,
                                     const ElementVariables& prevElemVars,
                                     const ElementVariables& curElemVars) const
    {}

    //! add additional flux and source contributions (e.g. hybrid CVFE schemes)
    void addToElementFluxAndSourceResidual(ElementResidualVector& residual,
                                           const Problem& problem,
                                           const Element& element,
                                           const FEElementGeometry& fvGeometry,
                                           const ElementVariables& curElemVars,
                                           const ElementDataCache& elemDataCache,
                                           const ElementBoundaryTypes &bcTypes) const
    {}

    /*!
     * \name Interfaces accessed by local residual implementations
     */
    // \{

    //! the problem
    const Problem& problem() const
    { return *problem_; }

    //! the timeloop for instationary problems
    //! calling this for stationary leads to undefined behaviour
    const TimeLoop& timeLoop() const
    { return *timeLoop_; }

    //! returns true if the residual is stationary
    bool isStationary() const
    { return !timeLoop_; }

    // \}
protected:

    Implementation &asImp()
    { return *static_cast<Implementation*>(this); }

    const Implementation &asImp() const
    { return *static_cast<const Implementation*>(this); }

private:
    const Problem* problem_; //!< the problem we are assembling this residual for
    const TimeLoop* timeLoop_; //!< the timeloop for instationary problems
};

} // end namespace Dumux

#endif
