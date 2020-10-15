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
 * \ingroup Assembly
 * \brief The base class for an element-wise local operator
 *        in finite element schemes.
 */
#ifndef DUMUX_FV_LOCAL_OPERATOR_HH
#define DUMUX_FV_LOCAL_OPERATOR_HH

#include <dune/common/fmatrix.hh>
#include <dune/common/exceptions.hh>
#include <dune/functions/functionspacebases/subentitydofs.hh>

#include <dune/geometry/quadraturerules.hh>
#include <dune/istl/bvector.hh>

#include <dumux/discretization/fem/ipdata.hh>
#include <dumux/discretization/fem/elementboundarytypes.hh>
#include <dumux/discretization/box/elementboundarytypes.hh>

namespace Dumux {

/*!
 * \ingroup Assembly
 * \brief The element-wise local operator for finite volume schemes.
 *        This allows for element-wise evaluation of individual terms
 *        of the equations to be solved.
 * \tparam ElementVariables TODO
 * \tparam Operators The model-specific operators
 * \todo TODO: Doc this!
 */
template<class ElementVariables, class Operators>
class FVLocalOperator
{
    // The variables required for the evaluation of the equation
    using GridVars = typename ElementVariables::GridVariables;
    using PrimaryVariables = typename GridVars::PrimaryVariables;

    // The grid geometry on which the scheme operates
    using GridGeometry = typename GridVars::GridGeometry;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    // user-input defined in the problem
    using Problem = typename GridVars::Problem;
    using NumEqVector = typename ProblemTraits<Problem>::NumEqVector;
    using BoundaryTypes = typename ProblemTraits<Problem>::BoundaryTypes;

    static constexpr int dim = GridView::dimension;
    static constexpr int numEq = NumEqVector::size();

public:
    //! export the grid variables type this residual requires a local view of
    using GridVariables = GridVars;

    //! export a grid-independent alias for compatibility with non grid-based schemes
    //! TODO: Necessary??
    using Variables = GridVars;

    //! export underlying scalar type
    using Scalar = typename PrimaryVariables::value_type;

    //! the container storing the residual on all dofs of an element
    using ElementResidualVector = Dune::BlockVector<NumEqVector>;

    //! export a grid-independent alias for compatibility with non grid-based schemes
    using Residual = ElementResidualVector;

    /*!
     * \brief The constructor
     * \note The grid geometry/grid variables local views are expected to
     *       be bound to the same element
     */
    FVLocalOperator(const Element& element,
                    const FVElementGeometry& fvGeometry,
                    const ElementVariables& elemVars)
    : element_(element)
    , fvGeometry_(fvGeometry)
    , elemVariables_(elemVars)
    // TODO: Operators class with different constructor?
    //       In particular, time loop should become obsolete
    //       Also, we could construct it with the local views
    //       and then call flux() etc with less arguments. For
    //       now, since we reuse FVLocalResidual here, we leave
    //       that class as it is.
    , operators_(&elemVars.gridVariables().problem())
    {}

    /*!
     * \name Main interface
     * \note Methods used by the assembler to compute derivatives and residual
     */
    // \{

    /*!
     * \brief Compute the terms of the local residual that do not appear in
     *        time derivatives. These are the sources and the fluxes.
     */
    ElementResidualVector evalFluxesAndSources() const
    {
        const auto& problem = elemVariables_.gridVariables().problem();
        const auto& evv = elemVariables_.elemVolVars();
        const auto& efc = elemVariables_.elemFluxVarsCache();

        // TODO: ELEMBCTYPES? SHOULD THIS BE DONE SOMEWHERE ELSE?
        BoxElementBoundaryTypes<BoundaryTypes> ebt;
        ebt.update(problem, element_, fvGeometry_);

        // TODO: The spirit of FEM local operator is such that the Operators
        //       class defines how to evaluate the terms at an integration point
        //       while this class here (local operator) assembles everything for
        //       the entire element. Thus, the equivalent here would be that the
        //       model operators define how to evaluate the fluxes/storage/sources
        //       for an individual scv/scvf, and here we put it all together...
        //       For now, we let the FVLocalResidual (operators_) do this...

        // TODO: use same return type as operators (FVLocalResidual uses ReservedVector)
        const auto values = operators_.evalFluxAndSource(element_, fvGeometry_, evv, efc, ebt);
        ElementResidualVector result(fvGeometry_.numScv());
        for (const auto& scv : scvs(fvGeometry_))
            result[scv.localDofIndex()] = values[scv.localDofIndex()];
        return result;
    }

    /*!
     * \brief Compute the storage term, i.e. the term appearing in the time derivative.
     */
    ElementResidualVector evalStorage() const
    {
        const auto& problem = elemVariables_.gridVariables().problem();
        const auto& evv = elemVariables_.elemVolVars();

        // TODO: Until now, FVLocalResidual defined storage as the entire
        //       time derivative. Now it means the term above the time derivative.
        //       We should think about the correct naming here...
        // TODO: Should compute storage already multiply with volume??
        ElementResidualVector result(fvGeometry_.numScv());
        for (const auto& scv : scvs(fvGeometry_))
        {
            auto storage = operators_.computeStorage(problem, scv, evv[scv]);
            storage *= scv.volume()*evv[scv].extrusionFactor();
            result[scv.localDofIndex()] += storage;
        }

        return result;
    }

    ElementResidualVector getEmptyResidual() const
    {
        ElementResidualVector res(fvGeometry_.numScv());
        res = 0.0;
        return res;
    }

    // \}

    /*!
     * \name Interfaces for analytic Jacobian computation
     */
    // \{

    //! \todo TODO: Add interfaces. Or, should this be here at all!?

    //\}

    // \}

protected:

    //! return reference to the underlying problem
    const Problem& problem_() const
    { return elemVariables_.gridVariables().problem(); }

private:

    const Element& element_;                 //!< pointer to the element for which the residual is computed
    const FVElementGeometry& fvGeometry_;    //!< the local view on the finite element grid geometry
    const ElementVariables& elemVariables_;  //!< the local view on the grid variables
    Operators operators_; //!< evaluates storage/flux operators of the actual equation
};

} // end namespace Dumux

#endif
