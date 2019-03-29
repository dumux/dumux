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
 * \brief Base class for all fully implicit Richards problems
 */
#ifndef DUMUX_ELRICHARDS_PROBLEM_HH
#define DUMUX_ELRICHARDS_PROBLEM_HH

//#include <dumux/porousmediumflow/implicit/problem.hh>
#include <dumux/geomechanics/elrichards/problem.hh>


#include "properties.hh"

namespace Dumux
{
/*!
 * \ingroup RichardsModel
 * \ingroup ImplicitBaseProblems
 * \brief Base class for all fully implicit Richards problems
 *
 * For a description of the Richards model, see RichardsModel
 */
template<class TypeTag>
class ElRichardsProblem : public ImplicitPorousMediaProblem<TypeTag>
{
    typedef ImplicitPorousMediaProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;

public:
    /*!
     * \brief The constructor.
     *
     * The overloaded class must allocate all data structures
     * required, but _must not_ do any calls to the model, the
     * jacobian assembler, etc inside the constructor.
     *
     * If the problem requires information from these, the
     * ImplicitProblem::init() method be overloaded.
     *
     * \param timeManager The TimeManager which keeps track of time
     * \param gridView The GridView used by the problem.
     */
    ElRichardsProblem(TimeManager &timeManager, const GridView &gridView)
    : ParentType(timeManager, gridView)
    {
        bool gravity_ = GET_PARAM_FROM_GROUP(TypeTag, bool, Problem, EnableGravity);
        static const bool useHead = GET_PROP_VALUE(TypeTag, UseHead);

        if ( (!gravity_) && (useHead) )
            DUNE_THROW(Dune::InvalidStateException, "Don't run the Richards model with pressure head without gravity!");
    }

    /*!
     * \name Problem parameters
     */
    // \{
    /*!
     * \brief Returns the reference pressure \f$\mathrm{[Pa]}\f$ of the non-wetting
     *        phase within a control volume.
     *
     * This method MUST be overwritten by the actual problem.
     *
     * \param element The DUNE Codim<0> enitiy which intersects with
     *                the finite volume.
     * \param fvGeometry The finite volume geometry of the element.
     * \param scvIdx     The local index of the sub control volume inside the element
     */
    Scalar referencePressure(const Element &element,
                             const FVElementGeometry fvGeometry,
                             const int scvIdx) const
    { DUNE_THROW(Dune::NotImplemented, "referencePressure() method not implemented by the actual problem"); };
    // \}
};

}

#endif
