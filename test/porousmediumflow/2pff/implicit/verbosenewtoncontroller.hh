// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/****************************************************************************
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
 * \brief Reference implementation of a controller class for the Newton solver.
 *
 * Usually this controller should be sufficient.
 */
#ifndef DUMUX_VERBOSE_NEWTON_CONTROLLER_HH
#define DUMUX_VERBOSE_NEWTON_CONTROLLER_HH

#include <dumux/nonlinear/newtoncontroller.hh>

namespace Dumux
{

/*!
 * \ingroup Newton
 * \brief A reference implementation of a Newton controller specific
 *        for the box scheme.
 *
 * If you want to specialize only some methods but are happy with the
 * defaults of the reference controller, derive your controller from
 * this class and simply overload the required methods.
 */
template <class TypeTag>
class VerboseNewtonController : public NewtonController<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, NewtonController) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;
    typedef typename GET_PROP_TYPE(TypeTag, NewtonMethod) NewtonMethod;
    typedef typename GET_PROP_TYPE(TypeTag, JacobianMatrix) JacobianMatrix;
    typedef typename GET_PROP_TYPE(TypeTag, JacobianAssembler) JacobianAssembler;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, VertexMapper) VertexMapper;

    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);

    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;

    typedef typename GET_PROP_TYPE(TypeTag, NewtonConvergenceWriter) ConvergenceWriter;

    typedef typename GET_PROP_TYPE(TypeTag, LinearSolver) LinearSolver;

    static constexpr int numEq = GET_PROP_VALUE(TypeTag, NumEq);

public:
    /*!
     * \brief Constructor
     */
    VerboseNewtonController(const Problem &problem)
        : NewtonController<TypeTag>(problem)
    {
        totalWastedIter_ = 0;
        totalSucceededIter_ = 0;
    }

    ~VerboseNewtonController()
    {
        std::cout << std::endl;
        std::cout << "##################################################" << std::endl;
        std::cout << "Total wasted Newton Iterations: " << totalWastedIter_ << std::endl;
        std::cout << "Total succeeded Newton Iterations: " << totalSucceededIter_ << std::endl;
        std::cout << "Total Newton Iterations: " << totalWastedIter_ + totalSucceededIter_ << std::endl;
    }


    /*!
     * \brief Called if the Newton method broke down.
     *
     * This method is called _after_ newtonEnd()
     */
    void newtonFail()
    {
        totalWastedIter_ += this->newtonNumSteps();
        NewtonController<TypeTag>::newtonFail();
    }

    /*!
     * \brief Called when the Newton method was successful.
     *
     * This method is called _after_ newtonEnd()
     */
    void newtonSucceed()
    {
        totalSucceededIter_ += this->newtonNumSteps();
    }


private:
unsigned int totalWastedIter_;
unsigned int totalSucceededIter_;
};
} // namespace Dumux

#endif
