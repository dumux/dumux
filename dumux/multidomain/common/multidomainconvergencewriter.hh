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
 * \brief Reference implementation of a newton convergence writer for coupled problems.
*/
#ifndef DUMUX_MULTIDOMAIN_CONVERGENCEWRITER_HH
#define DUMUX_MULTIDOMAIN_CONVERGENCEWRITER_HH

#include <dune/pdelab/backend/istlsolverbackend.hh>
#include <dumux/io/vtkmultiwriter.hh>
#include <dune/grid/multidomaingrid.hh>

#include "splitandmerge.hh"
#include "multidomainnewtoncontroller.hh"

namespace Dumux
{
/*!
 * \ingroup MultidomainModel
 * \brief Writes the intermediate solutions during
 *        the Newton scheme
 */
template <class TypeTag>
struct MultiDomainConvergenceWriter
{
    typedef typename GET_PROP_TYPE(TypeTag, NewtonController) NewtonController;

    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, SplitAndMerge) SplitAndMerge;

    typedef typename GET_PROP_TYPE(TypeTag, SubDomain1TypeTag) SubDomain1TypeTag;
    typedef typename GET_PROP_TYPE(TypeTag, SubDomain2TypeTag) SubDomain2TypeTag;

    typedef typename GET_PROP_TYPE(SubDomain1TypeTag, GridView) GridView1;
    typedef typename GET_PROP_TYPE(SubDomain2TypeTag, GridView) GridView2;

    typedef typename GET_PROP_TYPE(SubDomain1TypeTag, SolutionVector) SolutionVector1;
    typedef typename GET_PROP_TYPE(SubDomain2TypeTag, SolutionVector) SolutionVector2;

    typedef Dumux::VtkMultiWriter<GridView1> VtkMultiWriter1;
    typedef Dumux::VtkMultiWriter<GridView2> VtkMultiWriter2;

    /*!
    * \brief The constructor
    * \param ctl The newton controller
    */
    MultiDomainConvergenceWriter(NewtonController &ctl)
        : ctl_(ctl)
    {
        timeStepIndex_ = 0;
        iteration_ = 0;
        vtkMultiWriter1_ = 0;
        vtkMultiWriter2_ = 0;
    }

    //! \brief Destructor
    ~MultiDomainConvergenceWriter()
    {
        delete vtkMultiWriter1_;
        delete vtkMultiWriter2_;
    };

    /*!
     * \brief Start and advance in time
     */
    void beginTimestep()
    {
        ++timeStepIndex_;
        iteration_ = 0;
        if (!vtkMultiWriter1_)
            vtkMultiWriter1_ = new VtkMultiWriter1(problem_().sdProblem1().gridView(), "convergence1");

        if (!vtkMultiWriter2_)
            vtkMultiWriter2_ = new VtkMultiWriter2(problem_().sdProblem2().gridView(), "convergence2");
    }

    /*!
     * \brief Start and advance one iteration
     *
     * \param gridView1 The grid view of sub problem 1
     * \param gridView2 The grid view of sub problem 2
     */
    void beginIteration(const GridView1 &gridView1,
                        const GridView2 &gridView2)
    {
        ++ iteration_;
        vtkMultiWriter1_->beginWrite(timeStepIndex_ + iteration_ / 100.0);
        vtkMultiWriter2_->beginWrite(timeStepIndex_ + iteration_ / 100.0);
    }

    /*!
     * \brief Start and advance one iteration
     *
     * \param uLastIter The solution of the last iteration
     * \param deltaU The delta as calculated from solving the linear
     *               system of equations. This parameter also stores
     *               the updated solution.
     */
    void writeFields(const SolutionVector &uLastIter,
                     const SolutionVector &deltaU)
    {
            SolutionVector1 uLastIter1(ctl_.method().model().sdModel1().curSol());
            SolutionVector2 uLastIter2(ctl_.method().model().sdModel2().curSol());
            SolutionVector1 deltaU1(uLastIter1);
            SolutionVector2 deltaU2(uLastIter2);

            SplitAndMerge::splitSolVector(uLastIter, uLastIter1, uLastIter2);
            SplitAndMerge::splitSolVector(deltaU, deltaU1, deltaU2);

            std::cout << "\n writing convergence file of current Newton iteration \n";
            ctl_.method().model().sdModel1().addConvergenceVtkFields(*vtkMultiWriter1_, uLastIter1, deltaU1);
            ctl_.method().model().sdModel2().addConvergenceVtkFields(*vtkMultiWriter2_, uLastIter2, deltaU2);
    }

    //! \brief End of iteration
    void endIteration()
    {
        vtkMultiWriter1_->endWrite();
        vtkMultiWriter2_->endWrite();
    }

    //! \brief End of time step
    void endTimestep()
    {
        ++timeStepIndex_;
        iteration_ = 0;
    }

private:
    const Problem &problem_() const
    { return ctl_.method().problem(); }

    int timeStepIndex_;
    int iteration_;
    VtkMultiWriter1 *vtkMultiWriter1_;
    VtkMultiWriter2 *vtkMultiWriter2_;
    NewtonController &ctl_;
};

} // namespace Dumux

#endif // DUMUX_MULTIDOMAIN_CONVERGENCEWRITER_HH
