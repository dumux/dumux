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

/*!
 * \file
 * \brief Forward declaration of properties required for the coupled Newton convergence writer
 */
namespace Dumux
{
//! \cond INTERNAL
/*!
 * \brief Writes the intermediate solutions during
 *        the Newton scheme
 */
template <class TypeTag>
struct MultiDomainConvergenceWriter
{
    typedef typename GET_PROP_TYPE(TypeTag, NewtonController) NewtonController;

    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;

    typedef typename GET_PROP_TYPE(TypeTag, SubDomain1TypeTag) SubTypeTag1;
    typedef typename GET_PROP_TYPE(TypeTag, SubDomain2TypeTag) SubTypeTag2;

    typedef typename GET_PROP_TYPE(SubTypeTag1, GridView) GridView1;
    typedef typename GET_PROP_TYPE(SubTypeTag2, GridView) GridView2;

    typedef typename GET_PROP_TYPE(SubTypeTag1, SolutionVector) SolutionVector1;
    typedef typename GET_PROP_TYPE(SubTypeTag2, SolutionVector) SolutionVector2;

    typedef Dumux::VtkMultiWriter<GridView1> VtkMultiWriter1;
    typedef Dumux::VtkMultiWriter<GridView2> VtkMultiWriter2;

    MultiDomainConvergenceWriter(NewtonController &ctl)
        : ctl_(ctl)
    {
        timeStepIndex_ = 0;
        iteration_ = 0;
        vtkMultiWriter1_ = 0;
        vtkMultiWriter2_ = 0;
    }

    ~MultiDomainConvergenceWriter()
    {
        delete vtkMultiWriter1_;
        delete vtkMultiWriter2_;
    };

    void beginTimestep()
    {
        ++timeStepIndex_;
        iteration_ = 0;
        if (!vtkMultiWriter1_)
            vtkMultiWriter1_ = new VtkMultiWriter1(problem_().subProblem1().gridView(), "convergence1");

        if (!vtkMultiWriter2_)
            vtkMultiWriter2_ = new VtkMultiWriter2(problem_().subProblem2().gridView(), "convergence2");
    };

    void beginIteration(const GridView1 &gridView1,
                        const GridView2 &gridView2)
    {
        ++ iteration_;
        vtkMultiWriter1_->beginWrite(timeStepIndex_ + iteration_ / 100.0);
        vtkMultiWriter2_->beginWrite(timeStepIndex_ + iteration_ / 100.0);
    };

    void writeFields(const SolutionVector &uLastIter,
                     const SolutionVector &deltaU)
    {
            SolutionVector1 uLastIter1;
            SolutionVector2 uLastIter2;
            SolutionVector1 deltaU1;
            SolutionVector2 deltaU2;

            uLastIter1.resize(ctl_.method().model().subModel1().numDofs());
            uLastIter2.resize(ctl_.method().model().subModel2().numDofs());
            deltaU1.resize(ctl_.method().model().subModel1().numDofs());
            deltaU2.resize(ctl_.method().model().subModel2().numDofs());

            typedef Dumux::SplitAndMerge<TypeTag> Common;

            Common::splitSolVector(uLastIter, uLastIter1, uLastIter2);
            Common::splitSolVector(deltaU, deltaU1, deltaU2);


            std::cout << "\n writing convergence file of current Newton iteration \n";
            ctl_.method().model().subModel1().addConvergenceVtkFields(*vtkMultiWriter1_, uLastIter1, deltaU1);
            ctl_.method().model().subModel2().addConvergenceVtkFields(*vtkMultiWriter2_, uLastIter2, deltaU2);
    };

    void endIteration()
    {
        vtkMultiWriter1_->endWrite();
        vtkMultiWriter2_->endWrite();
    };

    void endTimestep()
    {
        ++timeStepIndex_;
        iteration_ = 0;
    };

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


#endif
