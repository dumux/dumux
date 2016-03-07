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
 *
 * \brief This class provides the infrastructure to write the
 *        convergence behaviour of the newton method into a VTK file.
 */
#ifndef DUMUX_NEWTON_CONVERGENCE_WRITER_HH
#define DUMUX_NEWTON_CONVERGENCE_WRITER_HH

#include <dumux/common/basicproperties.hh>
#include "newtoncontroller.hh"

namespace Dumux
{

namespace Properties
{
NEW_PROP_TAG(NewtonController);
NEW_PROP_TAG(SolutionVector);
}

/*!
 * \ingroup Newton
 * \brief Writes the intermediate solutions during
 *        the Newton scheme
 */
template <class TypeTag>
class NewtonConvergenceWriter
{
public:
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, NewtonController) NewtonController;

    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;

    typedef Dumux::VtkMultiWriter<GridView>  VtkMultiWriter;

    NewtonConvergenceWriter(NewtonController &ctl)
    : ctl_(ctl)
    {
        timeStepIndex_ = 0;
        iteration_ = 0;
    }

    void beginTimestep()
    {
        ++timeStepIndex_;
        iteration_ = 0;
    }

    void beginIteration(const GridView &gridView)
    {
        ++ iteration_;
        if (!vtkMultiWriter_)
            vtkMultiWriter_ = std::make_shared<VtkMultiWriter>(gridView, "convergence");
        vtkMultiWriter_->beginWrite(timeStepIndex_ + iteration_ / 100.0);
    }

    void writeFields(const SolutionVector &uLastIter,
                     const SolutionVector &deltaU)
    {
        ctl_.method().model().addConvergenceVtkFields(*vtkMultiWriter_, uLastIter, deltaU);
    }

    void endIteration()
    { vtkMultiWriter_->endWrite(); }

    void endTimestep()
    {
        iteration_ = 0;
    }

private:
    int timeStepIndex_;
    int iteration_;
    std::shared_ptr<VtkMultiWriter> vtkMultiWriter_;
    NewtonController &ctl_;
};

} // namespace Dumux

#endif
