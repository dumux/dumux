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
 * \brief A VTK output module to simplify writing dumux simulation data to VTK format
 */
#ifndef STAGGERED_VTK_OUTPUT_MODULE_HH
#define STAGGGERED_VTK_OUTPUT_MODULE_HH

#include <dune/common/fvector.hh>

#include <dumux/io/vtkoutputmodulebase.hh>
#include <dumux/io/staggeredvtkwriter.hh>

namespace Properties
{
NEW_PROP_TAG(VtkAddVelocity);
NEW_PROP_TAG(VtkAddProcessRank);
}

namespace Dumux
{

/*!
 * \ingroup InputOutput
 * \brief A VTK output module to simplify writing dumux simulation data to VTK format
 *
 * Handles the output of scalar and vector fields to VTK formatted file for multiple
 * variables and timesteps. Certain predefined fields can be registered on problem / model
 * initialization and/or be turned on/off using the designated properties. Additionally
 * non-standardized scalar and vector fields can be added to the writer manually.
 */
template<typename TypeTag>
class StaggeredVtkOutputModule : public VtkOutputModuleBase<TypeTag>
{
    friend class VtkOutputModuleBase<TypeTag>;
    using ParentType = VtkOutputModuleBase<TypeTag>;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);

    using DofTypeIndices = typename GET_PROP(TypeTag, DofTypeIndices);
    typename DofTypeIndices::CellCenterIdx cellCenterIdx;
    typename DofTypeIndices::FaceIdx faceIdx;

public:

    StaggeredVtkOutputModule(const Problem& problem,
                    Dune::VTK::DataMode dm = Dune::VTK::conforming) : ParentType(problem, dm), faceWriter_(problem)

    {}

    void write(double time, Dune::VTK::OutputType type = Dune::VTK::ascii)
    {
        ParentType::write(time, type);
        faceWriter_.write();
    }


protected:
    unsigned int numDofs_() const
    {
        return this->problem().model().numCellCenterDofs();
    }

    auto getPriVarData_(const std::size_t dofIdxGlobal, const std::size_t pvIdx)
    {
        return this->problem().model().curSol()[cellCenterIdx][dofIdxGlobal][pvIdx];
    }


private:
    StaggeredVtkWriter<TypeTag> faceWriter_;
};

} // end namespace Dumux

#endif
