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
 * \brief A VTK output module to simplify writing dumux simulation data to VTK format.
 *        This is a specialization for a staggered free-flow implementation on a regular grid.
 */
#ifndef FREEFLOW_STAGGERED_VTK_OUTPUT_MODULE_HH
#define FREEFLOW_STAGGERED_VTK_OUTPUT_MODULE_HH

#include <dune/common/fvector.hh>

#include <dumux/io/staggeredvtkoutputmodule.hh>
#include <dumux/discretization/staggered/freeflow/staggeredgeometryhelper.hh>

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
 *        This is a specialization for a staggered free-flow implementation on a regular grid.
 */
template<typename TypeTag>
class FreeFlowStaggeredVtkOutputModule : public StaggeredVtkOutputModule<TypeTag>
{
    friend class StaggeredVtkOutputModule<TypeTag>;
    using ParentType = StaggeredVtkOutputModule<TypeTag>;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);

    using DofTypeIndices = typename GET_PROP(TypeTag, DofTypeIndices);
    typename DofTypeIndices::FaceIdx faceIdx;

    using Data = std::vector<std::vector<Scalar>>;

public:
    FreeFlowStaggeredVtkOutputModule(const Problem& problem,
                    Dune::VTK::DataMode dm = Dune::VTK::conforming) : ParentType(problem, dm)

    {}

private:

     /*!
     * \brief Retrives vector-valued data from the face. This is a specialization for a free-flow implementation on a regular grid.
     *
     * \param priVarVectorData Container to store the data
     * \param face The face
     */
    template<class Face>
    void getPrivarVectorData_(Data& priVarVectorData, const Face& face)
    {
        const int dofIdxGlobal = face.dofIndex();
        const int dirIdx = directionIndex(face.unitOuterNormal());
        const Scalar velocity = this->problem().model().curSol()[faceIdx][dofIdxGlobal][0];
        for (int i = 0; i < this->priVarVectorDataInfo_.size(); ++i)
            priVarVectorData[i][dofIdxGlobal * this->priVarVectorDataInfo_[i].pvIdx.size() + dirIdx] = velocity;
    }
};

} // end namespace Dumux

#endif
