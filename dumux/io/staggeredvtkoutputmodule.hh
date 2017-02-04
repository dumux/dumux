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
#define STAGGERED_VTK_OUTPUT_MODULE_HH

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
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);

    enum { dim = GridView::dimension };
    enum { dimWorld = GridView::dimensionworld };

    using DofTypeIndices = typename GET_PROP(TypeTag, DofTypeIndices);
    typename DofTypeIndices::CellCenterIdx cellCenterIdx;
    typename DofTypeIndices::FaceIdx faceIdx;

    struct PriVarScalarDataInfo { unsigned int pvIdx; std::string name; };
    struct PriVarVectorDataInfo { std::vector<unsigned int> pvIdx; std::string name; };

public:

    StaggeredVtkOutputModule(const Problem& problem,
                    Dune::VTK::DataMode dm = Dune::VTK::conforming) : ParentType(problem, dm), faceWriter_(problem)

    {}

    //////////////////////////////////////////////////////////////////////////////////////////////
    //! Methods to conveniently add primary and secondary variables upon problem initialization
    //! Do not call these methods after initialization
    //////////////////////////////////////////////////////////////////////////////////////////////

    //! Output a scalar primary variable
    //! \param name The name of the vtk field
    //! \param pvIdx The index in the primary variables vector
    void addFacePrimaryVariable(const std::string& name, unsigned int pvIdx)
    {
        facePriVarScalarDataInfo_.push_back(PriVarScalarDataInfo{pvIdx, name});
    }

    //! Output a vector primary variable
    //! \param name The name of the vtk field
    //! \param pvIndices A vector of indices in the primary variables vector to group for vector visualization
    void addFacePrimaryVariable(const std::string& name, std::vector<unsigned int> pvIndices)
    {
        assert(pvIndices.size() < 4 && "Vtk doesn't support vector dimensions greater than 3!");
        facePriVarVectorDataInfo_.push_back(PriVarVectorDataInfo{pvIndices, name});
    }

    void write(double time, Dune::VTK::OutputType type = Dune::VTK::ascii)
    {
        ParentType::write(time, type);
        // faceWriter_.write(facePriVarScalarDataInfo_, facePriVarVectorDataInfo_);
        getData_();
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

    void getData_()
    {
        const int numPoints = this->problem().model().numFaceDofs();

        // get fields for all primary coordinates and variables
        std::vector<bool> dofVisited(numPoints, false);
        std::vector<Scalar> positions(numPoints*3);
        std::vector<std::vector<Scalar>> priVarScalarData(facePriVarScalarDataInfo_.size(), std::vector<Scalar>(numPoints));

        std::vector<std::vector<Scalar>> priVarVectorData(facePriVarVectorDataInfo_.size());
        for (std::size_t i = 0; i < facePriVarVectorDataInfo_.size(); ++i)
            priVarVectorData[i].resize(numPoints*facePriVarVectorDataInfo_[i].pvIdx.size());

        for(auto&& element : elements(this->problem().gridView()))
        {
            auto fvGeometry = localView(this->problem().model().globalFvGeometry());
            fvGeometry.bindElement(element);
            for(auto && scvf : scvfs(fvGeometry))
            {
                if(dofVisited[scvf.dofIndexSelf()])
                    continue;

                getPositions_(positions, scvf);
                getScalarData_(priVarScalarData, scvf);
                getVectorData_(priVarVectorData, scvf);
                dofVisited[scvf.dofIndexSelf()] = true;
            }
        }
        faceWriter_.write(priVarScalarData, priVarVectorData, positions, facePriVarScalarDataInfo_, facePriVarVectorDataInfo_);
    }

    template<class Facet>
    void getPositions_(std::vector<Scalar>& positions, const Facet& facet)
    {
        const int dofIdxGlobal = facet.dofIndexSelf();
        auto pos = facet.center();
        if(dim == 1)
        {
            positions[dofIdxGlobal*3] = pos[0];
            positions[dofIdxGlobal*3 + 1] = 0.0;
            positions[dofIdxGlobal*3 + 2] = 0.0;
        }
        else if(dim == 2)
        {
            positions[dofIdxGlobal*3] = pos[0];
            positions[dofIdxGlobal*3 + 1] = pos[1];
            positions[dofIdxGlobal*3 + 2] = 0.0;
        }
        else
        {
            positions[dofIdxGlobal*3] = pos[0];
            positions[dofIdxGlobal*3 + 1] = pos[1] ;
            positions[dofIdxGlobal*3 + 2] = pos[2] ;
        }
    }

    template<class Facet>
    void getScalarData_(std::vector<std::vector<Scalar>>& priVarScalarData, const Facet& facet)
    {
        const int dofIdxGlobal = facet.dofIndexSelf();
        for(int pvIdx = 0; pvIdx < facePriVarScalarDataInfo_.size(); ++pvIdx)
        {
            priVarScalarData[pvIdx][dofIdxGlobal] = this->problem().model().curSol()[faceIdx][dofIdxGlobal][pvIdx];
        }
    }

    template<class Facet>
    void getVectorData_(std::vector<std::vector<Scalar>>& priVarVectorData, const Facet& facet)
    {

        // const int dofIdxGlobal = this->problem().gridView().indexSet().index(facet);
        // for (int i = 0; i < facePriVarVectorDataInfo_.size(); ++i)
        //     for (int j = 0; j < facePriVarVectorDataInfo_[i].pvIdx.size(); ++j)
        //         priVarVectorData[i][dofIdxGlobal*facePriVarVectorDataInfo_[i].pvIdx.size() + j]
        //             = this->problem().model().curSol()[faceIdx][dofIdxGlobal][facePriVarVectorDataInfo_[i].pvIdx[j]];

        // TODO: put this into separate class, where the function is overloaded, this is a special case for staggered free flow
        const int dofIdxGlobal = facet.dofIndexSelf();
        const int dirIdx = directionIndex(facet.unitOuterNormal());
        const Scalar velocity = this->problem().model().curSol()[faceIdx][dofIdxGlobal][0];
        for (int i = 0; i < facePriVarVectorDataInfo_.size(); ++i)
            priVarVectorData[i][dofIdxGlobal*facePriVarVectorDataInfo_[i].pvIdx.size() + dirIdx] = velocity;
    }


    StaggeredVtkWriter<TypeTag, std::vector<PriVarScalarDataInfo>, std::vector<PriVarVectorDataInfo>> faceWriter_;
    std::vector<PriVarScalarDataInfo> facePriVarScalarDataInfo_;
    std::vector<PriVarVectorDataInfo> facePriVarVectorDataInfo_;
};

} // end namespace Dumux

#endif
