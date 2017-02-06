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
 * \brief A VTK output module to simplify writing dumux simulation data to VTK format. Specialization for staggered grids with dofs on faces.
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
 *        Specialization for staggered grids with dofs on faces.
 */
template<typename TypeTag>
class StaggeredVtkOutputModule : public VtkOutputModuleBase<TypeTag>
{
    friend class VtkOutputModuleBase<TypeTag>;
    using ParentType = VtkOutputModuleBase<TypeTag>;
    using Implementation = typename GET_PROP_TYPE(TypeTag, VtkOutputModule);
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

    using Positions = std::vector<Scalar>;
    using Data = std::vector<std::vector<Scalar>>;

    // a collection of types and variables to be passed to the staggered vtk writer
    struct WriterData
    {
        using Positions = StaggeredVtkOutputModule::Positions;
        using PriVarScalarDataInfo = StaggeredVtkOutputModule::PriVarScalarDataInfo;
        using PriVarVectorDataInfo = StaggeredVtkOutputModule::PriVarVectorDataInfo;
        using Data = StaggeredVtkOutputModule::Data;
        Data priVarScalarData;
        Data priVarVectorData;
        Data secondVarScalarData;
        Data secondVarVectorData;
        Positions positions;
        std::vector<PriVarScalarDataInfo> priVarScalarDataInfo;
        std::vector<PriVarVectorDataInfo> priVarVectorDataInfo;
        std::vector<PriVarScalarDataInfo> secondVarScalarDataInfo;
        std::vector<PriVarVectorDataInfo> secondVarVectorDataInfo;
    };

public:

    StaggeredVtkOutputModule(const Problem& problem,
                    Dune::VTK::DataMode dm = Dune::VTK::conforming) : ParentType(problem, dm), faceWriter_(problem)

    {
        writeFaceVars_ = GET_PARAM_FROM_GROUP(TypeTag, bool, Vtk, WriteFaceData);
    }

    //////////////////////////////////////////////////////////////////////////////////////////////
    //! Methods to conveniently add primary and secondary variables upon problem initialization
    //! Do not call these methods after initialization
    //////////////////////////////////////////////////////////////////////////////////////////////

    //! Output a scalar primary variable
    //! \param name The name of the vtk field
    //! \param pvIdx The index in the primary variables vector
    void addFacePrimaryVariable(const std::string& name, unsigned int pvIdx)
    {
        faceData_.priVarScalarDataInfo.push_back(PriVarScalarDataInfo{pvIdx, name});
    }

    //! Output a vector primary variable
    //! \param name The name of the vtk field
    //! \param pvIndices A vector of indices in the primary variables vector to group for vector visualization
    void addFacePrimaryVariable(const std::string& name, std::vector<unsigned int> pvIndices)
    {
        assert(pvIndices.size() < 4 && "Vtk doesn't support vector dimensions greater than 3!");
        faceData_.priVarVectorDataInfo.push_back(PriVarVectorDataInfo{pvIndices, name});
    }

    void write(double time, Dune::VTK::OutputType type = Dune::VTK::ascii)
    {
        ParentType::write(time, type);
        if(writeFaceVars_)
            getFaceDataAndWrite_();
    }

protected:

     /*!
     * \brief Returns the number of cell center dofs
     */
    unsigned int numDofs_() const
    {
        return this->problem().model().numCellCenterDofs();
    }

     /*!
     * \brief Returns priVar data from dofs not on the face.
     *
     * \param dofIdxGlobal The global dof index
     * \param pvIdx The primary variable index
     */
    auto getPriVarData_(const std::size_t dofIdxGlobal, const std::size_t pvIdx)
    {
        return this->problem().model().curSol()[cellCenterIdx][dofIdxGlobal][pvIdx];
    }

     /*!
     * \brief Returns a reference to the face data
     */
    const WriterData& faceData() const
    {
        return faceData_;
    }

private:

     /*!
     * \brief Gathers all face-related data and invokes the face vtk-writer using these data.
     */
    void getFaceDataAndWrite_()
    {
        const int numPoints = this->problem().model().numFaceDofs();

        // make sure not to iterate over the same dofs twice
        std::vector<bool> dofVisited(numPoints, false);

        // get fields for all primary coordinates and variables
        Positions positions(numPoints*3);

        Data priVarScalarData(faceData_.priVarScalarDataInfo.size(), std::vector<Scalar>(numPoints));

        Data priVarVectorData(faceData_.priVarVectorDataInfo.size());
        for (std::size_t i = 0; i < faceData_.priVarVectorDataInfo.size(); ++i)
            priVarVectorData[i].resize(numPoints*faceData_.priVarVectorDataInfo[i].pvIdx.size());

        for(auto&& element : elements(this->problem().gridView()))
        {
            auto fvGeometry = localView(this->problem().model().globalFvGeometry());
            fvGeometry.bindElement(element);
            for(auto && scvf : scvfs(fvGeometry))
            {
                if(dofVisited[scvf.dofIndexSelf()])
                    continue;

                asImp_().getPositions_(positions, scvf);
                asImp_().getScalarData_(priVarScalarData, scvf);
                asImp_().getVectorData_(priVarVectorData, scvf);
                dofVisited[scvf.dofIndexSelf()] = true;
            }
        }

        faceData_.priVarScalarData = std::move(priVarScalarData);
        faceData_.priVarVectorData = std::move(priVarVectorData);
        faceData_.positions = std::move(positions);
//         results.secondVarScalarData = xxx; //TODO: implemented secondVarData
//         results.secondVarScalarData = xxx;

        faceWriter_.write(faceData_);
    }

     /*!
     * \brief Retrives the position (center) of the face.
     *        ParaView expects three-dimensional coordinates, therefore we fill empty entries with zero for dim < 3.
     *
     * \param positions Container to store the positions
     * \param face The face
     */
    template<class Face>
    void getPositions_(Positions& positions, const Face& face)
    {
        const int dofIdxGlobal = face.dofIndexSelf();
        auto pos = face.center();
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

     /*!
     * \brief Retrives scalar-valued data from the face.
     *
     * \param priVarScalarData Container to store the data
     * \param face The face
     */
    template<class Face>
    void getScalarData_(Data& priVarScalarData, const Face& face)
    {
        const int dofIdxGlobal = face.dofIndexSelf();
        for(int pvIdx = 0; pvIdx < faceData_.priVarScalarDataInfo.size(); ++pvIdx)
        {
            priVarScalarData[pvIdx][dofIdxGlobal] = this->problem().model().curSol()[faceIdx][dofIdxGlobal][pvIdx];
        }
    }

     /*!
     * \brief Retrives vector-valued data from the face.
     *
     * \param priVarVectorData Container to store the data
     * \param face The face
     */
    template<class Face>
    void getVectorData_(Data& priVarVectorData, const Face& face)
    {
        const int dofIdxGlobal = face.dofIndexSelf();
        for (int i = 0; i < faceData_.priVarVectorDataInfo.size(); ++i)
            for (int j = 0; j < faceData_.priVarVectorDataInfo[i].pvIdx.size(); ++j)
                priVarVectorData[i][dofIdxGlobal*faceData_.priVarVectorDataInfo[i].pvIdx.size() + j]
                    = this->problem().model().curSol()[faceIdx][dofIdxGlobal][faceData_.priVarVectorDataInfo[i].pvIdx[j]];
    }

    StaggeredVtkWriter<TypeTag, WriterData> faceWriter_;
    WriterData faceData_;
    bool writeFaceVars_;

    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }
};

} // end namespace Dumux

#endif
