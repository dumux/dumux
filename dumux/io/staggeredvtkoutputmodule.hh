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

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

    using DofTypeIndices = typename GET_PROP(TypeTag, DofTypeIndices);
    typename DofTypeIndices::CellCenterIdx cellCenterIdx;
    typename DofTypeIndices::FaceIdx faceIdx;

    struct PriVarScalarDataInfo { unsigned int pvIdx; std::string name; };
    struct PriVarVectorDataInfo { std::vector<unsigned int> pvIdx; std::string name; };

    using Positions = std::vector<Scalar>;
    using Data = std::vector<std::vector<Scalar>>;

public:

    StaggeredVtkOutputModule(const Problem& problem,
                    Dune::VTK::DataMode dm = Dune::VTK::conforming) : ParentType(problem, dm), faceWriter_(coordinates_)

    {
        writeFaceVars_ = GET_PARAM_FROM_GROUP(TypeTag, bool, Vtk, WriteFaceData);
        coordinatesInitialized_ = false;
    }

    //////////////////////////////////////////////////////////////////////////////////////////////
    //! Methods to conveniently add primary and secondary variables upon problem initialization
    //! Do not call these methods after initialization
    //////////////////////////////////////////////////////////////////////////////////////////////


    //! Output a scalar field
    //! \param name The name of the vtk field
    //! \returns A reference to the resized scalar field to be filled with the actual data
    std::vector<Scalar>& createFaceScalarField(const std::string& name)
    {
        faceScalarFields_.emplace_back(std::make_pair(std::vector<Scalar>(this->problem().model().numFaceDofs()), name));
        return faceScalarFields_.back().first;
    }

    //! Output a vector field
    //! \param name The name of the vtk field
    //! \returns A reference to the resized vector field to be filled with the actual data
    std::vector<GlobalPosition>& createFaceVectorField(const std::string& name)
    {
        faceVectorFields_.emplace_back(std::make_pair(std::vector<GlobalPosition>(this->problem().model().numFaceDofs()), name));
        return faceVectorFields_.back().first;
    }

    //! Output a scalar primary variable
    //! \param name The name of the vtk field
    //! \param pvIdx The index in the primary variables vector
    void addFacePrimaryVariable(const std::string& name, unsigned int pvIdx)
    {
        priVarScalarDataInfo_.push_back(PriVarScalarDataInfo{pvIdx, name});
    }

    //! Output a vector primary variable
    //! \param name The name of the vtk field
    //! \param pvIndices A vector of indices in the primary variables vector to group for vector visualization
    void addFacePrimaryVariable(const std::string& name, std::vector<unsigned int> pvIndices)
    {
        assert(pvIndices.size() < 4 && "Vtk doesn't support vector dimensions greater than 3!");
        priVarVectorDataInfo_.push_back(PriVarVectorDataInfo{pvIndices, name});
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

    std::vector<PriVarScalarDataInfo> priVarScalarDataInfo_;
    std::vector<PriVarVectorDataInfo> priVarVectorDataInfo_;
    std::vector<PriVarScalarDataInfo> secondVarScalarDataInfo_;
    std::vector<PriVarVectorDataInfo> secondVarVectorDataInfo_;

private:

    void updateCoordinates_()
    {
        std::cout << "updating coordinates" << std::endl;
        coordinates_.resize(this->problem().model().numFaceDofs());
        for(auto&& facet : facets(this->problem().gridView()))
        {
            const int dofIdxGlobal = this->problem().gridView().indexSet().index(facet);
            coordinates_[dofIdxGlobal] = facet.geometry().center();
        }
        coordinatesInitialized_ = true;
    }

     /*!
     * \brief Gathers all face-related data and invokes the face vtk-writer using these data.
     */
    void getFaceDataAndWrite_()
    {
        const int numPoints = this->problem().model().numFaceDofs();

        // make sure not to iterate over the same dofs twice
        std::vector<bool> dofVisited(numPoints, false);

        // get fields for all primary coordinates and variables
        if(!coordinatesInitialized_)
            updateCoordinates_();

        Data priVarScalarData(priVarScalarDataInfo_.size(), std::vector<Scalar>(numPoints));

        Data priVarVectorData(priVarVectorDataInfo_.size());
        for (std::size_t i = 0; i < priVarVectorDataInfo_.size(); ++i)
            priVarVectorData[i].resize(numPoints*priVarVectorDataInfo_[i].pvIdx.size());

        for(auto&& element : elements(this->problem().gridView()))
        {
            auto fvGeometry = localView(this->problem().model().globalFvGeometry());
            fvGeometry.bindElement(element);
            for(auto && scvf : scvfs(fvGeometry))
            {
                if(dofVisited[scvf.dofIndex()])
                    continue;

                asImp_().getPrivarScalarData_(priVarScalarData, scvf);
                asImp_().getPrivarVectorData_(priVarVectorData, scvf);

                dofVisited[scvf.dofIndex()] = true;
            }
        }

        // transfer priVar scalar data to writer
        for(int i = 0; i < priVarScalarDataInfo_.size(); ++i)
            faceWriter_.addPointData(priVarScalarData[i], priVarScalarDataInfo_[i].name);

        // transfer priVar vector data to writer
        for(int i = 0; i < priVarVectorDataInfo_.size(); ++i)
            faceWriter_.addPointData(priVarVectorData[i], priVarVectorDataInfo_[i].name, priVarVectorDataInfo_[i].pvIdx.size());

        // transfer custom scalar data to writer
        for(auto&& scalarField : faceScalarFields_)
            faceWriter_.addPointData(scalarField.first, scalarField.second);

        // transfer custom vector data to writer
        for(auto&& vectorField : faceVectorFields_)
            faceWriter_.addPointData(vectorField.first, vectorField.second, 3);

        faceWriter_.write();
        faceScalarFields_.clear();
        faceVectorFields_.clear();
    }


     /*!
     * \brief Retrives scalar-valued data from the face.
     *
     * \param priVarScalarData Container to store the data
     * \param face The face
     */
    template<class Face>
    void getPrivarScalarData_(Data& priVarScalarData, const Face& face)
    {
        const int dofIdxGlobal = face.dofIndex();
        for(int pvIdx = 0; pvIdx < priVarScalarDataInfo_.size(); ++pvIdx)
            priVarScalarData[pvIdx][dofIdxGlobal] = this->problem().model().curSol()[faceIdx][dofIdxGlobal][pvIdx];
    }

     /*!
     * \brief Retrives vector-valued data from the face.
     *
     * \param priVarVectorData Container to store the data
     * \param face The face
     */
    template<class Face>
    void getPrivarVectorData_(Data& priVarVectorData, const Face& face)
    {
        const int dofIdxGlobal = face.dofIndex();
        for (int i = 0; i < priVarVectorDataInfo_.size(); ++i)
            for (int j = 0; j < priVarVectorDataInfo_[i].pvIdx.size(); ++j)
                priVarVectorData[i][dofIdxGlobal*priVarVectorDataInfo_[i].pvIdx.size() + j]
                    = this->problem().model().curSol()[faceIdx][dofIdxGlobal][priVarVectorDataInfo_[i].pvIdx[j]];
    }

    StaggeredVtkWriter<dimWorld> faceWriter_;
    bool writeFaceVars_;

    std::vector<GlobalPosition> coordinates_;
    bool coordinatesInitialized_;

    std::list<std::pair<std::vector<Scalar>, std::string>> faceScalarFields_;
    std::list<std::pair<std::vector<GlobalPosition>, std::string>> faceVectorFields_;

    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }
};

} // end namespace Dumux

#endif
