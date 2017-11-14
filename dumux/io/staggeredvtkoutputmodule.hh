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

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/pointcloudvtkwriter.hh>
#include <dumux/io/vtksequencewriter.hh>


namespace Dumux
{


/*!
 * \ingroup InputOutput
 * \brief A VTK output module to simplify writing dumux simulation data to VTK format
 *        Specialization for staggered grids with dofs on faces.
 */
template<typename TypeTag>
class StaggeredVtkOutputModule : public VtkOutputModule<TypeTag>
{
    friend class VtkOutputModule<TypeTag>;
    using ParentType = VtkOutputModule<TypeTag>;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using FaceVariables = typename GET_PROP_TYPE(TypeTag, FaceVariables);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);


    enum { dim = GridView::dimension };
    enum { dimWorld = GridView::dimensionworld };

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

    using DofTypeIndices = typename GET_PROP(TypeTag, DofTypeIndices);
    typename DofTypeIndices::CellCenterIdx cellCenterIdx;
    typename DofTypeIndices::FaceIdx faceIdx;

    struct FaceVarScalarDataInfo { std::function<Scalar(const FaceVariables&)> get; std::string name; };
    struct FaceVarVectorDataInfo { std::function<GlobalPosition(const SubControlVolumeFace& scvf, const FaceVariables&)> get; std::string name; };

    using Positions = std::vector<Scalar>;
    using Data = std::vector<std::vector<Scalar>>;

public:

    StaggeredVtkOutputModule(const Problem& problem,
                    const FVGridGeometry& fvGridGeometry,
                    const GridVariables& gridVariables,
                    const SolutionVector& sol,
                    const std::string& name,
                    bool verbose = true,
                    Dune::VTK::DataMode dm = Dune::VTK::conforming)
    : ParentType(problem, fvGridGeometry, gridVariables, sol, name, verbose, dm)
    , problem_(problem)
    , gridGeom_(fvGridGeometry)
    , gridVariables_(gridVariables)
    , sol_(sol)
    , faceWriter_(std::make_shared<PointCloudVtkWriter<Scalar, dim>>(coordinates_))
    , sequenceWriter_(faceWriter_, problem.name() + "-face", "","",
                      fvGridGeometry.gridView().comm().rank(),
                      fvGridGeometry.gridView().comm().size() )

    {
        writeFaceVars_ = getParamFromGroup<bool>(GET_PROP_VALUE(TypeTag, ModelParameterGroup), "Vtk.WriteFaceData", false);
        coordinatesInitialized_ = false;
    }

    //////////////////////////////////////////////////////////////////////////////////////////////
    //! Methods to conveniently add primary and secondary variables upon problem initialization
    //! Do not call these methods after initialization
    //////////////////////////////////////////////////////////////////////////////////////////////

    template<typename Vector>
    void addFaceField(const Vector& v, const std::string& name)
    {
        static_assert(std::is_same<Vector, std::vector<Scalar>>::value ||
                      std::is_same<Vector, std::vector<GlobalPosition>>::value,
                      "Only vectors of Scalar or GlobalPosition are supported");

        if (v.size() == this->gridGeom_.gridView().size(1))
        {
            if(!coordinatesInitialized_)
                updateCoordinates_();

            faceWriter_->addPointData(v, name);
        }
        else
            DUNE_THROW(Dune::RangeError, "Size mismatch of added field!");
    }

    void addFaceVariable(std::function<Scalar(const FaceVariables&)>&& f, const std::string& name)
    {
        faceVarScalarDataInfo_.push_back(FaceVarScalarDataInfo{f, name});
    }

    void addFaceVariable(std::function<GlobalPosition(const SubControlVolumeFace& scvf, const FaceVariables&)>&& f, const std::string& name)
    {
        faceVarVectorDataInfo_.push_back(FaceVarVectorDataInfo{f, name});
    }


    void write(double time, Dune::VTK::OutputType type = Dune::VTK::ascii)
    {
        ParentType::write(time, type);
        if(writeFaceVars_)
            getFaceDataAndWrite_(time);
    }


private:

    void updateCoordinates_()
    {
        coordinates_.resize(gridGeom_.numFaceDofs());
        for(auto&& facet : facets(gridGeom_.gridView()))
        {
            const int dofIdxGlobal = gridGeom_.gridView().indexSet().index(facet);
            coordinates_[dofIdxGlobal] = facet.geometry().center();
        }
        coordinatesInitialized_ = true;
    }

     /*!
     * \brief Gathers all face-related data and invokes the face vtk-writer using these data.
     */
    void getFaceDataAndWrite_(const Scalar time)
    {
        const auto numPoints = gridGeom_.numFaceDofs();

        // make sure not to iterate over the same dofs twice
        std::vector<bool> dofVisited(numPoints, false);

        // get fields for all primary coordinates and variables
        if(!coordinatesInitialized_)
            updateCoordinates_();

        std::vector<std::vector<Scalar>> faceVarScalarData;
        std::vector<std::vector<GlobalPosition>> faceVarVectorData;

        if(!faceVarScalarDataInfo_.empty())
            faceVarScalarData.resize(faceVarScalarDataInfo_.size(), std::vector<Scalar>(numPoints));

        if(!faceVarVectorDataInfo_.empty())
            faceVarVectorData.resize(faceVarVectorDataInfo_.size(), std::vector<GlobalPosition>(numPoints));

        for (const auto& element : elements(gridGeom_.gridView(), Dune::Partitions::interior))
        {
            auto fvGeometry = localView(gridGeom_);
            auto elemFaceVars = localView(gridVariables_.curGridFaceVars());

            if (!faceVarScalarDataInfo_.empty() || !faceVarVectorDataInfo_.empty())
            {
                fvGeometry.bind(element);
                elemFaceVars.bindElement(element, fvGeometry, sol_);

                for (auto&& scvf : scvfs(fvGeometry))
                {
                    const auto dofIdxGlobal = scvf.dofIndex();
                    if(dofVisited[dofIdxGlobal])
                        continue;

                    dofVisited[dofIdxGlobal] = true;

                    const auto& faceVars = elemFaceVars[scvf];

                    for (std::size_t i = 0; i < faceVarScalarDataInfo_.size(); ++i)
                        faceVarScalarData[i][dofIdxGlobal] = faceVarScalarDataInfo_[i].get(faceVars);

                    for (std::size_t i = 0; i < faceVarVectorDataInfo_.size(); ++i)
                            faceVarVectorData[i][dofIdxGlobal] = faceVarVectorDataInfo_[i].get(scvf, faceVars);
                }
            }
        }

        if(!faceVarScalarDataInfo_.empty())
            for (std::size_t i = 0; i < faceVarScalarDataInfo_.size(); ++i)
                faceWriter_->addPointData(faceVarScalarData[i], faceVarScalarDataInfo_[i].name);

        if(!faceVarVectorDataInfo_.empty())
            for (std::size_t i = 0; i < faceVarVectorDataInfo_.size(); ++i)
                faceWriter_->addPointData(faceVarVectorData[i], faceVarVectorDataInfo_[i].name);

        sequenceWriter_.write(time);
        coordinates_.clear();
        coordinates_.shrink_to_fit();
        coordinatesInitialized_ = false;
    }


    const Problem& problem_;
    const FVGridGeometry& gridGeom_;
    const GridVariables& gridVariables_;
    const SolutionVector& sol_;

    std::shared_ptr<PointCloudVtkWriter<Scalar, dim>> faceWriter_;

    VTKSequenceWriter<PointCloudVtkWriter<Scalar, dim>> sequenceWriter_;

    bool writeFaceVars_;

    std::vector<GlobalPosition> coordinates_;
    bool coordinatesInitialized_;

    std::vector<FaceVarScalarDataInfo> faceVarScalarDataInfo_;
    std::vector<FaceVarVectorDataInfo> faceVarVectorDataInfo_;

};

} // end namespace Dumux

#endif
