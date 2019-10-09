// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup InputOutput
 * \brief A VTK output module to simplify writing dumux simulation data to VTK format. Specialization for staggered grids with dofs on faces.
 */
#ifndef DUMUX_STAGGERED_VTK_OUTPUT_MODULE_HH
#define DUMUX_STAGGERED_VTK_OUTPUT_MODULE_HH

#include <dune/common/fvector.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/pointcloudvtkwriter.hh>
#include <dumux/io/vtksequencewriter.hh>
#include <dumux/discretization/staggered/freeflow/velocityoutput.hh>

namespace Dumux {

template<class Scalar, class GlobalPosition>
class PointCloudVtkWriter;

/*!
 * \ingroup InputOutput
 * \brief A VTK output module to simplify writing dumux simulation data to VTK format
 *        Specialization for staggered grids with dofs on faces.
 *
 * \tparam GridVariables The grid variables
 * \tparam SolutionVector The solution vector
 */
template<class GridVariables, class SolutionVector>
class StaggeredVtkOutputModule
: public VtkOutputModule<GridVariables, SolutionVector>
{
    using ParentType = VtkOutputModule<GridVariables, SolutionVector>;
    using GridGeometry = typename GridVariables::GridGeometry;
    using GridView = typename GridGeometry::GridView;
    using Scalar = typename GridVariables::Scalar;
    using FaceVariables = typename GridVariables::GridFaceVariables::FaceVariables;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    struct FaceVarScalarDataInfo { std::function<Scalar(const FaceVariables&)> get; std::string name; };
    struct FaceVarVectorDataInfo { std::function<GlobalPosition(const SubControlVolumeFace& scvf, const FaceVariables&)> get; std::string name; };

    struct FaceFieldScalarDataInfo
    {
        FaceFieldScalarDataInfo(const std::vector<Scalar>& f, const std::string& n) : data(f), name(n) {}
        const std::vector<Scalar>& data;
        const std::string name;
    };

    struct FaceFieldVectorDataInfo
    {
        FaceFieldVectorDataInfo(const std::vector<GlobalPosition>& f, const std::string& n) : data(f), name(n) {}
        const std::vector<GlobalPosition>& data;
        const std::string name;
    };

public:

    template<class Sol>
    StaggeredVtkOutputModule(const GridVariables& gridVariables,
                             const Sol& sol,
                             const std::string& name,
                             const std::string& paramGroup = "",
                             Dune::VTK::DataMode dm = Dune::VTK::conforming,
                             bool verbose = true)
    : ParentType(gridVariables, sol, name, paramGroup, dm, verbose)
    , faceWriter_(std::make_shared<PointCloudVtkWriter<Scalar, GlobalPosition>>(coordinates_))
    , sequenceWriter_(faceWriter_, name + "-face", "","",
                      gridVariables.curGridVolVars().problem().gridGeometry().gridView().comm().rank(),
                      gridVariables.curGridVolVars().problem().gridGeometry().gridView().comm().size() )

    {
        static_assert(std::is_same<Sol, SolutionVector>::value, "Make sure that sol has the same type as SolutionVector."
                                                                "Use StaggeredVtkOutputModule<GridVariables, decltype(sol)> when calling the constructor.");

        // enable velocity output per default
        this->addVelocityOutput(std::make_shared<StaggeredFreeFlowVelocityOutput<GridVariables, SolutionVector>>(gridVariables, sol));
        writeFaceVars_ = getParamFromGroup<bool>(paramGroup, "Vtk.WriteFaceData", false);
        coordinatesInitialized_ = false;
    }

    //////////////////////////////////////////////////////////////////////////////////////////////
    //! Methods to conveniently add face variables
    //! Do not call these methods after initialization
    //////////////////////////////////////////////////////////////////////////////////////////////

    //! Add a scalar valued field
    //! \param v The field to be added
    //! \param name The name of the vtk field
    void addFaceField(const std::vector<Scalar>& v, const std::string& name)
    {
        if (v.size() == this->gridGeometry().gridView().size(1))
            faceFieldScalarDataInfo_.emplace_back(v, name);
        else
            DUNE_THROW(Dune::RangeError, "Size mismatch of added field!");
    }

    //! Add a vector valued field
    //! \param v The field to be added
    //! \param name The name of the vtk field
    void addFaceField(const std::vector<GlobalPosition>& v, const std::string& name)
    {
        if (v.size() == this->gridGeometry().gridView().size(1))
            faceFieldVectorDataInfo_.emplace_back(v, name);
        else
            DUNE_THROW(Dune::RangeError, "Size mismatch of added field!");
    }

    //! Add a scalar-valued faceVarible
    //! \param f A function taking a FaceVariables object and returning the desired scalar
    //! \param name The name of the vtk field
    void addFaceVariable(std::function<Scalar(const FaceVariables&)>&& f, const std::string& name)
    {
        faceVarScalarDataInfo_.push_back(FaceVarScalarDataInfo{f, name});
    }

    //! Add a vector-valued faceVarible
    //! \param f A function taking a SubControlVolumeFace and FaceVariables object and returning the desired vector
    //! \param name The name of the vtk field
    void addFaceVariable(std::function<GlobalPosition(const SubControlVolumeFace& scvf, const FaceVariables&)>&& f, const std::string& name)
    {
        faceVarVectorDataInfo_.push_back(FaceVarVectorDataInfo{f, name});
    }

    //! Write the values to vtp files
    //! \param time The current time
    //! \param type The output type
    void write(double time, Dune::VTK::OutputType type = Dune::VTK::ascii)
    {
        ParentType::write(time, type);
        if(writeFaceVars_)
            getFaceDataAndWrite_(time);
    }


private:

    //! Update the coordinates (the face centers)
    void updateCoordinates_()
    {
        coordinates_.resize(this->gridGeometry().numFaceDofs());
        for(auto&& facet : facets(this->gridGeometry().gridView()))
        {
            const int dofIdxGlobal = this->gridGeometry().gridView().indexSet().index(facet);
            coordinates_[dofIdxGlobal] = facet.geometry().center();
        }
        coordinatesInitialized_ = true;
    }

     //! Gathers all face-related data and invokes the face vtk-writer using these data.
     //! \param time The current time
    void getFaceDataAndWrite_(const Scalar time)
    {
        const auto numPoints = this->gridGeometry().numFaceDofs();

        // make sure not to iterate over the same dofs twice
        std::vector<bool> dofVisited(numPoints, false);

        // get fields for all primary coordinates and variables
        if(!coordinatesInitialized_)
            updateCoordinates_();

        // prepare some containers to store the relevant data
        std::vector<std::vector<Scalar>> faceVarScalarData;
        std::vector<std::vector<GlobalPosition>> faceVarVectorData;

        if(!faceVarScalarDataInfo_.empty())
            faceVarScalarData.resize(faceVarScalarDataInfo_.size(), std::vector<Scalar>(numPoints));

        if(!faceVarVectorDataInfo_.empty())
            faceVarVectorData.resize(faceVarVectorDataInfo_.size(), std::vector<GlobalPosition>(numPoints));

        for (const auto& element : elements(this->gridGeometry().gridView(), Dune::Partitions::interior))
        {
            auto fvGeometry = localView(this->gridGeometry());
            auto elemFaceVars = localView(this->gridVariables().curGridFaceVars());

            if (!faceVarScalarDataInfo_.empty() || !faceVarVectorDataInfo_.empty())
            {
                fvGeometry.bind(element);
                elemFaceVars.bindElement(element, fvGeometry, this->sol());

                for (auto&& scvf : scvfs(fvGeometry))
                {
                    const auto dofIdxGlobal = scvf.dofIndex();
                    if(dofVisited[dofIdxGlobal])
                        continue;

                    dofVisited[dofIdxGlobal] = true;

                    const auto& faceVars = elemFaceVars[scvf];

                    // get the scalar-valued data
                    for (std::size_t i = 0; i < faceVarScalarDataInfo_.size(); ++i)
                        faceVarScalarData[i][dofIdxGlobal] = faceVarScalarDataInfo_[i].get(faceVars);

                    // get the vector-valued data
                    for (std::size_t i = 0; i < faceVarVectorDataInfo_.size(); ++i)
                            faceVarVectorData[i][dofIdxGlobal] = faceVarVectorDataInfo_[i].get(scvf, faceVars);
                }
            }
        }

        // transfer the data to the point writer
        if(!faceVarScalarDataInfo_.empty())
            for (std::size_t i = 0; i < faceVarScalarDataInfo_.size(); ++i)
                faceWriter_->addPointData(faceVarScalarData[i], faceVarScalarDataInfo_[i].name);

        if(!faceVarVectorDataInfo_.empty())
            for (std::size_t i = 0; i < faceVarVectorDataInfo_.size(); ++i)
                faceWriter_->addPointData(faceVarVectorData[i], faceVarVectorDataInfo_[i].name);

        // account for the custom fields
        for(const auto& field: faceFieldScalarDataInfo_)
            faceWriter_->addPointData(field.data, field.name);

        for(const auto& field: faceFieldVectorDataInfo_)
            faceWriter_->addPointData(field.data, field.name);

        // write for the current time step
        sequenceWriter_.write(time);

        // clear coordinates to save some memory
        coordinates_.clear();
        coordinates_.shrink_to_fit();
        coordinatesInitialized_ = false;
    }


    std::shared_ptr<PointCloudVtkWriter<Scalar, GlobalPosition>> faceWriter_;

    VTKSequenceWriter<PointCloudVtkWriter<Scalar, GlobalPosition>> sequenceWriter_;

    bool writeFaceVars_;

    std::vector<GlobalPosition> coordinates_;
    bool coordinatesInitialized_;

    std::vector<FaceVarScalarDataInfo> faceVarScalarDataInfo_;
    std::vector<FaceVarVectorDataInfo> faceVarVectorDataInfo_;

    std::vector<FaceFieldScalarDataInfo> faceFieldScalarDataInfo_;
    std::vector<FaceFieldVectorDataInfo> faceFieldVectorDataInfo_;

};

} // end namespace Dumux

#endif
