// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//

/*!
 * \file
 * \ingroup Geometry
 * \ingroup Volume Averaging
 * \brief Performs Volume Averaging for various pore scale quantities
 */

#ifndef DUMUX_GEOMETRY_VOLUME_AVERAGING_HH
#define DUMUX_GEOMETRY_VOLUME_AVERAGING_HH

#include <tuple>
#include <cmath>

#include <dumux/parallel/multithreading.hh>
#include <dumux/parallel/parallel_for.hh>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/geometry/type.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/geometry/referenceelements.hh>

#include <dumux/common/entitymap.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/geometry/geometricentityset.hh>
#include <dumux/geometry/intersectingentities.hh>

namespace Dumux::VolumeAveraging {

/////////////////////////////////////////////////////////
/////////////// Unit and Sub Averages ///////////////////
/////////////////////////////////////////////////////////

template<class SourceGridGeometry, class AveragedGridGeometry>
class UnitAveragingManager
{
public:
    using GridView = typename SourceGridGeometry::GridView;
    using Scalar = typename GridView::Grid::ctype;
    static constexpr int dim = GridView::Grid::dimension;
    using GlobalPosition = Dune::FieldVector<Scalar, dim>;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;

    UnitAveragingManager(std::shared_ptr<const SourceGridGeometry> sourceGridGeometry,
                         std::shared_ptr<const AveragedGridGeometry> averagedGridGeometry)
    : sourceGridGeometry_(sourceGridGeometry)
    , averagedGridGeometry_(averagedGridGeometry)
    {
        std::cout << "\nCreating a Unit Averaging Manager with the following parameters: \n";
        std::cout << "A Source Grid Geometry with size " << sourceGridGeometry_->elementMapper().size() << ",\n";
        std::cout << "an Averaged Grid Geometry with size " << averagedGridGeometry_->elementMapper().size() << ". \n\n";

        initializeUnitAveragingSets_();
    }

    //! Volume average the scalar solution on a base geometry to a coarser average geometry.
    template<class Content>
    const std::vector<Content> unitAverageSolution(const std::vector<Content>& solution,
                                                   const bool useFullVolume = true)
    {

        GridIndexType numTargetCells = averagedGridGeometry_->elementMapper().size();
        std::vector<Content> averagedSolution(numTargetCells, Content(0.0));
        for (int i = 0; i < numTargetCells; i++)
        {
            auto indexSet = indexMappingSet_[i];
            std::vector<Content> reducedSolVector(indexSet.size());
            // Use std::transform to copy the selected numbers
            std::transform(indexSet.begin(), indexSet.end(), reducedSolVector.begin(),
                           [&solution](int index){ return solution[index]; });

            // Sum the selected numbers using std::accumulate
            Content init(0.0);
            Content total = std::accumulate(reducedSolVector.begin(), reducedSolVector.end(), init) * sourceElementVolume_;
            Scalar volume = useFullVolume ? averagedElementVolume_ : averagingVolumeSet_[i] ;
            averagedSolution[i] = total / volume;
        }
        return averagedSolution;
    }

    template <class Content>
    std::vector<Content> perturbation(const std::vector<Content>& sourceSolution,
                                      const std::vector<Content>& averageSolution)
    {
        if ( (sourceSolution.size() != sourceGridGeometry_->elementMapper().size())
          || (averageSolution.size() != averagedGridGeometry_->elementMapper().size()))
            DUNE_THROW(Dune::InvalidStateException, "Solution sets must match the grid sizes");

        std::vector<Content> perturbation(sourceGridGeometry_->elementMapper().size(), Content(0.0));
        for (GridIndexType i = 0; i < averagedGridGeometry_->elementMapper().size(); i++)
        {
            auto indexSet = indexMappingSet_[i];
            for (GridIndexType eIdx = 0; eIdx < indexSet.size(); eIdx++)
                perturbation[eIdx] = sourceSolution[eIdx] - averageSolution[i];
        }
        return perturbation;
    }

    template <class Content>
    std::vector<Content> perturbationProduct(const std::vector<Scalar>& scalarPerturbation,
                                            const std::vector<Content>& vectorPerturbation)
    {
        if ( (scalarPerturbation.size() != sourceGridGeometry_->elementMapper().size())
          || (vectorPerturbation.size() != sourceGridGeometry_->elementMapper().size()))
            DUNE_THROW(Dune::InvalidStateException, "Perturbation sets must match the grid size");

        std::vector<Content> perturbationProduct(sourceGridGeometry_->elementMapper().size(), Content(0.0));
        for (const auto& sourceElement : elements(sourceGridGeometry_->gridView()))
        {
            GridIndexType eIdx = sourceGridGeometry_->elementMapper().index(sourceElement);
            perturbationProduct[eIdx] = (scalarPerturbation[eIdx] * vectorPerturbation[eIdx]);
        }
        return perturbationProduct;
    }

    //! Extend the scalar solution on a subgrid to a hostgrid with zeros for non-contained cells.
    template <class Content>
    std::vector<Content> extendSolutionToAveragedGrid(const std::vector<Content>& solution)
    { return unitAverageSolution(solution, true); }

private:

    void initializeUnitAveragingSets_()
    {
        sourceElementVolume_ = sourceGridGeometry_->element(0).geometry().volume();
        averagedElementVolume_ = averagedGridGeometry_->element(0).geometry().volume();
        checkGrids_();

        indexMappingSet_.resize(averagedGridGeometry_->elementMapper().size());
        averagingVolumeSet_.resize(averagedGridGeometry_->elementMapper().size());

        // Set up a mapping from the source gridGeometry to the averaged gridGeometry
        for (const auto& averageElement : elements(averagedGridGeometry_->gridView()))
        {
            int hostEIdx = averagedGridGeometry_->elementMapper().index(averageElement);
            const auto entities = intersectingEntities(averageElement.geometry(), sourceGridGeometry_->boundingBoxTree());
            for(int i=0; i<entities.size()/4; i++)
            {
                const auto entity = entities[4*i];
                const auto entityIdx = entity.second();
                const auto baseElement = sourceGridGeometry_->element(entityIdx);
                indexMappingSet_[hostEIdx].push_back(entityIdx);
                averagingVolumeSet_[hostEIdx] += baseElement.geometry().volume();
            }
        }
    }

    void checkGrids_()
    {
        // Check the grids first to make sure they are useable.
        for (const auto& sourceElement : elements(sourceGridGeometry_->gridView()))
        {
            if ( (sourceElement.geometry().volume() > sourceElementVolume_ + 1e-10) ||
                 (sourceElement.geometry().volume() < sourceElementVolume_ - 1e-10) )
                DUNE_THROW(Dune::InvalidStateException, "Non-uniform grids are not currently supported");
        }
        for (const auto& averageElement : elements(averagedGridGeometry_->gridView()))
        {
            if ( (averageElement.geometry().volume() > averagedElementVolume_ + 1e-10) ||
                 (averageElement.geometry().volume() < averagedElementVolume_ - 1e-10) )
                DUNE_THROW(Dune::InvalidStateException, "Non-uniform grids are not currently supported");
        }
    }

    std::shared_ptr<const SourceGridGeometry> sourceGridGeometry_;
    std::shared_ptr<const AveragedGridGeometry> averagedGridGeometry_;

    std::vector<std::vector<GridIndexType>> indexMappingSet_;
    std::vector<Scalar> averagingVolumeSet_;
    Scalar sourceElementVolume_;
    Scalar averagedElementVolume_;

};

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// Convolutional Averaging /////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////


template<class TargetGridGeometry, class SourceGridGeometry>
class ConvolutionalAveragingManager
{
public:
    using GridView = typename SourceGridGeometry::GridView;
    using Scalar = typename GridView::Grid::ctype;
    static constexpr int dim = GridView::Grid::dimension;
    using GlobalPosition = Dune::FieldVector<Scalar, dim>;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;

    ConvolutionalAveragingManager(std::shared_ptr<const TargetGridGeometry> targetGridGeometry,
                                  std::shared_ptr<const SourceGridGeometry> sourceGridGeometry,
                                  const Scalar& filterSize,
                                  const bool& useCircleFilter = false)
    : targetGridGeometry_(targetGridGeometry)
    , sourceGridGeometry_(sourceGridGeometry)
    , filterSize_(filterSize)
    , useCircleFilter_(useCircleFilter)
    {
        report_();

        if (useCircleFilter)
            initializeCircularConvolutionalFilter_();
        else
            initializeSquareConvolutionalFilter_();
    }

    ConvolutionalAveragingManager(const std::string& pathName,
                                  std::shared_ptr<const TargetGridGeometry> targetGridGeometry,
                                  std::shared_ptr<const SourceGridGeometry> sourceGridGeometry,
                                  const Scalar& filterSize,
                                  const bool& useCircleFilter = false)
    : targetGridGeometry_(targetGridGeometry)
    , sourceGridGeometry_(sourceGridGeometry)
    , filterSize_(filterSize)
    , useCircleFilter_(useCircleFilter)
    {
        report_(pathName);
        readIndexDataFromFile_(pathName);
    }

    template<class Content>
    const std::vector<Content> averageSolution(std::vector<Content>& solutionVector)
    {
        if (solutionVector.size() != sourceGridGeometry_->elementMapper().size())
            DUNE_THROW(Dune::InvalidStateException, "Size mismatch, the solution vector should have the same size as the sourceGrid");

        GridIndexType numTargetCells = targetGridGeometry_->elementMapper().size();
        std::vector<Content> averagedSolutionSet(numTargetCells);
        for (int i = 0; i < numTargetCells; i++)
        {
            auto indexSet = indexMapping_[i];
            std::vector<Content> reducedSolVector(indexSet.size());
            // Use std::transform to copy the selected numbers
            std::transform(indexSet.begin(), indexSet.end(), reducedSolVector.begin(),
                           [&solutionVector](int index){ return solutionVector[index]; });

            // Sum the selected numbers using std::accumulate
            Content init(0.0);
            Content total = std::accumulate(reducedSolVector.begin(), reducedSolVector.end(), init);
            averagedSolutionSet[i] = total / (indexSet.size());
        }

        return averagedSolutionSet;
    }

    const std::vector<GridIndexType> indexCountPerCell()
    {
        GridIndexType numTargetCells = targetGridGeometry_->elementMapper().size();
        std::vector<GridIndexType> indexCount(numTargetCells, 0);

        for (int i = 0; i < numTargetCells; i++)
            indexCount[i] = indexMapping_[i].size();

        return indexCount;
    }

    template <class Content>
    std::vector<Content> perturbation(const std::vector<Content>& sourceSolution,
                                      const std::vector<Content>& averageSolution)
    {
        if (sourceSolution.size() != averageSolution.size() ||
            averageSolution.size() != targetGridGeometry_->elementMapper().size())
            DUNE_THROW(Dune::InvalidStateException, "To evaluate the perturbation via the convolutionalManager,"
                                                 << " the source and average solutions should be the same size");

        std::vector<Content> perturbation(targetGridGeometry_->elementMapper().size(), Content(0.0));
        for (GridIndexType i = 0; i < targetGridGeometry_->elementMapper().size(); i++)
            perturbation[i] = sourceSolution[i] - averageSolution[i];

        return perturbation;
    }

    template <class Content>
    std::vector<Content> perturbationProduct(const std::vector<Scalar>& scalarPerturbation,
                                             const std::vector<Content>& vectorPerturbation)
    {
        if (scalarPerturbation.size() != vectorPerturbation.size() ||
            scalarPerturbation.size() != targetGridGeometry_->elementMapper().size())
            DUNE_THROW(Dune::InvalidStateException, "To evaluate the perturbation via the convolutionalManager,"
                                                 << " the source and average solutions should be the same size");

        std::vector<Content> perturbationProduct(targetGridGeometry_->elementMapper().size(), Content(0.0));
        for (GridIndexType i = 0; i < targetGridGeometry_->elementMapper().size(); i++)
            perturbationProduct[i] = scalarPerturbation[i] * vectorPerturbation[i];

        return perturbationProduct;
    }

    const std::vector<std::vector<GridIndexType>> indexMapping()
    { return indexMapping_; }

    const Scalar filterSize()
    { return filterSize_; }

    const std::string filterType()
    { return (useCircleFilter_) ? "Circle" : "Square"; }

    void writeIndexDataToFile(const std::string& filename)
    {
        std::ofstream file(filename, std::ios::binary);
        if (!file.is_open())
            DUNE_THROW(Dune::InvalidStateException, "File opening error");

        // Iterate over the index data and write each element
        for (const auto& row : indexMapping_)
        {
            // Write the size of the row
            std::size_t rowSize = row.size();
            file.write(reinterpret_cast<const char*>(&rowSize), sizeof(rowSize));

            // Write the elements for the current row
            for (const auto& element : row)
                file.write(reinterpret_cast<const char*>(&element), sizeof(element));
        }
    }

private:

    void readIndexDataFromFile_(const std::string& filename)
    {
        // Create one set of source grid indexes per target Grid Element
        indexMapping_.resize(targetGridGeometry_->elementMapper().size());

        std::ifstream file(filename, std::ios::binary);
        if (!file.is_open())
            DUNE_THROW(Dune::InvalidStateException, "File opening error");

        GridIndexType element;
        int index = 0;
        std::size_t rowSize;

        while (file.read(reinterpret_cast<char*>(&rowSize), sizeof(rowSize)))
        {
            std::vector<GridIndexType> row(rowSize);

            // Read the elements for the current row
            for (std::size_t i = 0; i < rowSize; ++i)
            {
                if (!file.read(reinterpret_cast<char*>(&element), sizeof(element)))
                    DUNE_THROW(Dune::InvalidStateException, "error at index: " + std::to_string(index));
                row[i] = element;
            }

            indexMapping_[index] = row;
            index++;
        }
    }

    void initializeSquareConvolutionalFilter_()
    {
        // Create one set of source grid indexes per target Grid Element
        indexMapping_.resize(targetGridGeometry_->elementMapper().size());

        // Extend the radius check by one half cell.
        GlobalPosition halfSquareOffset = GlobalPosition((filterSize_ / 2.0));
        GlobalPosition eps = GlobalPosition(1e-6);
        elementVolume_ = targetGridGeometry_->element(0).geometry().volume(); // Store the size one one element

        //for all averaged elements (averaged grid)
        Dumux::parallelFor(targetGridGeometry_->gridView().size(0), [&](const std::size_t targetElementIdx)
        {
            const auto targetElement = targetGridGeometry_->element(targetElementIdx);
            const auto bBoxMax = targetElement.geometry().center() + halfSquareOffset - eps;
            const auto bBoxMin = targetElement.geometry().center() - halfSquareOffset + eps;

            std::vector<GridIndexType> indexSet;
            for (const auto& sourceElement : elements(sourceGridGeometry_->gridView()))
            {
                bool withinFilter = true;
                for (int dimIdx = 0; dimIdx < dim; dimIdx++)
                {
                    GlobalPosition sourceCenter = sourceElement.geometry().center();
                    if ((sourceCenter[dimIdx] < bBoxMin[dimIdx]) || (sourceCenter[dimIdx] > bBoxMax[dimIdx]))
                        withinFilter = false;
                }

                if (withinFilter)
                    indexSet.push_back(sourceGridGeometry_->elementMapper().index(sourceElement));
            }
            indexMapping_[targetElementIdx] = indexSet;
        });
    }

    void initializeCircularConvolutionalFilter_()
    {
        // Create one set of source grid indexes per target Grid Element
        indexMapping_.resize(targetGridGeometry_->elementMapper().size());

        // Extend the radius check by one half cell.
        auto cellGeometry = targetGridGeometry_->element(0).geometry();
        elementVolume_ = cellGeometry.volume(); // Store the size one one element
        Scalar halfCellLength = (cellGeometry.corner(0) - cellGeometry.corner(1)).two_norm() * 0.5;

        //for all averaged elements (averaged grid)
        Dumux::parallelFor(targetGridGeometry_->gridView().size(0), [&](const std::size_t targetElementIdx)
        {
            const auto targetElement = targetGridGeometry_->element(targetElementIdx);
            std::vector<GridIndexType> indexSet;
            for (const auto& sourceElement : elements(sourceGridGeometry_->gridView()))
            {
                if ( (targetElement.geometry().center() - sourceElement.geometry().center()).two_norm() > (filterSize_ + halfCellLength))
                    continue;
                else
                    indexSet.push_back(sourceGridGeometry_->elementMapper().index(sourceElement));
            }
            indexMapping_[targetElementIdx] = indexSet;
        });
    }

    void report_(const std::string pathName)
    {
        report_();
        std::cout << "The source file is stored here: " << pathName << ".\n\n";
    }

    void report_()
    {
        std::cout << "\nCreating a Convolutional Averaging Manager from a file path with the following parameters: \n";
        std::cout << "A Target Grid Geometry with size " << targetGridGeometry_->elementMapper().size() << ",\n";
        std::cout << "a Source Grid Geometry with size " << sourceGridGeometry_->elementMapper().size() << ",\n";
        std::cout << "and a " <<  (useCircleFilter_ ? "Circular" : "Square") << " Filter that includes " << filterSize_ << "in each direction. \n\n";
    }

    std::shared_ptr<const TargetGridGeometry> targetGridGeometry_;
    std::shared_ptr<const SourceGridGeometry> sourceGridGeometry_;
    Scalar filterSize_;
    bool useCircleFilter_;

    Scalar elementVolume_;
    std::vector<std::vector<GridIndexType>> indexMapping_;
};

}  // end namespace Dumux::VolumeAveraging

#endif
