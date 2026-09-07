// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup VETest
 * \brief The immiscible 2p coarse-level VE test problem.
 */

#ifndef DUMUX_TEST_TWOPVE_PROBLEM_HH
#define DUMUX_TEST_TWOPVE_PROBLEM_HH

#include <memory>
#include <string>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/fvproblem.hh>
#include <dumux/porousmediumflow/2pve/finelevel_view.hh>

#include "problem_fine.hh"

namespace Dumux {

template<class TypeTag>
class TwoPVETestProblem : public FVProblem<TypeTag>
{
    using ParentType = FVProblem<TypeTag>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    enum {
        pressureH2OIdx = Indices::pressureIdx,
        saturationGasIdx = Indices::saturationIdx,
    };
    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;
    using ElementFluxVariablesCache = typename GridVariables::GridFluxVariablesCache::LocalView;

    using FineProblem = TwoPVEFineProblem<TypeTag>;

public:
    using SpatialParams = GetPropType<TypeTag, Properties::SpatialParams>;
    using FineLevelView = TwoPVEFineLevelView<TypeTag,FineProblem>;

    TwoPVETestProblem(std::shared_ptr<const GridGeometry> gridGeometryCoarse,
                      std::shared_ptr<FineLevelView> fineLevelView,
                      const std::string& modelParamGroup = "")
        : ParentType(gridGeometryCoarse, modelParamGroup),
          problemName_(getParamFromGroup<std::string>(modelParamGroup, "Problem.Name")),
          fineLevelView_(fineLevelView),
          spatialParams_(std::make_shared<SpatialParams>(gridGeometryCoarse, fineLevelView_->columnMap(), fineLevelView_->spatialParamsPtr(), fineLevelView_->fineCellHeight()))
    {
        calcPermeabilityAndPorosityCoarseInSpatialParams_();
    }

    /*!
     * \brief Getter function for the problem name
     */
    const std::string& name() const
    {
        return problemName_;
    }

    /*!
     * \brief Specifies which kind of boundary condition should be used for which equation on a given boundary segment
     *
     * \param globalPos the global position
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        if (onRightBoundary_(globalPos))
            values.setAllDirichlet();
        else
            values.setAllNeumann();
        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet boundary segment.
     *
     * \param element the element for which the Dirichlet boundary condition is set
     * \param scvf    the boundary sub control volume face
     */
    PrimaryVariables dirichlet(const Element &element, const SubControlVolumeFace &scvf) const
    {
        PrimaryVariables values(0.0);

        const unsigned int columnIdx = this->gridGeometry().elementMapper().index(element);
        const Scalar deltaZ = fineLevelView_->fineCellHeight();
        const auto& column = fineLevelView_->columnMap().column(columnIdx);

        // evaluate pressure at bottom of column
        GlobalPosition globalPosFineElementBottom = column[0].geometry().center();
        // shift by 0.5*cellHeight downwards, to evaluate p at the bottom of column and not center of fine cell
        globalPosFineElementBottom[dim-1] -= 0.5*deltaZ;
        values[pressureH2OIdx] = fineLevelView_->problem().dirichletAtPos(globalPosFineElementBottom)[pressureH2OIdx];

        for(const auto& fineElement : column)
        {
            GlobalPosition globalPosFineElement = fineElement.geometry().center();
            globalPosFineElement[0] = scvf.center()[0];
            // integration of fine-level saturation to obtain coarse-level saturation
            values[saturationGasIdx] += (fineLevelView_->problem().dirichletAtPos(globalPosFineElement)[saturationGasIdx] * fineLevelView_->spatialParams().porosityAtElement(fineElement))*deltaZ;
        }

        //average saturation by coarse-level porosity
        values[saturationGasIdx] /= this->spatialParams().porosityCoarseAtElement(element);

        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann control volume.
     *
     * \param element           the element for which the Neumann boundary condition is set
     * \param fvGeometry        the fvGeometry
     * \param elemVolVars       the element volume variables
     * \param elemFluxVarsCache flux variables caches for all faces in stencil
     * \param scvf              the boundary sub control volume face
     */
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVariablesCache& elemFluxVarsCache,
                        const SubControlVolumeFace& scvf) const
    {
        NumEqVector values(0.0);
        GlobalPosition globalPosCoarseScvf = scvf.center();
        const Scalar deltaZ = fineLevelView_->fineCellHeight();
        const unsigned int columnIdx = this->gridGeometry().elementMapper().index(element);
        const auto& column = fineLevelView_->columnMap().column(columnIdx);

        if(onLeftBoundary_(globalPosCoarseScvf))
        {
            for(const auto& fineElement : column)
            {
                GlobalPosition globalPosFineElement = fineElement.geometry().center();
                globalPosFineElement[0] = scvf.center()[0];
                values += fineLevelView_->problem().neumannAtPos(globalPosFineElement)*deltaZ;
            }
        }

        return values;
    }

    /*!
     * \brief Evaluates the initial values for a control volume.
     *
     * \param element the element for which the initial condition is set
     */
    PrimaryVariables initial(const Element &element) const
    {
        PrimaryVariables values(0.0);
        const Scalar deltaZ = fineLevelView_->fineCellHeight();
        const unsigned int columnIdx = this->gridGeometry().elementMapper().index(element);
        const auto& column = fineLevelView_->columnMap().column(columnIdx);

        // evaluate pressure at bottom of column
        GlobalPosition globalPosFineElementBottom = column[0].geometry().center();
        // shift by 0.5*cellHeight downwards, to evaluate p at the bottom of column and not center of fine cell
        globalPosFineElementBottom[dim-1] -= 0.5*deltaZ;
        values[pressureH2OIdx] = fineLevelView_->problem().dirichletAtPos(globalPosFineElementBottom)[pressureH2OIdx];

        for(const auto& fineElement : column)
        {
            GlobalPosition globalPosFineElement = fineElement.geometry().center();
            // integration of fine-level saturation to obtain coarse-level saturation
            values[saturationGasIdx] += (fineLevelView_->problem().initialAtPos(globalPosFineElement)[saturationGasIdx] * fineLevelView_->spatialParams().porosityAtElement(fineElement))*deltaZ;
        }

        //average saturation by coarse-level porosity
        values[saturationGasIdx] /= this->spatialParams().porosityCoarseAtElement(element);

        return values;
    }

    /*!
     * \brief Getter function for pointer to fine-level view
     */
    const std::shared_ptr<FineLevelView> getFineLevelView() const
    {
        return fineLevelView_;
    }

    /*!
     * \brief Return a reference to the coarse-level spatial parameters
     */
    const SpatialParams& spatialParams() const
    { return *spatialParams_; }

    /*!
     * \brief Return a reference to the coarse-level spatial parameters
     */
    SpatialParams& spatialParams()
    { return *spatialParams_; }


private:

    static constexpr Scalar eps_ = 1e-6;

    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_;
    }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_;
    }

    /*!
     * \brief Computes the coarse-level permeability and porosity
     */
    void calcPermeabilityAndPorosityCoarseInSpatialParams_()
    {
        this->spatialParams().calcSpatialParamsPermeabilityAndPorosityCoarse();
    }

    std::string problemName_;
    std::shared_ptr<FineLevelView> fineLevelView_;
    std::shared_ptr<SpatialParams> spatialParams_;
};

} // end namespace Dumux

#endif
