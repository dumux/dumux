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
 * \ingroup RANSTests
 */
#ifndef DUMUX_VERIFICATION_KOMEGA_PROBLEM_HH
#define DUMUX_VERIFICATION_KOMEGA_PROBLEM_HH

#include <math.h>

#include <dumux/common/properties.hh>
#include <dumux/freeflow/navierstokes/boundarytypes.hh>
#include <dumux/freeflow/turbulenceproperties.hh>

#include "komegaanalytical.hh"

namespace Dumux {

/*!
 * \ingroup RANSTests
 * \brief  Verification of the KOmega Model.
 */
template <class TypeTag>
class AnalyticalProblem : public RANSProblem<TypeTag>
{
    using ParentType = RANSProblem<TypeTag>;

    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using FluidState = GetPropType<TypeTag, Properties::FluidState>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;
    using BoundaryTypes = Dumux::NavierStokesBoundaryTypes<ModelTraits::numEq()>;

    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

    static constexpr auto dimWorld = GridGeometry::GridView::dimensionworld;

    using VelocityVector = Dune::FieldVector<Scalar, dimWorld>;
    using AnalyticalModel = typename Dumux::KOmegaAnalytical<TypeTag>;

public:
    AnalyticalProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , eps_(1e-6)
    {
        FluidSystem::init();
        density_ = getParam<Scalar>("Component.LiquidDensity", 1.0);
        kinematicViscosity_ = getParam<Scalar>("Component.LiquidKinematicViscosity", 1.0);

        //Set the naming convention here
        name_ = getParam<std::string>("Problem.Name");
        using CellVector = std::vector<unsigned int>;
        const auto numCellsXVec = getParam<CellVector>("Grid.Cells0");
        const auto numCellsYVec = getParam<CellVector>("Grid.Cells1");
//         name_ = name_ + std::to_string(accumulate(numCellsXVec.begin(),numCellsXVec.end(),0)) + "x"
//                       + std::to_string(accumulate(numCellsYVec.begin(),numCellsYVec.end(),0)) + "y";

        turbulenceModelName_ = turbulenceModelToString(ModelTraits::turbulenceModel());
        std::cout << "Using the "<< turbulenceModelName_ << " Turbulence Model. \n";
        std::cout << std::endl;

        createManufacturedSolution_();
    }

   // Returns the temperature within the domain in [K]
    Scalar temperature() const
    { return 283.15; }

    // Set the name via parameter in input-file
    const std::string& name() const
    { return name_; }

    bool isOnWallAtPos(const GlobalPosition &globalPos) const
    {
        return globalPos[1] < this->gridGeometry().bBoxMin()[1] + eps_
            || globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_
            || globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_
            || globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_;
    }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary control volume.
     *
     * \param globalPos The position of the center of the finite volume
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        values.setDirichlet(Indices::velocityXIdx);
        values.setDirichlet(Indices::velocityYIdx);
        values.setDirichlet(Indices::turbulentKineticEnergyIdx);
        values.setDirichlet(Indices::dissipationIdx);
        return values;
    }

    /*!
     * \brief Returns whether a fixed Dirichlet value shall be used at a given cell.
     *
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param scv The sub control volume
     * \param pvIdx The primary variable index in the solution vector
     */
    bool isDirichletCell(const Element& element,
                         const typename GridGeometry::LocalView& fvGeometry,
                         const typename GridGeometry::SubControlVolume& scv,
                         int pvIdx) const
    { return isDirichletCell_(element, fvGeometry, scv, pvIdx); }

    // Evaluate the boundary conditions for a dirichlet values at the boundary.
    PrimaryVariables dirichlet(const Element& element,
                               const SubControlVolumeFace& scvf) const
    { return analyticalModel_.analyticalSolutionAtPos(scvf.center()); }

    // Evaluate the boundary conditions for fixed values at cell centers
    PrimaryVariables dirichlet([[maybe_unused]] const Element& element,
                               const SubControlVolume& scv) const
    { return analyticalModel_.analyticalSolutionAtPos(scv.dofPosition()); }

    PrimaryVariables sourceAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables source(0.0);
        source = analyticalModel_.analyticalSourceTermAtPos(globalPos, density_, kinematicViscosity_);
        return source;
    }

    // Evaluate the initial value for a control volume.
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables initial(1.0);
        initial[Indices::velocityXIdx] = analyticalModel_.analyticalSolutionAtPos(globalPos)[Indices::velocityXIdx];
        initial[Indices::velocityYIdx] = analyticalModel_.analyticalSolutionAtPos(globalPos)[Indices::velocityYIdx];
        initial[Indices::turbulentKineticEnergyIdx] = analyticalModel_.analyticalSolutionAtPos(globalPos)[Indices::turbulentKineticEnergyIdx];
        initial[Indices::dissipationIdx] = analyticalModel_.analyticalSolutionAtPos(globalPos)[Indices::dissipationIdx];
        return initial;
    }

    auto& getAnalyticalPressureSolution() const
    { return analyticalPressureSolution_; }

    auto& getAnalyticalTKESolution() const
    { return analyticalTKESolution_; }

    auto& getAnalyticalDissipationSolution() const
    { return analyticalDissipationSolution_; }

    auto& getAnalyticalVelocitySolution() const
    { return analyticalVelocitySolution_; }

    auto& getAnalyticalVelocitySolutionOnFace() const
    { return analyticalVelocityOnFaceSolution_; }

private:

    void createManufacturedSolution_()
    {
        analyticalPressureSolution_.resize(this->gridGeometry().numCellCenterDofs());
        analyticalTKESolution_.resize(this->gridGeometry().numCellCenterDofs());
        analyticalDissipationSolution_.resize(this->gridGeometry().numCellCenterDofs());
        analyticalVelocitySolution_.resize(this->gridGeometry().numCellCenterDofs());
        analyticalVelocityOnFaceSolution_.resize(this->gridGeometry().numFaceDofs());

        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            auto fvGeometry = localView(this->gridGeometry());
            fvGeometry.bindElement(element);
            for (auto&& scv : scvs(fvGeometry))
            {
                auto ccDofIdx = scv.dofIndex();
                auto ccDofPosition = scv.dofPosition();
                auto analyticalSolutionAtCc = analyticalModel_.analyticalSolutionAtPos(ccDofPosition);
                analyticalPressureSolution_[ccDofIdx] = analyticalSolutionAtCc[Indices::pressureIdx];
                analyticalTKESolution_[ccDofIdx] = analyticalSolutionAtCc[Indices::turbulentKineticEnergyIdx];
                analyticalDissipationSolution_[ccDofIdx] = analyticalSolutionAtCc[Indices::dissipationIdx];

                for(int dirIdx = 0; dirIdx < ModelTraits::dim(); ++dirIdx)
                    analyticalVelocitySolution_[ccDofIdx][dirIdx] = analyticalSolutionAtCc[Indices::velocity(dirIdx)];

                // velocities on faces
                for (auto&& scvf : scvfs(fvGeometry))
                {
                    const auto faceDofIdx = scvf.dofIndex();
                    const auto faceDofPosition = scvf.center();
                    const auto dirIdx = scvf.directionIndex();
                    const auto analyticalSolutionAtFace = analyticalModel_.analyticalSolutionAtPos(faceDofPosition);
                    analyticalVelocityOnFaceSolution_[faceDofIdx][dirIdx] = analyticalSolutionAtFace[Indices::velocity(dirIdx)];
                }
            }
        }
    }

    //! Should the cell be dirichlet constrained for a certain eq?
    template<class Element, class FVElementGeometry, class SubControlVolume>
    bool isDirichletCell_([[maybe_unused]] const Element& element,
                          const FVElementGeometry& fvGeometry,
                          [[maybe_unused]] const SubControlVolume& scv,
                          const int& pvIdx) const
    {
        for (const auto& scvf : scvfs(fvGeometry))
        {
            if (isOnWallAtPos(scvf.center()) && (pvIdx == Indices::pressureIdx
                                              || pvIdx == Indices::turbulentKineticEnergyEqIdx
                                              || pvIdx == Indices::dissipationIdx))
                return true;
        }
        return false;
    }

    Scalar eps_;
    Scalar density_;
    Scalar kinematicViscosity_;
    AnalyticalModel analyticalModel_;

    std::vector<Scalar> analyticalPressureSolution_;
    std::vector<Scalar> analyticalTKESolution_;
    std::vector<Scalar> analyticalDissipationSolution_;
    std::vector<VelocityVector> analyticalVelocitySolution_;
    std::vector<VelocityVector> analyticalVelocityOnFaceSolution_;

    std::vector<Scalar> analyticalPressureSource_;
    std::vector<Scalar> analyticalTKESource_;
    std::vector<Scalar> analyticalDissipationSource_;
    std::vector<VelocityVector> analyticalVelocitySource_;
    std::vector<VelocityVector> analyticalVelocityOnFaceSource_;

    std::string name_;
    std::string turbulenceModelName_;
};
} // end namespace Dumux

#endif
