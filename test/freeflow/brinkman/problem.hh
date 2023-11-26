// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief Darcy Brinkman model for a single-domain evaluation of coupled freeflow and porous medium flows
 */

#ifndef DUMUX_BRINKMAN_TEST_PROBLEM_HH
#define DUMUX_BRINKMAN_TEST_PROBLEM_HH

#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/freeflow/navierstokes/boundarytypes.hh>

#include <dumux/freeflow/navierstokes/momentum/fluxhelper.hh>
#include <dumux/freeflow/navierstokes/scalarfluxhelper.hh>
#include <dumux/freeflow/navierstokes/mass/1p/advectiveflux.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesTests
 * \brief  Darcy Brinkman model for a single-domain evaluation of coupled freeflow and porous medium flows
 *
 *  TODO: ADD TEXT!
 *
 */

template <class TypeTag, class BaseProblem>
class BrinkmanProblem : public BaseProblem
{
    using ParentType = BaseProblem;
    using BoundaryTypes = typename ParentType::BoundaryTypes;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using InitialValues = typename ParentType::InitialValues;
    using Sources = typename ParentType::Sources;
    using DirichletValues = typename ParentType::DirichletValues;
    using BoundaryFluxes = typename ParentType::BoundaryFluxes;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;

    static constexpr auto dimWorld = GridGeometry::GridView::dimensionworld;
    using Element = typename FVElementGeometry::Element;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    using SpatialParams = GetPropType<TypeTag, Properties::SpatialParams>;
    using PermeabilityType = typename SpatialParams::PermeabilityType;

public:
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    BrinkmanProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                    std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, couplingManager)
    {
        deltaP_ = getParam<Scalar>("Problem.PressureDifference");
        dynamicViscosity_ = getParam<Scalar>("Component.LiquidKinematicViscosity", 1.0)
                            * getParam<Scalar>("Component.LiquidDensity", 1.0);
        problemName_  =  getParam<std::string>("Problem.Name");

        storeBrinkmanParamsForOutput_();
    }


    /*!
     * \brief The problem name.
     */
    const std::string& name() const
    { return problemName_ ; }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary control volume.
     *
     * \param globalPos The position of the center of the finite volume
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;

        if constexpr (ParentType::isMomentumProblem())
        {
            if (onLowerBoundary_(globalPos) || onUpperBoundary_(globalPos))
                values.setAllDirichlet();
            else
                values.setAllNeumann();
        }
        else
            values.setNeumann(Indices::conti0EqIdx);

        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet control volume.
     *
     * \param globalPos The center of the finite volume which ought to be set.
     */
    DirichletValues dirichletAtPos(const GlobalPosition &globalPos) const
    {
        // no-flow/no-slip
        return DirichletValues(0.0);
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann control volume.
     *
     * \param element The element for which the Neumann boundary condition is set
     * \param fvGeometry The fvGeometry
     * \param elemVolVars The element volume variables
     * \param elemFaceVars The element face variables
     * \param scvf The boundary sub control volume face
     */
    template<class ElementVolumeVariables, class ElementFluxVariablesCache>
    BoundaryFluxes neumann(const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const ElementVolumeVariables& elemVolVars,
                           const ElementFluxVariablesCache& elemFluxVarsCache,
                           const SubControlVolumeFace& scvf) const
    {
        BoundaryFluxes values(0.0);
        const auto& globalPos = scvf.ipGlobal();

        if constexpr (ParentType::isMomentumProblem())
        {
            const auto p = isInlet_(globalPos) ? referencePressure(element, fvGeometry, scvf) + deltaP_
                                               : referencePressure(element, fvGeometry, scvf);
            values = NavierStokesMomentumBoundaryFluxHelper::fixedPressureMomentumFlux(
                *this, fvGeometry, scvf, elemVolVars, elemFluxVarsCache, p);
        }
        else
        {
            values = NavierStokesScalarBoundaryFluxHelper<AdvectiveFlux<ModelTraits>>::scalarOutflowFlux(
                *this, element, fvGeometry, scvf, elemVolVars);
        }

        return values;
    }


    /*!
     * \brief Returns a reference pressure at a given sub control volume face.
     *        This pressure is subtracted from the actual pressure for the momentum balance
     *        which potentially helps to improve numerical accuracy by avoiding issues related do floating point arithmetic.
     */
    Scalar referencePressure(const Element& element,
                             const FVElementGeometry& fvGeometry,
                             const SubControlVolumeFace& scvf) const
    { return 0.0; }

    const std::vector<PermeabilityType> permeabilityOutput() const
    { return outputPermeability_; }

    const std::vector<Scalar> brinkmanEpsilon() const
    { return outputBrinkmanEpsilon_; }


private:
    bool isInlet_(const GlobalPosition& globalPos) const
    { return globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_; }

    bool isOutlet_(const GlobalPosition& globalPos) const
    { return globalPos[0] > this->gridGeometry().bBoxMin()[0] + eps_; }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] < this->gridGeometry().bBoxMin()[1] + eps_; }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_; }

    void storeBrinkmanParamsForOutput_()
    {
        if constexpr (!ParentType::isMomentumProblem())
        {
            outputPermeability_.resize(this->gridGeometry().numDofs());
            outputBrinkmanEpsilon_.resize(this->gridGeometry().numDofs());
            for (const auto& element : elements(this->gridGeometry().gridView()))
            {
                const auto& eIdx = this->gridGeometry().elementMapper().index(element);
                outputPermeability_[eIdx] = this->spatialParams().permeabilityAtPos(element.geometry().center());
                outputBrinkmanEpsilon_[eIdx] = this->spatialParams().brinkmanEpsilonAtPos(element.geometry().center());
            }
        }
    }

    static constexpr Scalar eps_ = 1e-6;
    Scalar dynamicViscosity_;
    Scalar deltaP_;

    std::vector<PermeabilityType> outputPermeability_;
    std::vector<Scalar> outputBrinkmanEpsilon_;

    std::string problemName_;

};
} // end namespace Dumux
#endif
