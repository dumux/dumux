// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief Channel flow test for the staggered grid (Navier-)Stokes model.
 */

#ifndef DUMUX_CHANNEL_TEST_PROBLEM_HH
#define DUMUX_CHANNEL_TEST_PROBLEM_HH

#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/timeloop.hh>

#include <dumux/freeflow/navierstokes/boundarytypes.hh>

#include <dumux/freeflow/navierstokes/momentum/fluxhelper.hh>
#include <dumux/freeflow/navierstokes/scalarfluxhelper.hh>
#include <dumux/freeflow/navierstokes/mass/1p/advectiveflux.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesTests
 * \brief  Test problem for the one-phase (Navier-) Stokes problem in a channel.
 *
 * Flow from left to right in a two-dimensional channel is considered. At the inlet (left),
 * fixed values for velocity are set, while at the outlet (right), a fixed pressure
 * boundary condition is used. The channel is confined by solid walls at the top and bottom
 * of the domain which corresponds to no-slip/no-flow conditions.
 * For the non-isothermal test, water of increased temperature is injected at the inlet
 * while the walls are fully isolating.
 */
template <class TypeTag, class BaseProblem>
class ChannelTestProblem : public BaseProblem
{
    using ParentType = BaseProblem;

    using BoundaryTypes = typename ParentType::BoundaryTypes;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using InitialValues = typename ParentType::InitialValues;
    using Sources = typename ParentType::Sources;
    using DirichletValues = typename ParentType::DirichletValues;
    using BoundaryFluxes = typename ParentType::BoundaryFluxes;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;

    static constexpr auto dimWorld = GridGeometry::GridView::dimensionworld;
    using Element = typename FVElementGeometry::Element;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;


public:
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    ChannelTestProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, couplingManager)
    {
        temperatureBot_ = getParam<Scalar>("Problem.TemperatureBottom");
        temperatureTop_ = getParam<Scalar>("Problem.TemperatureTop");
        thermalExpansion_ = getParam<Scalar>("Problem.ThermalExpansion");
    }

    template<class ElementVolumeVariables>
    Sources source(const Element& element,
                   const FVElementGeometry& fvGeometry,
                   const ElementVolumeVariables& elemVolVars,
                   const SubControlVolume& scv) const
    {
        Sources source(0.0);
        if constexpr (ParentType::isMomentumProblem())
        {
            source[Indices::velocityYIdx] = this->gravity()[scv.dofAxis()] * thermalExpansion_ * (this->temperature(element, scv) - 283.15);
        }
        return source;
    }

    /*!
     * \name Boundary conditions
     */
    // \{

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
            values.setDirichlet(Indices::velocityXIdx);
            values.setDirichlet(Indices::velocityYIdx);
        }
        else
        {
            values.setAllNeumann();

            if (isPlate_(globalPos))
            {
                values.setDirichlet(Indices::pressureIdx);
#if NONISOTHERMAL
                values.setDirichlet(Indices::temperatureIdx);
#endif
            }
        }

        return values;
    }


    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet control volume.
     *
     * \param globalPos The center of the finite volume which ought to be set.
     */
    DirichletValues dirichlet(const Element& element, const SubControlVolumeFace& scvf) const
    {
        const auto& globalPos = scvf.ipGlobal();
        DirichletValues values = initialAtPos(globalPos);

        if constexpr (ParentType::isMomentumProblem())
        {
            values[Indices::velocityXIdx] = 0.0;
            values[Indices::velocityXIdx] = 0.0;
        }
        else
        {
            values[Indices::pressureIdx] = this->couplingManager().cellPressure(element, scvf);

#if NONISOTHERMAL
            if (scvf.ipGlobal()[1] < this->gridGeometry().bBoxMin()[1] + eps_)
                values[Indices::temperatureIdx] = temperatureBot_;
            if (scvf.ipGlobal()[1] > this->gridGeometry().bBoxMax()[1] - eps_)
                values[Indices::temperatureIdx] = temperatureTop_;
#endif
        }

        return values;
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

        if constexpr (ParentType::isMomentumProblem())
        {
        }
        else
        {
        }

        return values;
    }

    /*!
     * \brief Return the analytical solution of the problem at a given position
     *
     * \param globalPos The global position
     */
    DirichletValues analyticalSolution(const GlobalPosition& globalPos, Scalar time = 0.0) const
    {
        DirichletValues values(0.0);

        if constexpr (ParentType::isMomentumProblem())
        {
        }
        else
        {
        }

        return values;
    }

    // \}

   /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluates the initial value for a control volume.
     *
     * \param globalPos The global position
     */
    InitialValues initialAtPos(const GlobalPosition& globalPos) const
    {
        InitialValues values;

        if constexpr (ParentType::isMomentumProblem())
        {
            values[Indices::velocityYIdx] = 0.0;
            values[Indices::velocityXIdx] = 0.0;
        }
        else
        {
            values[Indices::pressureIdx] = 1.0e5;
#if NONISOTHERMAL
            //values[Indices::temperatureIdx] = temperatureBot_ + (temperatureTop_ - temperatureBot_)*globalPos[1]/this->gridGeometry().bBoxMax()[1];
            values[Indices::temperatureIdx] = (temperatureBot_ + temperatureTop_)/2.0;
#endif
        }



        return values;
    }

    // \}

    /*!
     * \brief Returns a reference pressure at a given sub control volume face.
     *        This pressure is subtracted from the actual pressure for the momentum balance
     *        which potentially helps to improve numerical accuracy by avoiding issues related do floating point arithmetic.
     */
    Scalar referencePressure(const Element& element,
                             const FVElementGeometry& fvGeometry,
                             const SubControlVolumeFace& scvf) const
    { return 1.0e5; }

    void setTime(Scalar time)
    {
        time_ = time;
    }

    Scalar time() const
    {
        return time_;
    }

    bool hasAnalyticalSolution() const
    { return false; }

private:
    bool isPlate_(const GlobalPosition& globalPos) const
    {
        return(globalPos[1] < eps_ || globalPos[1] > this->gridGeometry().bBoxMax()[0] - eps_);
    }

    static constexpr Scalar eps_=1e-6;
    Scalar time_;
    Scalar temperatureBot_;
    Scalar temperatureTop_;
    Scalar thermalExpansion_;
};
} // end namespace Dumux

#endif
