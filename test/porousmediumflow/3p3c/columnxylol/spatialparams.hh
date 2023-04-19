// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ThreePThreeCTests
 * \brief Definition of the spatial parameters for the column problem.
 */
#ifndef DUMUX_COLUMNXYLOL_SPATIAL_PARAMS_HH
#define DUMUX_COLUMNXYLOL_SPATIAL_PARAMS_HH

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/porousmediumflow/fvspatialparamsmp.hh>
#include <dumux/material/fluidmatrixinteractions/3p/parkervangenuchten.hh>

namespace Dumux {

/*!
 * \ingroup ThreePThreeCModel
 * \brief Definition of the spatial parameters for the column problem.
 */
template<class GridGeometry, class Scalar>
class ColumnSpatialParams
: public FVPorousMediumFlowSpatialParams<GridGeometry, Scalar,
                                         ColumnSpatialParams<GridGeometry, Scalar>>
{
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;

    using Element = typename GridView::template Codim<0>::Entity;
    using ParentType = FVPorousMediumFlowSpatialParams<GridGeometry, Scalar,
                                                       ColumnSpatialParams<GridGeometry, Scalar>>;

    using GlobalPosition = typename SubControlVolume::GlobalPosition;

    using ThreePhasePcKrSw = FluidMatrix::ParkerVanGenuchten3PDefault<Scalar>;

public:
    using PermeabilityType = Scalar;

    ColumnSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    , pcKrSwCurveFine_("SpatialParams.Fine")
    , pcKrSwCurveCoarse_("SpatialParams.Coarse")
    {
        // intrinsic permeabilities
        fineK_ = 1.4e-11;
        coarseK_ = 1.4e-8;

        // porosities
        porosity_ = 0.46;

        // specific heat capacities
        fineHeatCap_ = getParam<Scalar>("Component.SolidHeatCapacityFine", 850.0);
        coarseHeatCap_ = getParam<Scalar>("Component.SolidHeatCapacityCoarse", 84000.0);
    }

    /*!
     * \brief Function for defining the (intrinsic) permeability \f$[m^2]\f$
     * \note  It is possibly solution dependent.
     *
     * \param element The current element
     * \param scv The sub-control volume inside the element.
     * \param elemSol The solution at the dofs connected to the element.
     * \return The permeability
     */
    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol) const
    {
        const auto& globalPos = scv.dofPosition();
        if (isFineMaterial_(globalPos))
            return fineK_;
        return coarseK_;
    }

    /*! \brief Defines the porosity in [-].
    *
    * \param globalPos The global position where we evaluate
    */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    {
        return porosity_;
    }

    /*!
     * \brief Returns the fluid-matrix interaction law.
     *
     * \param element The current element
     * \param scv The sub-control volume inside the element.
     * \param elemSol The solution at the dofs connected to the element.
     */
     template<class ElementSolution>
     auto fluidMatrixInteraction(const Element& element,
                                 const SubControlVolume& scv,
                                 const ElementSolution& elemSol) const
    {
        const auto& globalPos = scv.dofPosition();
        if (isFineMaterial_(globalPos))
            return makeFluidMatrixInteraction(pcKrSwCurveFine_);
        else
            return makeFluidMatrixInteraction(pcKrSwCurveCoarse_);
    }

    /*!
     * \brief Returns the user-defined solid heat capacity.
     *
     * \param element The current element
     * \param scv The sub-control volume inside the element.
     * \param elemSol The solution at the dofs connected to the element.
     * \param solidState The solid state
     * \return The solid heat capacity
     */
    template <class ElementSolution, class SolidState>
    Scalar solidHeatCapacity(const Element& element,
                             const SubControlVolume& scv,
                             const ElementSolution& elemSol,
                             const SolidState& solidState) const
    {
        const auto& globalPos = scv.dofPosition();
        if (isFineMaterial_(globalPos))
            return fineHeatCap_;
        else
            return coarseHeatCap_;
    }

private:
    bool isFineMaterial_(const GlobalPosition &globalPos) const
    {
        return (0.90 - eps_ <= globalPos[1]);
    }

    Scalar fineK_;
    Scalar coarseK_;

    Scalar porosity_;

    Scalar fineHeatCap_;
    Scalar coarseHeatCap_;

    const ThreePhasePcKrSw pcKrSwCurveFine_;
    const ThreePhasePcKrSw pcKrSwCurveCoarse_;

    static constexpr Scalar eps_ = 1e-6;
};

} // end namespace Dumux

#endif
