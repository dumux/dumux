// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_WOOD_TEST_PROBLEM_HH
#define DUMUX_WOOD_TEST_PROBLEM_HH

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/fvproblemwithspatialparams.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/properties.hh>

namespace Dumux {

template<class TypeTag>
class WoodProblem : public FVProblemWithSpatialParams<TypeTag>
{
    using ParentType = FVProblemWithSpatialParams<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    static constexpr int dim = GridView::dimension;
    using Tensor = Dune::FieldMatrix<Scalar, dim, dim>;

public:
    WoodProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        yMin_ = this->gridGeometry().bBoxMin()[1];
        yMax_ = this->gridGeometry().bBoxMax()[1];
        xMin_ = this->gridGeometry().bBoxMin()[0];
        xMax_ = this->gridGeometry().bBoxMax()[0];
        xMid_ = 0.5*(xMin_ + xMax_);

        // Cyclic humidity parameters (read from input or use defaults).
        mAirMin_ = getParam<Scalar>("Problem.CyclicMoistureMin", 0.05);
        mAirMax_ = getParam<Scalar>("Problem.CyclicMoistureMax", 0.30);
        cyclePeriod_ = getParam<Scalar>("Problem.CyclePeriodSeconds", 1e6);  // very long by default (effectively static)
        useCyclic_ = getParam<bool>("Problem.UseCyclicMoisture", false);
    }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        // Symmetry axis (x = xMid): fix horizontal displacement to prevent
        // rigid-body rotation about the vertical axis.
        if (onSymmetryAxis_(globalPos))
            values.setDirichlet(Indices::displacementIdx(0));
        // Bottom-centre point: fix vertical displacement to prevent
        // rigid-body translation. Corners are free to lift -> cusping.
        if (onBottomCentre_(globalPos))
            values.setDirichlet(Indices::displacementIdx(1));
        return values;
    }

    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        // moisture stays at the initial value at the bottom (won't be used unless marked Dirichlet)
        PrimaryVariables values(0.0);
        values[Indices::moistureIdx] = this->spatialParams().initialMoisture();
        return values;
    }

    template<class ElementVolumeVariables, class ElementFluxVariablesCache>
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVariablesCache& elemFluxVarsCache,
                        const SubControlVolumeFace& scvf) const
    {
        NumEqVector flux(0.0);
        const auto pos = scvf.ipGlobal();
        if (onTop_(pos))
        {
            // boundary-layer evaporation: q_m = h_m (m - m_air), outflow positive
            const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
            const Scalar m = elemVolVars[insideScv].moistureContent();
            const auto hM = this->spatialParams().massTransferCoefficient();
            const Scalar mAir = cyclicAirMoisture_(this->time());
            flux[Indices::moistureEqIdx] = hM * (m - mAir);
        }
        return flux;
    }

    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values(0.0);
        values[Indices::moistureIdx] = this->spatialParams().initialMoisture();
        return values;
    }

    void setTime(Scalar t)
    {
        time_ = t;
    }

    Scalar time() const
    {
        return time_;
    }

private:
    //! Compute time-dependent air moisture: oscillates between mAirMin and mAirMax if useCyclic is true.
    Scalar cyclicAirMoisture_(Scalar t) const
    {
        if (!useCyclic_ || cyclePeriod_ <= 0.0)
            return this->spatialParams().airMoisture();  // use static value from params

        // Smooth triangular wave: ramps up and down over the period.
        const Scalar phase = std::fmod(t / cyclePeriod_, 1.0);
        if (phase < 0.5)
            return mAirMin_ + 2.0 * (mAirMax_ - mAirMin_) * phase;
        else
            return mAirMax_ - 2.0 * (mAirMax_ - mAirMin_) * (phase - 0.5);
    }

    bool onSymmetryAxis_(const GlobalPosition& p) const
    { return std::abs(p[0] - xMid_) < eps_*(xMax_ - xMin_); }

    bool onBottomCentre_(const GlobalPosition& p) const
    { return p[1] < yMin_ + eps_*(yMax_ - yMin_)
          && std::abs(p[0] - xMid_) < eps_*(xMax_ - xMin_); }

    bool onTop_(const GlobalPosition& p) const
    { return p[1] > yMax_ - eps_*(yMax_ - yMin_); }

    static constexpr Scalar eps_ = 1e-6;
    Scalar yMin_, yMax_, xMin_, xMax_, xMid_;
    Scalar mAirMin_, mAirMax_, cyclePeriod_;
    bool useCyclic_;

    Scalar time_ = 0.0;
};

} // end namespace Dumux

#endif
