// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Bulk and facet problems of the complex-valued facet-coupled Helmholtz test
 *
 * The exact solution is the one of the analytical single-phase facet-coupling test, multiplied by a
 * complex amplitude \f$ \alpha \f$: \f$ u = \alpha (k_f \cos x \cosh y + (1 - k_f) \cos x \cosh(a/2)) \f$ in the bulk and
 * \f$ u_f = \alpha \cos x \f$ in the facet, with the facet coefficient \f$ k_f \f$ and the aperture \f$ a \f$. The
 * reaction term \f$ -k^2 u \f$ is compensated by the sources, so that for \f$ \alpha = 1 \f$ and \f$ k^2 = 0 \f$ the
 * discrete problem is the one of the single-phase test.
 */
#ifndef DUMUX_TEST_COMPLEX_FACET_HELMHOLTZ_PROBLEMS_HH
#define DUMUX_TEST_COMPLEX_FACET_HELMHOLTZ_PROBLEMS_HH

#include <cmath>
#include <complex>
#include <memory>
#include <string>
#include <vector>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/fvproblem.hh>

namespace Dumux {

template<class TypeTag>
class FacetHelmholtzProblemBase : public FVProblem<TypeTag>
{
    using ParentType = FVProblem<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
public:
    using Complex = std::complex<Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using GlobalPosition = typename GridGeometry::GridView::template Codim<0>::Entity::Geometry::GlobalCoordinate;
    using BoundaryTypes = Dumux::BoundaryTypes<1>;

    FacetHelmholtzProblemBase(std::shared_ptr<const GridGeometry> gridGeometry,
                              std::shared_ptr<CouplingManager> couplingManager,
                              const std::string& paramGroup)
    : ParentType(gridGeometry, paramGroup)
    , couplingManager_(couplingManager)
    {
        facetCoefficient_ = getParam<Scalar>("LowDim.SpatialParams.Permeability");
        aperture_ = getParam<Scalar>("LowDim.SpatialParams.Aperture");
        const auto alpha = getParam<std::vector<Scalar>>("Problem.Amplitude", std::vector<Scalar>{1.0, 0.0});
        const auto kSquared = getParam<std::vector<Scalar>>("Problem.WaveNumberSquared", std::vector<Scalar>{0.0, 0.0});
        amplitude_ = Complex(alpha[0], alpha[1]);
        waveNumberSquared_ = Complex(kSquared[0], kSquared[1]);
    }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllDirichlet();
        return values;
    }

    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(this->asImp_().exact(globalPos)); }

    const Complex& waveNumberSquared() const { return waveNumberSquared_; }
    const Complex& amplitude() const { return amplitude_; }
    const CouplingManager& couplingManager() const { return *couplingManager_; }

protected:
    const auto& asImp_() const { return static_cast<const GetPropType<TypeTag, Properties::Problem>&>(*this); }

    Scalar facetCoefficient_;
    Scalar aperture_;
    Complex amplitude_;
    Complex waveNumberSquared_;

private:
    std::shared_ptr<CouplingManager> couplingManager_;
};

template<class TypeTag>
class FacetHelmholtzBulkProblem : public FacetHelmholtzProblemBase<TypeTag>
{
    using ParentType = FacetHelmholtzProblemBase<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
public:
    using typename ParentType::Complex;
    using typename ParentType::NumEqVector;
    using typename ParentType::GlobalPosition;
    using typename ParentType::BoundaryTypes;
    using ParentType::ParentType;

    BoundaryTypes interiorBoundaryTypes(const Element& element, const SubControlVolumeFace& scvf) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        return values;
    }

    NumEqVector sourceAtPos(const GlobalPosition& globalPos) const
    {
        const Complex f = this->amplitude_*(1.0 - this->facetCoefficient_)*std::cos(globalPos[0])*std::cosh(0.5*this->aperture_);
        return NumEqVector(f - this->waveNumberSquared_*exact(globalPos));
    }

    Scalar exactProfile(const GlobalPosition& globalPos) const
    {
        const auto x = globalPos[0], y = globalPos[1];
        return this->facetCoefficient_*std::cos(x)*std::cosh(y) + (1.0 - this->facetCoefficient_)*std::cos(x)*std::cosh(0.5*this->aperture_);
    }

    Complex exact(const GlobalPosition& globalPos) const
    { return this->amplitude_*exactProfile(globalPos); }

    Scalar extrusionFactorAtPos(const GlobalPosition& globalPos) const { return 1.0; }
    Scalar coefficientAtPos(const GlobalPosition& globalPos) const { return 1.0; }
};

template<class TypeTag>
class FacetHelmholtzLowDimProblem : public FacetHelmholtzProblemBase<TypeTag>
{
    using ParentType = FacetHelmholtzProblemBase<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
public:
    using typename ParentType::Complex;
    using typename ParentType::NumEqVector;
    using typename ParentType::GlobalPosition;
    using ParentType::ParentType;

    //! the fluxes from the bulk domain enter as source, the reaction term is compensated
    template<class ElementVolumeVariables>
    NumEqVector source(const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolume& scv) const
    {
        auto source = this->couplingManager().evalSourcesFromBulk(element, fvGeometry, elemVolVars, scv);
        source /= scv.volume()*elemVolVars[scv].extrusionFactor();
        source[0] -= this->waveNumberSquared_*exact(scv.center());
        return source;
    }

    Scalar exactProfile(const GlobalPosition& globalPos) const
    { return std::cos(globalPos[0])*std::cosh(globalPos[1]); }

    Complex exact(const GlobalPosition& globalPos) const
    { return this->amplitude_*exactProfile(globalPos); }

    Scalar extrusionFactorAtPos(const GlobalPosition& globalPos) const { return this->aperture_; }
    Scalar coefficientAtPos(const GlobalPosition& globalPos) const { return this->facetCoefficient_; }
};

} // end namespace Dumux

#endif
