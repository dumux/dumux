// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief A complex-valued Helmholtz model for the box scheme that can be used in the bulk and in the
 *        lower-dimensional domain of a facet-coupled problem
 *
 * Solves \f$ -\nabla\cdot(a \nabla u) - k^2 u = f \f$ with complex \f$ u, k^2, f \f$. The extrusion factor and
 * the coefficient \f$ a \f$ come from the problem. On interior boundaries of the facet coupling the flux
 * uses the lower-dimensional solution in the same way as the box facet-coupling Darcy law.
 */
#ifndef DUMUX_TEST_COMPLEX_FACET_HELMHOLTZ_MODEL_HH
#define DUMUX_TEST_COMPLEX_FACET_HELMHOLTZ_MODEL_HH

#include <complex>

#include <dune/common/fvector.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/properties/model.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/discretization/defaultlocaloperator.hh>

namespace Dumux {

template<class TypeTag>
class FacetHelmholtzLocalResidual : public DiscretizationDefaultLocalOperator<TypeTag>
{
    using ParentType = DiscretizationDefaultLocalOperator<TypeTag>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using PrimaryVariable = typename PrimaryVariables::value_type;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using VolumeVariables = typename GridVariables::GridVolumeVariables::VolumeVariables;
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;
    using ElementFluxVariablesCache = typename GridVariables::GridFluxVariablesCache::LocalView;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Extrusion = Extrusion_t<GridGeometry>;
    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;

    static_assert(DiscretizationMethods::isCVFE<typename GridGeometry::DiscretizationMethod>,
                  "This local residual is implemented for control-volume finite element schemes");

public:
    using ParentType::ParentType;

    NumEqVector computeStorage(const Problem& problem,
                               const SubControlVolume& scv,
                               const VolumeVariables& volVars) const
    { return NumEqVector(0.0); }

    NumEqVector computeFlux(const Problem& problem,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolumeFace& scvf,
                            const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        if constexpr (requires { scvf.interiorBoundary(); })
            if (scvf.interiorBoundary())
                return interiorBoundaryFlux_(problem, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);

        const auto& fluxVarCache = elemFluxVarsCache[scvf];
        Dune::FieldVector<PrimaryVariable, dimWorld> gradU(0.0);
        for (const auto& scv : scvs(fvGeometry))
            gradU.axpy(elemVolVars[scv].priVar(0), fluxVarCache.gradN(scv.indexInElement()));

        const auto& insideVolVars = elemVolVars[fvGeometry.scv(scvf.insideScvIdx())];
        NumEqVector flux(0.0);
        flux[0] = -insideVolVars.coefficient()*(gradU*scvf.unitOuterNormal())
                  *Extrusion::area(fvGeometry, scvf)*insideVolVars.extrusionFactor();
        return flux;
    }

    NumEqVector computeSource(const Problem& problem,
                              const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const ElementVolumeVariables& elemVolVars,
                              const SubControlVolume& scv) const
    {
        NumEqVector source = problem.source(element, fvGeometry, elemVolVars, scv);
        source[0] += problem.waveNumberSquared()*elemVolVars[scv].priVar(0);
        return source;
    }

private:
    //! flux into the facet: the normal gradient between the bulk value at the face and the facet value over half the aperture
    NumEqVector interiorBoundaryFlux_(const Problem& problem,
                                      const Element& element,
                                      const FVElementGeometry& fvGeometry,
                                      const ElementVolumeVariables& elemVolVars,
                                      const SubControlVolumeFace& scvf,
                                      const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        static_assert(dim == dimWorld, "The interior boundary flux is implemented for dim == dimWorld");

        const auto& shapeValues = elemFluxVarsCache[scvf].shapeValues();
        PrimaryVariable u(0.0);
        for (const auto& scv : scvs(fvGeometry))
            u += elemVolVars[scv].priVar(0)*shapeValues[scv.indexInElement()][0];

        const auto& facetVolVars = problem.couplingManager().getLowDimVolVars(element, scvf);
        const auto aperture = facetVolVars.extrusionFactor();
        const auto normalGradient = (facetVolVars.priVar(0) - u)/(0.5*aperture);

        const auto& insideVolVars = elemVolVars[fvGeometry.scv(scvf.insideScvIdx())];
        NumEqVector flux(0.0);
        flux[0] = -facetVolVars.coefficient()*normalGradient
                  *Extrusion::area(fvGeometry, scvf)*insideVolVars.extrusionFactor();
        return flux;
    }
};

//! primary variables, extrusion factor and coefficient of the gradient term
template<class PV>
class FacetHelmholtzVolumeVariables
{
    using Scalar = typename Dune::FieldTraits<typename PV::value_type>::real_type;
public:
    using PrimaryVariables = PV;

    template<class ElementSolution, class Problem, class Element, class SubControlVolume>
    void update(const ElementSolution& elemSol, const Problem& problem, const Element& element, const SubControlVolume& scv)
    {
        priVars_ = elemSol[scv.indexInElement()];
        extrusionFactor_ = problem.extrusionFactorAtPos(scv.center());
        coefficient_ = problem.coefficientAtPos(scv.center());
    }

    typename PV::value_type priVar(int pvIdx) const { return priVars_[pvIdx]; }
    const PrimaryVariables& priVars() const { return priVars_; }
    Scalar extrusionFactor() const { return extrusionFactor_; }
    Scalar coefficient() const { return coefficient_; }

private:
    PrimaryVariables priVars_;
    Scalar extrusionFactor_;
    Scalar coefficient_;
};

struct FacetHelmholtzModelTraits
{
    struct Indices { static constexpr int uIdx = 0; };
    static constexpr int numEq() { return 1; }
};

} // end namespace Dumux

namespace Dumux::Properties::TTag {
struct FacetHelmholtzModel { using InheritsFrom = std::tuple<ModelProperties>; };
} // end namespace Dumux::Properties::TTag

namespace Dumux::Properties {

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::FacetHelmholtzModel>
{ using type = FacetHelmholtzLocalResidual<TypeTag>; };

template<class TypeTag>
struct ModelTraits<TypeTag, TTag::FacetHelmholtzModel>
{ using type = FacetHelmholtzModelTraits; };

template<class TypeTag>
struct PrimaryVariables<TypeTag, TTag::FacetHelmholtzModel>
{ using type = Dune::FieldVector<std::complex<GetPropType<TypeTag, Properties::Scalar>>, 1>; };

template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::FacetHelmholtzModel>
{ using type = FacetHelmholtzVolumeVariables<GetPropType<TypeTag, Properties::PrimaryVariables>>; };

} // end namespace Dumux::Properties

#endif
