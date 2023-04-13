// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Linear
 * \brief Scalar Helmholtz operator to be used as a solver component
 */
#ifndef DUMUX_LINEAR_HELMHOLTZ_OPERATOR_HH
#define DUMUX_LINEAR_HELMHOLTZ_OPERATOR_HH

#include <dumux/common/math.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/discretization/cellcentered/tpfa/computetransmissibility.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/fvproblem.hh>
#include <dumux/assembly/fvassembler.hh>

namespace Dumux::Detail::HelmholtzOperator {

template <class Traits>
class HelmholtzModelVolumeVariables
{
    using Scalar = typename Traits::PrimaryVariables::value_type;
public:
    using PrimaryVariables = typename Traits::PrimaryVariables;

    template<class ElementSolution, class Problem, class Element, class SubControlVolume>
    void update(const ElementSolution& elemSol, const Problem&, const Element&, const SubControlVolume& scv)
    { priVars_ = elemSol[scv.indexInElement()]; }

    Scalar priVar(const int pvIdx) const { return priVars_[pvIdx]; }
    const PrimaryVariables& priVars() const { return priVars_; }
    Scalar extrusionFactor() const { return 1.0; }

private:
    PrimaryVariables priVars_;
};

template<class TypeTag>
class HelmholtzModelLocalResidual
: public GetPropType<TypeTag, Properties::BaseLocalResidual>
{
    using ParentType = GetPropType<TypeTag, Properties::BaseLocalResidual>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;
    using VolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::VolumeVariables;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using ElementFluxVariablesCache = typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    static constexpr int dimWorld = GridView::dimensionworld;
public:
    using ParentType::ParentType;

    NumEqVector computeStorage(const Problem&,
                               const SubControlVolume&,
                               const VolumeVariables&) const
    { return NumEqVector(0.0); }

    NumEqVector computeFlux(const Problem& problem,
                            const Element&, const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolumeFace& scvf,
                            const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        NumEqVector flux(0.0);

        // CVFE schemes
        if constexpr (GridGeometry::discMethod == DiscretizationMethods::box
            || GridGeometry::discMethod == DiscretizationMethods::fcdiamond
            || GridGeometry::discMethod == DiscretizationMethods::pq1bubble)
        {
            const auto& fluxVarCache = elemFluxVarsCache[scvf];
            Dune::FieldVector<Scalar, dimWorld> gradConcentration(0.0);
            for (const auto& scv : scvs(fvGeometry))
            {
                const auto& volVars = elemVolVars[scv];
                gradConcentration.axpy(volVars.concentration(), fluxVarCache.gradN(scv.indexInElement()));
            }

            NumEqVector flux;
            flux[0] = -1.0*vtmv(scvf.unitOuterNormal(), problem.a(), gradConcentration)*scvf.area();
        }

        // CCTpfa
        else if constexpr (GridGeometry::discMethod == DiscretizationMethods::cctpfa)
        {
            const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
            const auto& insideV = elemVolVars[insideScv];
            const auto& outsideV = elemVolVars[scvf.outsideScvIdx()];
            const auto ti = computeTpfaTransmissibility(
                fvGeometry, scvf, insideScv, problem.a(), insideV.extrusionFactor()
            );

            const auto tij = [&]{
                if (scvf.boundary())
                    return scvf.area()*ti;
                else
                {
                    const auto& outsideScv = fvGeometry.scv(scvf.outsideScvIdx());
                    const Scalar tj = -computeTpfaTransmissibility(
                        fvGeometry, scvf, outsideScv, problem.a(), outsideV.extrusionFactor()
                    );

                    return scvf.area()*(ti * tj)/(ti + tj);
                }
            }();

            const auto u = insideV.priVar(0);
            const auto v = outsideV.priVar(0);
            flux[0] = tij*(u - v);
        }
        else
            DUNE_THROW(Dune::NotImplemented, "Discretization method " << GridGeometry::discMethod);

        return flux;
    }

    NumEqVector computeSource(const Problem& problem,
                              const Element&, const FVElementGeometry&,
                              const ElementVolumeVariables& elemVolVars,
                              const SubControlVolume& scv) const
    {
        NumEqVector source(0.0);
        source[0] = -problem.b()*elemVolVars[scv].priVar(0);
        return source;
    }
};

template<class TypeTag>
class HelmholtzModelHomogeneousNeumannProblem : public FVProblem<TypeTag>
{
    using ParentType = FVProblem<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
public:
    HelmholtzModelHomogeneousNeumannProblem(std::shared_ptr<const GridGeometry> gridGeometry, Scalar a, Scalar b)
    : ParentType(gridGeometry)
    , a_(a), b_(b)
    {}

    template<class GlobalPosition>
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    { BoundaryTypes values; values.setAllNeumann(); return values; }

    Scalar a() const { return a_; }
    Scalar b() const { return b_; }
private:
    Scalar a_, b_;
};

template<class G, class D, class S>
struct Policy
{
    using Scalar = S;
    using Grid = G;
    using Discretization = D;
};

template<class P>
struct TTag {
    using InheritsFrom = std::tuple<typename P::Discretization>;
};

} // end namespace Dumux::Detail::HelmholtzOperator

namespace Dumux::Properties {

//! Set the default type of scalar values to double
template<class TypeTag, class P>
struct Scalar<TypeTag, Dumux::Detail::HelmholtzOperator::TTag<P>>
{ using type = typename P::Scalar; };

//! Set the default primary variable vector to a vector of size of number of equations
template<class TypeTag, class P>
struct PrimaryVariables<TypeTag, Dumux::Detail::HelmholtzOperator::TTag<P>>
{
    using type = Dune::FieldVector<
        GetPropType<TypeTag, Properties::Scalar>,
        GetPropType<TypeTag, Properties::ModelTraits>::numEq()
    >;
};

template<class TypeTag, class P>
struct ModelTraits<TypeTag, Dumux::Detail::HelmholtzOperator::TTag<P>>
{
    struct Traits { static constexpr int numEq() { return 1; } };
    using type = Traits;
};

template<class TypeTag, class P>
struct LocalResidual<TypeTag, Dumux::Detail::HelmholtzOperator::TTag<P>>
{ using type = Dumux::Detail::HelmholtzOperator::HelmholtzModelLocalResidual<TypeTag>; };

template<class TypeTag, class P>
struct Problem<TypeTag, Dumux::Detail::HelmholtzOperator::TTag<P>>
{ using type = Dumux::Detail::HelmholtzOperator::HelmholtzModelHomogeneousNeumannProblem<TypeTag>; };

template<class TypeTag, class P>
struct VolumeVariables<TypeTag, Dumux::Detail::HelmholtzOperator::TTag<P>>
{
    struct Traits { using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>; };
    using type = Dumux::Detail::HelmholtzOperator::HelmholtzModelVolumeVariables<Traits>;
};

template<class TypeTag, class P>
struct EnableGridVolumeVariablesCache<TypeTag, Dumux::Detail::HelmholtzOperator::TTag<P>>
{ static constexpr bool value = true; };

template<class TypeTag, class P>
struct EnableGridFluxVariablesCache<TypeTag, Dumux::Detail::HelmholtzOperator::TTag<P>>
{ static constexpr bool value = true; };

template<class TypeTag, class P>
struct EnableGridGeometryCache<TypeTag, Dumux::Detail::HelmholtzOperator::TTag<P>>
{ static constexpr bool value = true; };

template<class TypeTag, class P>
struct Grid<TypeTag, Dumux::Detail::HelmholtzOperator::TTag<P>>
{ using type = typename P::Grid; };

} // end namespace Dumux::Properties

namespace Dumux {

/*!
 * \brief make a Helmholtz matrix operator (aΔ + bI)
 * \tparam Discretization the discretization model type tag
 * \param gridView the gridView on which to assemble the discrete operator
 * \param a the weight for the Laplacian
 * \param b the weight for the mass matrix
 *
 * For example: a=-1.0 and b=0.0 will yield a negative Laplacian operator
 * For example: a=0.0 and b=1.0 will yield the mass matrix
 * For example: a=1.0 and b=1.0 will yield the unweighted Helmholtz operator
 *
 * Useful as a preconditioner component or for testing solvers
 */
template<class Discretization, class GridView, class Scalar>
auto makeHelmholtzMatrix(const GridView& gridView, const Scalar a = 1.0, const Scalar b = 1.0)
{
    using TypeTag = Detail::HelmholtzOperator::TTag<
        Detail::HelmholtzOperator::Policy<typename GridView::Grid, Discretization, Scalar>
    >;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(gridView);

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry, a, b);

    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    SolutionVector x(gridGeometry->numDofs());
    gridVariables->init(x);

    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables);
    using A = typename Assembler::JacobianMatrix;
    using V = typename Assembler::ResidualType;
    auto jacobian = std::make_shared<A>();
    auto residual = std::make_shared<V>();
    assembler->setLinearSystem(jacobian, residual);
    assembler->assembleJacobian(x);
    return jacobian;
};

/*!
 * \brief make a Laplace matrix operator aΔ
 */
template<class Discretization, class GridView, class Scalar>
auto makeLaplaceMatrix(const GridView& gridView, const Scalar a = 1.0)
{ return createHelmholtzMatrix(gridView, a, 0.0); };

template<class LinearOperator, class Discretization, class GridView, class Scalar>
std::shared_ptr<LinearOperator> makeHelmholtzLinearOperator(const GridView& gridView, const Scalar a, const Scalar b)
{ return std::make_shared<LinearOperator>(makeHelmholtzMatrix<Discretization>(gridView, a, b)); }

template<class LinearOperator, class Discretization, class GridView, class Comm, class Scalar>
std::shared_ptr<LinearOperator> makeHelmholtzLinearOperator(const GridView& gridView, const Comm& comm, const Scalar a, const Scalar b)
{ return std::make_shared<LinearOperator>(makeHelmholtzMatrix<Discretization>(gridView, a, b), comm); }

template<class LinearOperator, class Discretization, class GridView, class Scalar>
std::shared_ptr<LinearOperator> makeLaplaceLinearOperator(const GridView& gridView, const Scalar a)
{ return std::make_shared<LinearOperator>(makeLaplaceMatrix<Discretization>(gridView, a)); }

template<class LinearOperator, class Discretization, class GridView, class Comm, class Scalar>
std::shared_ptr<LinearOperator> makeLaplaceLinearOperator(const GridView& gridView, const Comm& comm, const Scalar a)
{ return std::make_shared<LinearOperator>(makeLaplaceMatrix<Discretization>(gridView, a), comm); }

} // end namespace Dumux

#endif
