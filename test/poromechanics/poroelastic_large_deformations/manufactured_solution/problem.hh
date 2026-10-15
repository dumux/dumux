// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoromechanicsTests
 * \brief Spatial parameters and subdomain problems of the plane-strain manufactured
 *        solution for the large-deformation poroelastic model.
 *
 * Displacement and fluid pressure are prescribed on the whole boundary of the unit square,
 * where both vanish. The body force and the fluid source are those of the exact solution
 * in solution.hh.
 */
#ifndef DUMUX_TEST_POROELASTIC_LARGE_DEF_MMS_PROBLEM_HH
#define DUMUX_TEST_POROELASTIC_LARGE_DEF_MMS_PROBLEM_HH

#include <cmath>
#include <memory>
#include <utility>
#include <vector>

#include <dune/common/fmatrix.hh>

#include <dumux/common/boundarytypes_.hh>
#include <dumux/common/constraintinfo.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/problemwithspatialparams.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/dirichletconstraints.hh>

#include "solution.hh"

namespace Dumux {

template<class GridDiscretization, class Scalar>
class PoroElasticLargeDefMMSSpatialParams
{
    using Element = typename GridDiscretization::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename GridDiscretization::GlobalCoordinate;

public:
    PoroElasticLargeDefMMSSpatialParams(std::shared_ptr<const GridDiscretization>)
    : G_(getParam<Scalar>("SpatialParams.ShearModulus"))
    , K_(getParam<Scalar>("SpatialParams.BulkModulus"))
    , alphaB_(getParam<Scalar>("SpatialParams.BiotCoefficient"))
    , Sp_(getParam<Scalar>("SpatialParams.StorageCoefficient"))
    , kappa0_(getParam<Scalar>("SpatialParams.ReferencePermeability"))
    , muf_(getParam<Scalar>("SpatialParams.FluidViscosity"))
    {}

    Scalar shearModulus() const { return G_; }
    Scalar bulkModulus() const { return K_; }
    Scalar biotCoefficient() const { return alphaB_; }
    Scalar storageCoefficient() const { return Sp_; }
    Scalar referencePermeability() const { return kappa0_; }
    Scalar fluidViscosity() const { return muf_; }

    template<class ElementDiscretization>
    Scalar biotCoefficient(const Element&, const ElementDiscretization&, const GlobalPosition&) const
    { return alphaB_; }

    template<class ElementDiscretization>
    Scalar storageCoefficient(const Element&, const ElementDiscretization&, const GlobalPosition&) const
    { return Sp_; }

    //! permeability growing with the pore volume, kappa = kappa_0 J
    template<class ElementDiscretization>
    Scalar permeability(const Element&, const ElementDiscretization&, const GlobalPosition&, Scalar J) const
    { return kappa0_*J; }

    template<class ElementDiscretization>
    Scalar fluidViscosity(const Element&, const ElementDiscretization&) const
    { return muf_; }

private:
    Scalar G_, K_, alphaB_, Sp_, kappa0_, muf_;
};

template<class TypeTag>
class PoroElasticLargeDefMMSMomentumProblem : public Experimental::ProblemWithSpatialParams<TypeTag>
{
    using ParentType = Experimental::ProblemWithSpatialParams<TypeTag>;
    using GridDiscretization = GetPropType<TypeTag, Properties::GridGeometry>;
    using ElementDiscretization = typename GridDiscretization::LocalView;
    using GridView = typename GridDiscretization::GridView;
    using GlobalPosition = typename GridDiscretization::GlobalCoordinate;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    static constexpr int dim = GridView::dimension;
    using Tensor = Dune::FieldMatrix<Scalar, dim, dim>;

    static constexpr int numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq();
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using ConstraintData = DirichletConstraintData<DirichletConstraintInfo<numEq>, Dune::FieldVector<Scalar, numEq>, GridIndexType>;

public:
    PoroElasticLargeDefMMSMomentumProblem(std::shared_ptr<const GridDiscretization> gridDiscretization,
                                          std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridDiscretization)
    , couplingManager_(couplingManager)
    { appendDirichletConstraints_(); }

    const CouplingManager& couplingManager() const { return *couplingManager_; }
    const std::vector<ConstraintData>& constraints() const { return constraints_; }
    void setTime(Scalar t) const { time_ = t; }

    NumEqVector sourceAtPos(const GlobalPosition& pos) const
    {
        const auto& params = this->spatialParams();
        const auto f = PoroElasticLargeDefMMS::momentumSource(pos[0], pos[1], time_,
            params.shearModulus(), params.bulkModulus(), params.biotCoefficient());
        NumEqVector source;
        source[0] = f[0];
        source[1] = f[1];
        return source;
    }

    //! Plane-strain isochoric Neo-Hookean stress, the out-of-plane stretch is one
    Tensor isochoricFirstPiolaKirchhoffStressTensor(const Tensor& F,
                                                    const ElementDiscretization&,
                                                    const GlobalPosition&) const
    {
        const Scalar G = this->spatialParams().shearModulus();
        const Scalar Jm23 = std::pow(F.determinant(), -2.0/3.0);

        Scalar I1 = 1.0;
        for (int i = 0; i < dim; ++i)
            for (int k = 0; k < dim; ++k)
                I1 += F[k][i]*F[k][i];

        Tensor Finv(F);
        Finv.invert();

        Tensor P(0.0);
        for (int i = 0; i < dim; ++i)
            for (int j = 0; j < dim; ++j)
                P[i][j] = G*Jm23*(F[i][j] - I1/3.0*Finv[j][i]);
        return P;
    }

private:
    //! The exact solution vanishes on the boundary
    void appendDirichletConstraints_()
    {
        auto elemDisc = localView(this->gridDiscretization());
        for (const auto& element : elements(this->gridDiscretization().gridView()))
        {
            elemDisc.bind(element);
            for (const auto& boundaryFace : boundaryFaces(elemDisc))
            {
                for (const auto& localDof : localDofs(elemDisc, boundaryFace))
                {
                    typename ConstraintData::ConstraintInfo info;
                    info.setAll();
                    constraints_.push_back(ConstraintData{std::move(info), typename ConstraintData::ConstraintValues(0.0), localDof.dofIndex()});
                }
            }
        }
    }

    std::shared_ptr<CouplingManager> couplingManager_;
    std::vector<ConstraintData> constraints_;
    mutable Scalar time_ = 0.0;
};

template<class TypeTag>
class PoroElasticLargeDefMMSSolidPressureProblem : public Experimental::ProblemWithSpatialParams<TypeTag>
{
    using ParentType = Experimental::ProblemWithSpatialParams<TypeTag>;
    using GridDiscretization = GetPropType<TypeTag, Properties::GridGeometry>;
    using GlobalPosition = typename GridDiscretization::GlobalCoordinate;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    using BoundaryTypes = Experimental::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;

public:
    PoroElasticLargeDefMMSSolidPressureProblem(std::shared_ptr<const GridDiscretization> gridDiscretization,
                                               std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridDiscretization)
    , couplingManager_(couplingManager)
    {}

    const CouplingManager& couplingManager() const { return *couplingManager_; }

    //! the solid pressure is algebraic and needs no boundary conditions
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition&) const
    {
        BoundaryTypes bcTypes;
        bcTypes.setAllFluxBoundary();
        return bcTypes;
    }

    Scalar volumetricPressure(Scalar J) const
    {
        const Scalar K = this->spatialParams().bulkModulus();
        return K*(J*J - 1.0)/(2.0*J);
    }

private:
    std::shared_ptr<CouplingManager> couplingManager_;
};

template<class TypeTag>
class PoroElasticLargeDefMMSFluidPressureProblem : public Experimental::ProblemWithSpatialParams<TypeTag>
{
    using ParentType = Experimental::ProblemWithSpatialParams<TypeTag>;
    using GridDiscretization = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridDiscretization::GridView;
    using GlobalPosition = typename GridDiscretization::GlobalCoordinate;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

    static constexpr int numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq();
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using ConstraintData = DirichletConstraintData<DirichletConstraintInfo<numEq>, Dune::FieldVector<Scalar, numEq>, GridIndexType>;

public:
    PoroElasticLargeDefMMSFluidPressureProblem(std::shared_ptr<const GridDiscretization> gridDiscretization,
                                               std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridDiscretization)
    , couplingManager_(couplingManager)
    { appendDirichletConstraints_(); }

    const CouplingManager& couplingManager() const { return *couplingManager_; }
    const std::vector<ConstraintData>& constraints() const { return constraints_; }
    void setTime(Scalar t) const { time_ = t; }

    NumEqVector sourceAtPos(const GlobalPosition& pos) const
    {
        const auto& params = this->spatialParams();
        return NumEqVector(PoroElasticLargeDefMMS::fluidSource(pos[0], pos[1], time_,
            params.biotCoefficient(), params.storageCoefficient(),
            params.referencePermeability(), params.fluidViscosity()));
    }

private:
    //! The exact solution vanishes on the boundary
    void appendDirichletConstraints_()
    {
        auto elemDisc = localView(this->gridDiscretization());
        for (const auto& element : elements(this->gridDiscretization().gridView()))
        {
            elemDisc.bind(element);
            for (const auto& boundaryFace : boundaryFaces(elemDisc))
            {
                for (const auto& localDof : localDofs(elemDisc, boundaryFace))
                {
                    typename ConstraintData::ConstraintInfo info;
                    info.setAll();
                    constraints_.push_back(ConstraintData{std::move(info), typename ConstraintData::ConstraintValues(0.0), localDof.dofIndex()});
                }
            }
        }
    }

    std::shared_ptr<CouplingManager> couplingManager_;
    std::vector<ConstraintData> constraints_;
    mutable Scalar time_ = 0.0;
};

} // end namespace Dumux

#endif
