// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_COOKS_MEMBRANE_MOMENTUM_MIXED_PROBLEM_HH
#define DUMUX_COOKS_MEMBRANE_MOMENTUM_MIXED_PROBLEM_HH

#include <vector>

#include <dune/common/fvector.hh>

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/constraintinfo.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/common/math.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/problemwithspatialparams.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/dirichletconstraints.hh>

namespace Dumux {

/*!
 * \brief Momentum subproblem of Cook's membrane with the mixed u–p formulation.
 *
 * Compressible neo-Hookean material
 * \f$ \mathbf{P} = \mu (\mathbf{F} - \mathbf{F}^{-\top}) + \frac{\lambda}{2} (J^2 - 1) \mathbf{F}^{-\top} \f$,
 * split into the part \f$ \mu (\mathbf{F} - \mathbf{F}^{-\top}) \f$ supplied here and the volumetric
 * part \f$ p_s J \mathbf{F}^{-\top} \f$ with \f$ p_s = \lambda (J^2 - 1)/(2J) \f$ from the pressure
 * subproblem. The left face \f$ x = 0 \f$ is clamped, the right face \f$ x = 48 \f$ carries an upward
 * shear traction `Problem.ShearTraction`; in 3D, the face \f$ z = 0 \f$ has \f$ u_z = 0 \f$.
 */
template<class TypeTag>
class CooksMembraneMomentumMixedProblem : public Experimental::ProblemWithSpatialParams<TypeTag>
{
    using ParentType = Experimental::ProblemWithSpatialParams<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::Experimental::BoundaryTypes<PrimaryVariables::size()>;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

    static constexpr int dim = GridView::dimension;
    static constexpr int numEq = GetPropType<TypeTag, Properties::ModelTraits>::numEq();
    using ConstraintInfo = Dumux::DirichletConstraintInfo<numEq>;
    using ConstraintValues = Dune::FieldVector<Scalar, numEq>;
    using DirichletConstraintData = Dumux::DirichletConstraintData<ConstraintInfo, ConstraintValues, GridIndexType>;

public:
    CooksMembraneMomentumMixedProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                                      std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry)
    , couplingManager_(couplingManager)
    , shearTraction_(getParam<Scalar>("Problem.ShearTraction"))
    {
        appendDirichletConstraints_();
    }

    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

    template<class Tensor, class FVElementGeometry_, class GlobalPosition_>
    Tensor isochoricFirstPiolaKirchhoffStressTensor(const Tensor& F,
                                                    const FVElementGeometry_&,
                                                    const GlobalPosition_&) const
    {
        auto FinvT = F;
        FinvT.invert();
        FinvT = transpose(FinvT);

        Tensor P = F;
        P.axpy(-1.0, FinvT);
        P *= this->spatialParams().shearModulus();
        return P;
    }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        if (globalPos[0] > 48.0 - eps_)
            values.setAllFluxBoundary();
        return values;
    }

    //! The negative traction: upward shear on the right face
    template<class ElementVariables, class IpData>
    NumEqVector boundaryFlux(const FVElementGeometry&, const ElementVariables&, const IpData& ipData) const
    {
        NumEqVector flux(0.0);
        if (ipData.global()[0] > 48.0 - eps_)
            flux[1] = -shearTraction_;
        return flux;
    }

    const auto& constraints() const
    { return constraints_; }

private:
    void appendDirichletConstraints_()
    {
        auto fvGeometry = localView(this->gridGeometry());
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            fvGeometry.bind(element);
            for (const auto& boundaryFace : boundaryFaces(fvGeometry))
            {
                const auto& center = boundaryFace.center();
                ConstraintInfo info;
                if (center[0] < eps_)
                    info.setAll();
                else if (dim == 3 && center[dim-1] < eps_)
                    info.set(dim-1, dim-1);
                else
                    continue;

                for (const auto& localDof : localDofs(fvGeometry, boundaryFace))
                    constraints_.push_back(DirichletConstraintData{info, ConstraintValues(0.0), localDof.dofIndex()});
            }
        }
    }

    static constexpr Scalar eps_ = 1e-6;
    std::shared_ptr<CouplingManager> couplingManager_;
    Scalar shearTraction_;
    std::vector<DirichletConstraintData> constraints_;
};

} // end namespace Dumux

#endif
