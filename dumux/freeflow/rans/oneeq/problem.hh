// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OneEqModel
 * \brief One-equation turbulence problem base class.
 */
#ifndef DUMUX_ONEEQ_PROBLEM_HH
#define DUMUX_ONEEQ_PROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/common/staggeredfvproblem.hh>
#include <dumux/discretization/localview.hh>
#include <dumux/discretization/staggered/elementsolution.hh>
#include <dumux/discretization/method.hh>
#include <dumux/freeflow/rans/problem.hh>
#include <dumux/freeflow/turbulencemodel.hh>

#include "model.hh"

namespace Dumux {

/*!
 * \ingroup OneEqModel
 * \brief One-equation turbulence problem base class.
 *
 * This implements some base functionality for one-equation Spalart-Allmaras model.
 */
template<class TypeTag>
class RANSProblemImpl<TypeTag, TurbulenceModel::oneeq> : public RANSProblemBase<TypeTag>
{

    using ParentType = RANSProblemBase<TypeTag>;
    using Implementation = GetPropType<TypeTag, Properties::Problem>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Grid = typename GridView::Grid;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using DimVector = Dune::FieldVector<Scalar, Grid::dimension>;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using CellCenterPrimaryVariables = GetPropType<TypeTag, Properties::CellCenterPrimaryVariables>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

public:
    //! The constructor sets the gravity, if desired by the user.
    RANSProblemImpl(std::shared_ptr<const GridGeometry> gridGeometry, const std::string& paramGroup = "")
    : ParentType(gridGeometry, paramGroup)
    { }

    /*!
     * \brief Correct size of the static (solution independent) wall variables
     */
    void updateStaticWallProperties()
    {
        ParentType::updateStaticWallProperties();

        // update size and initial values of the global vectors
        storedDynamicEddyViscosity_.resize(this->gridGeometry().elementMapper().size(), 0.0);
        storedViscosityTilde_.resize(this->gridGeometry().elementMapper().size(), 0.0);
        storedViscosityTildeGradient_.resize(this->gridGeometry().elementMapper().size(), DimVector(0.0));
    }

    /*!
     * \brief Update the dynamic (solution dependent) relations to the walls
     *
     * \param curSol The solution vector.
     */
    template<class SolutionVector>
    void updateDynamicWallProperties(const SolutionVector& curSol)
    {
        ParentType::updateDynamicWallProperties(curSol);

        auto fvGeometry = localView(this->gridGeometry());
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            unsigned int elementIdx = this->gridGeometry().elementMapper().index(element);
            fvGeometry.bindElement(element);

            for (auto&& scv : scvs(fvGeometry))
            {
                const int dofIdx = scv.dofIndex();
                const auto& cellCenterPriVars = curSol[GridGeometry::cellCenterIdx()][dofIdx];
                PrimaryVariables priVars = makePriVarsFromCellCenterPriVars<PrimaryVariables>(cellCenterPriVars);
                auto elemSol = elementSolution<typename GridGeometry::LocalView>(std::move(priVars));
                // NOTE: first update the turbulence quantities
                storedViscosityTilde_[elementIdx] = elemSol[0][Indices::viscosityTildeIdx];
                // NOTE: then update the volVars
                VolumeVariables volVars;
                volVars.update(elemSol, asImp_(), element, scv);
                storedDynamicEddyViscosity_[elementIdx] = volVars.calculateEddyViscosity();
            }
        }

        // calculate cell-center-averaged gradient
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            const unsigned int elementIdx = this->gridGeometry().elementMapper().index(element);
            fvGeometry.bindElement(element);

            for (unsigned int axisIdx = 0; axisIdx < Grid::dimension; ++axisIdx)
            {
                const unsigned int neighborIndex0 = ParentType::neighborIndex(elementIdx, axisIdx, 0);
                const unsigned int neighborIndex1 = ParentType::neighborIndex(elementIdx, axisIdx, 1);

                // calculate cell-centered turbulentEddyViscosity (viscosityTilde) gradient
                storedViscosityTildeGradient_[elementIdx][axisIdx]
                    = (storedViscosityTilde(neighborIndex1) - storedViscosityTilde(neighborIndex0))
                    / (ParentType::cellCenter(neighborIndex1)[axisIdx] - ParentType::cellCenter(neighborIndex0)[axisIdx]);
            }

            // Adjust for dirichlet boundary conditions
            for (auto&& scvf : scvfs(fvGeometry))
            {
                const unsigned int normDim = scvf.directionIndex();
                if (scvf.boundary() && asImp_().boundaryTypes(element, scvf).isDirichlet(Indices::viscosityTildeIdx))
                {
                    // face Value
                    Scalar dirichletViscosityTilde = asImp_().dirichlet(element, scvf)[Indices::viscosityTildeIdx];

                    unsigned int neighborIndex = ParentType::neighborIndex(elementIdx, normDim, 0);
                    if (scvf.center()[normDim] < ParentType::cellCenter(elementIdx)[normDim])
                        neighborIndex = ParentType::neighborIndex(elementIdx, normDim, 1);

                    storedViscosityTildeGradient_[elementIdx][normDim]
                        = (storedViscosityTilde(neighborIndex) - dirichletViscosityTilde)
                        / (ParentType::cellCenter(neighborIndex)[normDim] - scvf.center()[normDim]);
                }
            }
        }
    }

    bool useStoredEddyViscosity() const
    {
        static const bool useStoredEddyViscosity = getParamFromGroup<bool>(this->paramGroup(), "RANS.UseStoredEddyViscosity", false);
        return useStoredEddyViscosity;
    }

    Scalar storedDynamicEddyViscosity(const int elementIdx) const
    { return storedDynamicEddyViscosity_[elementIdx]; }

    Scalar storedViscosityTilde(const int elementIdx) const
    { return storedViscosityTilde_[elementIdx]; }

    DimVector storedViscosityTildeGradient(const int elementIdx) const
    { return storedViscosityTildeGradient_[elementIdx]; }

private:
    std::vector<Scalar> storedDynamicEddyViscosity_;
    std::vector<Scalar> storedViscosityTilde_;
    std::vector<DimVector> storedViscosityTildeGradient_;

    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }
};

} // end namespace Dumux

#endif
