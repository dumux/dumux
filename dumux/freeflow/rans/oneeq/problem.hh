// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup OneEqModel
 * \copydoc Dumux::OneEqProblem
 */
#ifndef DUMUX_ONEEQ_PROBLEM_HH
#define DUMUX_ONEEQ_PROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/common/staggeredfvproblem.hh>
#include <dumux/discretization/localview.hh>
#include <dumux/discretization/staggered/elementsolution.hh>
#include <dumux/discretization/methods.hh>
#include <dumux/freeflow/rans/problem.hh>

#include "model.hh"

namespace Dumux
{

/*!
 * \ingroup OneEqModel
 * \brief One-equation turbulence problem base class.
 *
 * This implements some base functionality for one-equation Spalart-Allmaras model.
 */
template<class TypeTag>
class OneEqProblem : public RANSProblem<TypeTag>
{
    using ParentType = RANSProblem<TypeTag>;
    using Implementation = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Grid = typename GridView::Grid;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using DimVector = Dune::FieldVector<Scalar, Grid::dimension>;

    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using CellCenterPrimaryVariables = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using Indices = typename GET_PROP_TYPE(TypeTag, ModelTraits)::Indices;

public:
    //! The constructor sets the gravity, if desired by the user.
    OneEqProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
        useStoredEddyViscosity_ = getParamFromGroup<bool>(this->paramGroup(),
                                                          "RANS.UseStoredEddyViscosity", false);
    }

    /*!
     * \brief Correct size of the static (solution independent) wall variables
     */
    void updateStaticWallProperties()
    {
        ParentType::updateStaticWallProperties();

        // update size and initial values of the global vectors
        storedDynamicEddyViscosity_.resize(this->fvGridGeometry().elementMapper().size(), 0.0);
        storedViscosityTilde_.resize(this->fvGridGeometry().elementMapper().size(), 0.0);
        storedViscosityTildeGradient_.resize(this->fvGridGeometry().elementMapper().size(), DimVector(0.0));
    }

    /*!
     * \brief Update the dynamic (solution dependent) relations to the walls
     *
     * \param curSol The solution vector.
     */
    void updateDynamicWallProperties(const SolutionVector& curSol)
    {
        ParentType::updateDynamicWallProperties(curSol);

        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            unsigned int elementIdx = this->fvGridGeometry().elementMapper().index(element);

            auto fvGeometry = localView(this->fvGridGeometry());
            fvGeometry.bindElement(element);
            for (auto&& scv : scvs(fvGeometry))
            {
                const int dofIdx = scv.dofIndex();
                const auto& cellCenterPriVars = curSol[FVGridGeometry::cellCenterIdx()][dofIdx];
                PrimaryVariables priVars = makePriVarsFromCellCenterPriVars<PrimaryVariables>(cellCenterPriVars);
                auto elemSol = elementSolution<typename FVGridGeometry::LocalView>(std::move(priVars));
                // NOTE: first update the turbulence quantities
                storedViscosityTilde_[elementIdx] = elemSol[0][Indices::viscosityTildeIdx];
                // NOTE: then update the volVars
                VolumeVariables volVars;
                volVars.update(elemSol, asImp_(), element, scv);
                storedDynamicEddyViscosity_[elementIdx] = volVars.calculateEddyViscosity();
            }
        }

        // calculate cell-center-averaged velocity gradients, maximum, and minimum values
        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            unsigned int elementIdx = this->fvGridGeometry().elementMapper().index(element);

            for (unsigned int dimIdx = 0; dimIdx < Grid::dimension; ++dimIdx)
            {
                storedViscosityTildeGradient_[elementIdx][dimIdx]
                    = (storedViscosityTilde_[ParentType::neighborIdx_[elementIdx][dimIdx][1]]
                          - storedViscosityTilde_[ParentType::neighborIdx_[elementIdx][dimIdx][0]])
                      / (ParentType::cellCenter_[ParentType::neighborIdx_[elementIdx][dimIdx][1]][dimIdx]
                          - ParentType::cellCenter_[ParentType::neighborIdx_[elementIdx][dimIdx][0]][dimIdx]);
            }

            auto fvGeometry = localView(this->fvGridGeometry());
            fvGeometry.bindElement(element);
            for (auto&& scvf : scvfs(fvGeometry))
            {
                unsigned int normDim = scvf.directionIndex();
                if (scvf.boundary() && asImp_().boundaryTypes(element, scvf).isDirichlet(Indices::viscosityTildeIdx))
                {
                    // face Value
                    Scalar dirichletViscosityTilde = asImp_().dirichlet(element, scvf)[Indices::viscosityTildeIdx];

                    unsigned int neighborIdx = ParentType::neighborIdx_[elementIdx][normDim][0];
                    if (scvf.center()[normDim] < ParentType::cellCenter_[elementIdx][normDim])
                        neighborIdx = ParentType::neighborIdx_[elementIdx][normDim][1];

                    storedViscosityTildeGradient_[elementIdx][normDim]
                        = (storedViscosityTilde_[neighborIdx] - dirichletViscosityTilde)
                          / (ParentType::cellCenter_[neighborIdx][normDim] - scvf.center()[normDim]);
                }
            }
        }
    }

public:
    std::vector<Scalar> storedDynamicEddyViscosity_;
    std::vector<Scalar> storedViscosityTilde_;
    std::vector<DimVector> storedViscosityTildeGradient_;
    bool useStoredEddyViscosity_;

private:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }
};

}

#endif
