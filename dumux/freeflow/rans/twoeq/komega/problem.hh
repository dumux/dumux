// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup KOmegaModel
 * \brief K-Omega turbulence model problem base class.
 */
#ifndef DUMUX_KOMEGA_PROBLEM_HH
#define DUMUX_KOMEGA_PROBLEM_HH

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
 * \ingroup KOmegaModel
 * \brief K-Omega turbulence model problem base class.
 *
 * This implements the 2-equation k-omega turbulence model developed in Wilcox08 and Wilcox88
 */
template<class TypeTag>
class RANSProblemImpl<TypeTag, TurbulenceModel::komega> : public RANSProblemBase<TypeTag>
{
    using ParentType = RANSProblemBase<TypeTag>;
    using Implementation = GetPropType<TypeTag, Properties::Problem>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;

    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using CellCenterPrimaryVariables = GetPropType<TypeTag, Properties::CellCenterPrimaryVariables>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using DimVector = typename Element::Geometry::GlobalCoordinate;

public:
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
        storedDissipation_.resize(this->gridGeometry().elementMapper().size(), 0.0);
        storedDissipationGradient_.resize(this->gridGeometry().elementMapper().size(), DimVector(0.0));
        storedTurbulentKineticEnergy_.resize(this->gridGeometry().elementMapper().size(), 0.0);
        storedTurbulentKineticEnergyGradient_.resize(this->gridGeometry().elementMapper().size(), DimVector(0.0));
    }

    /*!
     * \brief Update the dynamic (solution dependent) relations to the walls
     *
     * \param curSol The solution vector.
     */
    void updateDynamicWallProperties(const SolutionVector& curSol)
    {
        ParentType::updateDynamicWallProperties(curSol);

        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            unsigned int elementIdx = this->gridGeometry().elementMapper().index(element);

            auto fvGeometry = localView(this->gridGeometry());
            fvGeometry.bindElement(element);
            for (auto&& scv : scvs(fvGeometry))
            {
                const int dofIdx = scv.dofIndex();
                const auto& cellCenterPriVars = curSol[GridGeometry::cellCenterIdx()][dofIdx];
                PrimaryVariables priVars = makePriVarsFromCellCenterPriVars<PrimaryVariables>(cellCenterPriVars);
                auto elemSol = elementSolution<typename GridGeometry::LocalView>(std::move(priVars));
                // NOTE: first update the turbulence quantities
                storedDissipation_[elementIdx] = elemSol[0][Indices::dissipationEqIdx];
                storedTurbulentKineticEnergy_[elementIdx] = elemSol[0][Indices::turbulentKineticEnergyEqIdx];
                // NOTE: then update the volVars
                VolumeVariables volVars;
                volVars.update(elemSol, asImp_(), element, scv);
                storedDynamicEddyViscosity_[elementIdx] = volVars.calculateEddyViscosity(*this);
            }
        }

        // calculate cell-centered gradients
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            const unsigned int elementIdx = this->gridGeometry().elementMapper().index(element);

            for (unsigned int dimIdx = 0; dimIdx < DimVector::dimension; ++dimIdx)
            {
                const unsigned neighborIdx0 = ParentType::neighborIndex(elementIdx, dimIdx, 0);
                const unsigned neighborIdx1 = ParentType::neighborIndex(elementIdx, dimIdx, 1);

                // Cell centered TKE Gradient
                storedTurbulentKineticEnergyGradient_[elementIdx][dimIdx]
                    = (storedTurbulentKineticEnergy(neighborIdx1) - storedTurbulentKineticEnergy(neighborIdx0))
                    / (ParentType::cellCenter(neighborIdx1)[dimIdx] - ParentType::cellCenter(neighborIdx0)[dimIdx]);
                // Cell centered Omega Gradient
                storedDissipationGradient_[elementIdx][dimIdx]
                    = (storedDissipation(neighborIdx1) - storedDissipation(neighborIdx0))
                    / (ParentType::cellCenter(neighborIdx1)[dimIdx] - ParentType::cellCenter(neighborIdx0)[dimIdx]);
            }
        }
    }

    //! \brief Returns the \f$ \beta_{\omega} \f$ constant
    const Scalar betaOmega() const
    { return 0.0708; }

    bool useStoredEddyViscosity() const
    {
        static const bool useStoredEddyViscosity = getParamFromGroup<bool>(this->paramGroup(), "RANS.UseStoredEddyViscosity", false);
        return useStoredEddyViscosity;
    }

    Scalar storedDynamicEddyViscosity(const int elementIdx) const
    { return storedDynamicEddyViscosity_[elementIdx]; }

    Scalar storedTurbulentKineticEnergy(const int elementIdx) const
    { return storedTurbulentKineticEnergy_[elementIdx]; }

    Scalar storedDissipation(const int elementIdx) const
    { return storedDissipation_[elementIdx]; }

    DimVector storedTurbulentKineticEnergyGradient(const int elementIdx) const
    { return storedTurbulentKineticEnergyGradient_[elementIdx]; }

    DimVector storedDissipationGradient(const int elementIdx) const
    { return storedDissipationGradient_[elementIdx]; }

private:
    std::vector<Scalar> storedDynamicEddyViscosity_;
    std::vector<Scalar> storedTurbulentKineticEnergy_;
    std::vector<Scalar> storedDissipation_;
    std::vector<DimVector> storedDissipationGradient_;
    std::vector<DimVector> storedTurbulentKineticEnergyGradient_;

    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }
};

} // end namespace Dumux

#endif
