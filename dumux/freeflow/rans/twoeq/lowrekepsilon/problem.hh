// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup LowReKEpsilonModel
 * \brief Low-Re k-epsilon turbulence problem base class.
 */
#ifndef DUMUX_LOWREKEPSILON_PROBLEM_HH
#define DUMUX_LOWREKEPSILON_PROBLEM_HH

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
 * \ingroup LowReKEpsilonModel
 * \brief Low-Re k-epsilon turbulence problem base class.
 *
 * This implements some base functionality for low-Re k-epsilon models.
 */
template<class TypeTag>
class RANSProblemImpl<TypeTag, TurbulenceModel::lowrekepsilon> : public RANSProblemBase<TypeTag>
{
    using ParentType = RANSProblemBase<TypeTag>;
    using Implementation = GetPropType<TypeTag, Properties::Problem>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;

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
        storedDissipationTilde_.resize(this->gridGeometry().elementMapper().size(), 0.0);
        storedDynamicEddyViscosity_.resize(this->gridGeometry().elementMapper().size(), 0.0);
        storedTurbulentKineticEnergy_.resize(this->gridGeometry().elementMapper().size(), 0.0);
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
                storedDissipationTilde_[elementIdx] = elemSol[0][Indices::dissipationEqIdx];
                storedTurbulentKineticEnergy_[elementIdx] = elemSol[0][Indices::turbulentKineticEnergyEqIdx];
                // NOTE: then update the volVars
                VolumeVariables volVars;
                volVars.update(elemSol, asImp_(), element, scv);
                storedDynamicEddyViscosity_[elementIdx] = volVars.calculateEddyViscosity();
            }
        }
    }

    bool useStoredEddyViscosity() const
    {
        static const bool useStoredEddyViscosity = getParamFromGroup<bool>(this->paramGroup(), "RANS.UseStoredEddyViscosity", true);
        return useStoredEddyViscosity;
    }

    Scalar storedDissipationTilde(const int elementIdx) const
    { return storedDissipationTilde_[elementIdx]; }

    Scalar storedDynamicEddyViscosity(const int elementIdx) const
    { return storedDynamicEddyViscosity_[elementIdx]; }

    Scalar storedTurbulentKineticEnergy(const int elementIdx) const
    { return storedTurbulentKineticEnergy_[elementIdx]; }

private:
    std::vector<Scalar> storedDissipationTilde_;
    std::vector<Scalar> storedDynamicEddyViscosity_;
    std::vector<Scalar> storedTurbulentKineticEnergy_;

    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }
};

}

#endif
