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
 * \ingroup ZeroEqModel
 * \copydoc Dumux::ZeroEqProblem
 */
#ifndef DUMUX_ZEROEQ_PROBLEM_HH
#define DUMUX_ZEROEQ_PROBLEM_HH

#include <string>

#include <dumux/common/properties.hh>
#include <dumux/common/staggeredfvproblem.hh>
#include <dumux/discretization/localview.hh>
#include <dumux/discretization/staggered/elementsolution.hh>
#include <dumux/discretization/methods.hh>
#include <dumux/freeflow/rans/problem.hh>

#include "model.hh"

namespace Dumux {

/*!
 * \ingroup ZeroEqModel
 * \brief Zero-equation turbulence problem base class.
 *
 * This implements some base functionality for zero-equation models
 * and a routine for the determining the eddy viscosity of the Baldwin-Lomax model.
 */
template<class TypeTag>
class ZeroEqProblem : public RANSProblem<TypeTag>
{
    using ParentType = RANSProblem<TypeTag>;
    using Implementation = typename GET_PROP_TYPE(TypeTag, Problem);

    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Grid = typename GridView::Grid;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);

    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using CellCenterPrimaryVariables = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using Indices = typename GET_PROP_TYPE(TypeTag, ModelTraits)::Indices;

    enum {
        dim = Grid::dimension,
      };
    using DimVector = Dune::FieldVector<Scalar, dim>;
    using DimMatrix = Dune::FieldMatrix<Scalar, dim, dim>;

public:
    /*!
     * \brief The constructor
     * \param fvGridGeometry The finite volume grid geometry
     * \param paramGroup The parameter group in which to look for runtime parameters first (default is "")
     */
    ZeroEqProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry, const std::string& paramGroup = "")
    : ParentType(fvGridGeometry, paramGroup)
    {
        eddyViscosityModel_ = getParamFromGroup<std::string>(paramGroup, "RANS.EddyViscosityModel", "vanDriest");
    }

    /*!
     * \brief Correct size of the static (solution independent) wall variables
     */
    void updateStaticWallProperties()
    {
        ParentType::updateStaticWallProperties();

        // update size and initial values of the global vectors
        kinematicEddyViscosity_.resize(this->fvGridGeometry().elementMapper().size(), 0.0);
        additionalRoughnessLength_.resize(this->fvGridGeometry().elementMapper().size(), 0.0);
    }

    /*!
     * \brief Update the dynamic (solution dependent) relations to the walls
     *
     * This calculates the roughness related properties
     *
     * \param curSol The solution vector.
     */

    void updateDynamicWallProperties(const SolutionVector& curSol)
    {
        updateDynamicWallProperties(curSol, [](auto& x){});
    }
    /*!
     * \brief Update the dynamic (solution dependent) relations to the walls
     *
     * This calculates the roughness related properties
     *
     * \param curSol The solution vector.
     */
    template<class ContexHelper>
    void updateDynamicWallProperties(const SolutionVector& curSol, ContexHelper contextHelper)
    {
        ParentType::updateDynamicWallProperties(curSol, contextHelper);

        // calculate additional roughness
        bool printedRangeWarning = false;
        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            unsigned int elementIdx = this->fvGridGeometry().elementMapper().index(element);

            auto fvGeometry = localView(this->fvGridGeometry());
            fvGeometry.bindElement(element);
            for (auto&& scv : scvs(fvGeometry))
            {
                using std::sqrt;
                using std::exp;

                const int dofIdx = scv.dofIndex();

                // construct a privars object from the cell center solution vector
                const auto& cellCenterPriVars = curSol[FVGridGeometry::cellCenterIdx()][dofIdx];
                PrimaryVariables priVars = makePriVarsFromCellCenterPriVars<PrimaryVariables>(cellCenterPriVars);
                auto elemSol = elementSolution<typename FVGridGeometry::LocalView>(std::move(priVars));

                VolumeVariables volVars;
                volVars.update(elemSol, asImp_(), element, scv);

                Scalar ksPlus = this->sandGrainRoughness_[elementIdx] * volVars.uStar() / volVars.kinematicViscosity();
                if (ksPlus > 0 && eddyViscosityModel_.compare("baldwinLomax") == 0)
                {
                    DUNE_THROW(Dune::NotImplemented, "Roughness is not implemented for the Baldwin-Lomax model.");
                }
                if (ksPlus > 2000.)
                {
                    std::cout << "info: equivalent sand grain roughness ks+=" << ksPlus << " at " << this->cellCenter_[this->wallElementIdx_[elementIdx]]
                              << " is not in the valid range (ksPlus < 2000),"
                              << " for high ksPlus values the roughness function reaches a turning point."<< std::endl;
                    DUNE_THROW(Dune::InvalidStateException, "Unphysical roughness behavior.");
                }
                else if (ksPlus > 0.0 && ksPlus < 4.535 && !printedRangeWarning)
                {
                    Dune::dinfo << "info: equivalent sand grain roughness ks+=" << ksPlus << " at " << this->cellCenter_[this->wallElementIdx_[elementIdx]]
                                << " is not in the valid range (ksPlus > 4.535) and now set to 0.0"<< std::endl;
                    ksPlus = 0.0;
                    printedRangeWarning = true;
                }
                additionalRoughnessLength_[elementIdx] = 0.9 / (volVars.uStar() / volVars.kinematicViscosity())
                                                        * (sqrt(ksPlus) - ksPlus * exp(-ksPlus / 6.0));
            }
        }

        // update routine for specfic models
        if (eddyViscosityModel_.compare("baldwinLomax") == 0)
            updateBaldwinLomaxProperties();
    }

    /*!
     * \brief Update the relations and coefficients for the Baldwin-Lomax turbulence model
     */
    void updateBaldwinLomaxProperties()
    {
        std::vector<Scalar> kinematicEddyViscosityInner(this->fvGridGeometry().elementMapper().size(), 0.0);
        std::vector<Scalar> kinematicEddyViscosityOuter(this->fvGridGeometry().elementMapper().size(), 0.0);
        std::vector<Scalar> kinematicEddyViscosityDifference(this->fvGridGeometry().elementMapper().size(), 0.0);
        std::vector<Scalar> switchingPosition(this->fvGridGeometry().elementMapper().size(), std::numeric_limits<Scalar>::max());

        using std::abs;
        using std::exp;
        using std::min;
        using std::pow;
        using std::sqrt;
        const Scalar aPlus = 26.0;
        const Scalar k = 0.0168;
        const Scalar cCP = 1.6;
        const Scalar cWake = 0.25;
        const Scalar cKleb = 0.3;

        std::vector<Scalar> storedFMax;
        std::vector<Scalar> storedYFMax;
        storedFMax.resize(this->fvGridGeometry().elementMapper().size(), 0.0);
        storedYFMax.resize(this->fvGridGeometry().elementMapper().size(), 0.0);

        // (1) calculate inner viscosity and Klebanoff function
        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            unsigned int elementIdx = this->fvGridGeometry().elementMapper().index(element);
            unsigned int wallElementIdx = this->wallElementIdx_[elementIdx];
            Scalar wallDistance = this->wallDistance_[elementIdx] + additionalRoughnessLength_[elementIdx];
            unsigned int flowNormalAxis = this->flowNormalAxis_[elementIdx];
            unsigned int wallNormalAxis = this->wallNormalAxis_[elementIdx];

            Scalar omegaAbs = abs(this->velocityGradients_[elementIdx][flowNormalAxis][wallNormalAxis]
                                  - this->velocityGradients_[elementIdx][wallNormalAxis][flowNormalAxis]);
            Scalar uStar = sqrt(this->kinematicViscosity_[wallElementIdx]
                                * abs(this->velocityGradients_[wallElementIdx][flowNormalAxis][wallNormalAxis]));
            Scalar yPlus = wallDistance * uStar / this->kinematicViscosity_[elementIdx];
            Scalar mixingLength = this->karmanConstant() * wallDistance * (1.0 - exp(-yPlus / aPlus));
            kinematicEddyViscosityInner[elementIdx] = mixingLength * mixingLength * omegaAbs;

            Scalar f = wallDistance * omegaAbs * (1.0 - exp(-yPlus / aPlus));
            if (f > storedFMax[wallElementIdx])
            {
                storedFMax[wallElementIdx] = f;
                storedYFMax[wallElementIdx] = wallDistance;
            }
        }

        // (2) calculate outer viscosity
        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            unsigned int elementIdx = this->fvGridGeometry().elementMapper().index(element);
            unsigned int wallElementIdx = this->wallElementIdx_[elementIdx];
            Scalar wallDistance = this->wallDistance_[elementIdx] + additionalRoughnessLength_[elementIdx];

            Scalar maxVelocityNorm = 0.0;
            Scalar minVelocityNorm = 0.0;
            for (unsigned dimIdx = 0; dimIdx < dim; ++dimIdx)
            {
                maxVelocityNorm += this->velocityMaximum_[wallElementIdx][dimIdx]
                                   * this->velocityMaximum_[wallElementIdx][dimIdx];
                minVelocityNorm += this->velocityMinimum_[wallElementIdx][dimIdx]
                                   * this->velocityMinimum_[wallElementIdx][dimIdx];
            }
            Scalar deltaU = sqrt(maxVelocityNorm) - sqrt(minVelocityNorm);
            Scalar yFMax = storedYFMax[wallElementIdx];
            Scalar fMax = storedFMax[wallElementIdx];
            Scalar fWake = min(yFMax * fMax, cWake * yFMax * deltaU * deltaU / fMax);
            Scalar fKleb = 1.0 / (1.0 + 5.5 * pow(cKleb * wallDistance / yFMax, 6.0));
            kinematicEddyViscosityOuter[elementIdx] = k * cCP * fWake * fKleb;

            kinematicEddyViscosityDifference[elementIdx]
              = kinematicEddyViscosityInner[elementIdx] - kinematicEddyViscosityOuter[elementIdx];
        }

        // (3) switching point
        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            unsigned int elementIdx = this->fvGridGeometry().elementMapper().index(element);
            unsigned int wallElementIdx = this->wallElementIdx_[elementIdx];
            Scalar wallDistance = this->wallDistance_[elementIdx] + additionalRoughnessLength_[elementIdx];

            // checks if sign switches, by multiplication
            Scalar check = kinematicEddyViscosityDifference[wallElementIdx] * kinematicEddyViscosityDifference[elementIdx];
            if (check < 0 // means sign has switched
                && switchingPosition[wallElementIdx] > wallDistance)
            {
                switchingPosition[wallElementIdx] = wallDistance;
            }
        }

        // (4) finally determine eddy viscosity
        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            unsigned int elementIdx = this->fvGridGeometry().elementMapper().index(element);
            unsigned int wallElementIdx = this->wallElementIdx_[elementIdx];
            Scalar wallDistance = this->wallDistance_[elementIdx] + additionalRoughnessLength_[elementIdx];

            kinematicEddyViscosity_[elementIdx] = kinematicEddyViscosityInner[elementIdx];
            if (wallDistance >= switchingPosition[wallElementIdx])
            {
                kinematicEddyViscosity_[elementIdx] = kinematicEddyViscosityOuter[elementIdx];
            }
        }
    }

public:
    std::string eddyViscosityModel_;
    std::vector<Scalar> kinematicEddyViscosity_;
    std::vector<Scalar> additionalRoughnessLength_;

private:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }
};

} // end namespace Dumux

#endif
