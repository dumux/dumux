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
 * \ingroup ZeroEqModel
 * \brief Zero-equation turbulence problem base class.
 */
#ifndef DUMUX_ZEROEQ_PROBLEM_HH
#define DUMUX_ZEROEQ_PROBLEM_HH

#include <string>

#include <dune/common/math.hh>
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
 * \ingroup ZeroEqModel
 * \brief Zero-equation turbulence problem base class.
 *
 * This implements some base functionality for zero-equation models
 * and a routine for the determining the eddy viscosity of the Baldwin-Lomax model.
 */
template<class TypeTag>
class RANSProblemImpl<TypeTag, TurbulenceModel::zeroeq> : public RANSProblemBase<TypeTag>
{
    using ParentType = RANSProblemBase<TypeTag>;
    using Implementation = GetPropType<TypeTag, Properties::Problem>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Grid = typename GridView::Grid;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using CellCenterPrimaryVariables = GetPropType<TypeTag, Properties::CellCenterPrimaryVariables>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    enum {
        dim = Grid::dimension,
      };
    using DimVector = Dune::FieldVector<Scalar, dim>;
    using DimMatrix = Dune::FieldMatrix<Scalar, dim, dim>;

public:
    /*!
     * \brief The constructor
     * \param gridGeometry The finite volume grid geometry
     * \param paramGroup The parameter group in which to look for runtime parameters first (default is "")
     */
    RANSProblemImpl(std::shared_ptr<const GridGeometry> gridGeometry, const std::string& paramGroup = "")
    : ParentType(gridGeometry, paramGroup)
    {}

    /*!
     * \brief Correct size of the static (solution independent) wall variables
     */
    void updateStaticWallProperties()
    {
        if (!ParentType::isFlatWallBounded())
        {
            DUNE_THROW(Dune::NotImplemented, "\n Due to grid/geometric concerns, zero-eq models should only be used for flat channel geometries. "
                                          << "\n If your geometry is a flat channel, please set the runtime parameter RANS.IsFlatWallBounded to true. \n");
        }

        ParentType::updateStaticWallProperties();

        // update size and initial values of the global vectors
        kinematicEddyViscosity_.resize(this->gridGeometry().elementMapper().size(), 0.0);
        additionalRoughnessLength_.resize(this->gridGeometry().elementMapper().size(), 0.0);
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
        ParentType::updateDynamicWallProperties(curSol);

        // correct roughness lengths if a sand grain roughness is specified
        if (hasParam("Problem.SandGrainRoughness"))
            calculateRoughnessLength_(curSol);

        // update routine for specfic models
        if (eddyViscosityModel().compare("baldwinLomax") == 0)
            updateBaldwinLomaxProperties();
    }

    /*!
     * \brief Update the relations and coefficients for the Baldwin-Lomax turbulence model
     */
    void updateBaldwinLomaxProperties()
    {
        std::vector<Scalar> kinematicEddyViscosityInner(this->gridGeometry().elementMapper().size(), 0.0);
        std::vector<Scalar> kinematicEddyViscosityOuter(this->gridGeometry().elementMapper().size(), 0.0);
        std::vector<Scalar> kinematicEddyViscosityDifference(this->gridGeometry().elementMapper().size(), 0.0);
        std::vector<Scalar> switchingPosition(this->gridGeometry().elementMapper().size(), std::numeric_limits<Scalar>::max());

        using std::abs;
        using std::exp;
        using std::min;
        using std::sqrt;
        using Dune::power;
        const Scalar aPlus = 26.0;
        const Scalar k = 0.0168;
        const Scalar cCP = 1.6;
        const Scalar cWake = 0.25;
        const Scalar cKleb = 0.3;

        std::vector<Scalar> storedFMax;
        std::vector<Scalar> storedYFMax;
        storedFMax.resize(this->gridGeometry().elementMapper().size(), 0.0);
        storedYFMax.resize(this->gridGeometry().elementMapper().size(), 0.0);

        // (1) calculate inner viscosity and Klebanoff function
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            unsigned int elementIdx = this->gridGeometry().elementMapper().index(element);
            Scalar effectiveWallDistance = asImp_().wallDistance(elementIdx) + additionalRoughnessLength(elementIdx);
            unsigned int flowDirectionAxis = this->flowDirectionAxis(elementIdx);
            unsigned int wallNormalAxis = this->wallNormalAxis(elementIdx);

            Scalar omegaAbs = abs(this->velocityGradient(elementIdx, flowDirectionAxis, wallNormalAxis)
                                  - this->velocityGradient(elementIdx, wallNormalAxis, flowDirectionAxis));
            Scalar uStar = sqrt(this->kinematicViscosity(asImp_().wallElementIndex(elementIdx))
                                * abs(this->velocityGradient(asImp_().wallElementIndex(elementIdx), flowDirectionAxis, wallNormalAxis)));
            Scalar yPlus = effectiveWallDistance * uStar / this->kinematicViscosity(elementIdx);
            Scalar mixingLength = this->karmanConstant() * effectiveWallDistance * (1.0 - exp(-yPlus / aPlus));
            kinematicEddyViscosityInner[elementIdx] = mixingLength * mixingLength * omegaAbs;

            Scalar f = effectiveWallDistance * omegaAbs * (1.0 - exp(-yPlus / aPlus));
            if (f > storedFMax[asImp_().wallElementIndex(elementIdx)])
            {
                storedFMax[asImp_().wallElementIndex(elementIdx)] = f;
                storedYFMax[asImp_().wallElementIndex(elementIdx)] = effectiveWallDistance;
            }
        }

        // (2) calculate outer viscosity
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            unsigned int elementIdx = this->gridGeometry().elementMapper().index(element);
            Scalar effectiveWallDistance = asImp_().wallDistance(elementIdx) + additionalRoughnessLength(elementIdx);

            Scalar maxVelocityNorm = 0.0;
            Scalar minVelocityNorm = 0.0;
            for (unsigned dimIdx = 0; dimIdx < dim; ++dimIdx)
            {
                maxVelocityNorm += asImp_().velocityMaximum(asImp_().wallElementIndex(elementIdx))[dimIdx]
                                   * asImp_().velocityMaximum(asImp_().wallElementIndex(elementIdx))[dimIdx];
                minVelocityNorm += asImp_().velocityMinimum(asImp_().wallElementIndex(elementIdx))[dimIdx]
                                   * asImp_().velocityMinimum(asImp_().wallElementIndex(elementIdx))[dimIdx];
            }

            Scalar deltaU = sqrt(maxVelocityNorm) - sqrt(minVelocityNorm);
            Scalar yFMax = storedYFMax[asImp_().wallElementIndex(elementIdx)];
            Scalar fMax = storedFMax[asImp_().wallElementIndex(elementIdx)];
            Scalar fWake = min(yFMax * fMax, cWake * yFMax * deltaU * deltaU / fMax);
            Scalar fKleb = 1.0 / (1.0 + 5.5 * power(cKleb * effectiveWallDistance / yFMax, 6));
            kinematicEddyViscosityOuter[elementIdx] = k * cCP * fWake * fKleb;

            kinematicEddyViscosityDifference[elementIdx]
              = kinematicEddyViscosityInner[elementIdx] - kinematicEddyViscosityOuter[elementIdx];
        }

        // (3) switching point
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            unsigned int elementIdx = this->gridGeometry().elementMapper().index(element);
            Scalar effectiveWallDistance = asImp_().wallDistance(elementIdx) + additionalRoughnessLength(elementIdx);

            // checks if sign switches, by multiplication
            Scalar check = kinematicEddyViscosityDifference[asImp_().wallElementIndex(elementIdx)] * kinematicEddyViscosityDifference[elementIdx];
            if (check < 0 // means sign has switched
                && switchingPosition[asImp_().wallElementIndex(elementIdx)] > effectiveWallDistance)
            {
                switchingPosition[asImp_().wallElementIndex(elementIdx)] = effectiveWallDistance;
            }
        }

        // (4) finally determine eddy viscosity
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            unsigned int elementIdx = this->gridGeometry().elementMapper().index(element);
            Scalar effectiveWallDistance = asImp_().wallDistance(elementIdx) + additionalRoughnessLength(elementIdx);

            kinematicEddyViscosity_[elementIdx] = kinematicEddyViscosityInner[elementIdx];
            if (effectiveWallDistance >= switchingPosition[asImp_().wallElementIndex(elementIdx)])
            {
                kinematicEddyViscosity_[elementIdx] = kinematicEddyViscosityOuter[elementIdx];
            }
        }
    }

    std::string eddyViscosityModel() const
    {
        static const std::string eddyViscosityModel = getParamFromGroup<std::string>(this->paramGroup(), "RANS.EddyViscosityModel", "vanDriest");
        return eddyViscosityModel;
    }

    int additionalRoughnessLength(const int elementIdx) const
    { return additionalRoughnessLength_[elementIdx]; }

    Scalar kinematicEddyViscosity(const int elementIdx) const
    { return kinematicEddyViscosity_[elementIdx]; }

private:

    void calculateRoughnessLength_(const SolutionVector& curSol)
    {
        bool printedRangeWarning = false;
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            static const Scalar sandGrainRoughness = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.SandGrainRoughness");
            unsigned int elementIdx = this->gridGeometry().elementMapper().index(element);

            auto fvGeometry = localView(this->gridGeometry());
            fvGeometry.bindElement(element);
            for (auto&& scv : scvs(fvGeometry))
            {
                using std::sqrt;
                using std::exp;

                const int dofIdx = scv.dofIndex();

                // construct a privars object from the cell center solution vector
                const auto& cellCenterPriVars = curSol[GridGeometry::cellCenterIdx()][dofIdx];
                PrimaryVariables priVars = makePriVarsFromCellCenterPriVars<PrimaryVariables>(cellCenterPriVars);
                auto elemSol = elementSolution<typename GridGeometry::LocalView>(std::move(priVars));

                VolumeVariables volVars;
                volVars.update(elemSol, asImp_(), element, scv);

                Scalar ksPlus = sandGrainRoughness * volVars.uStar() / volVars.kinematicViscosity();
                if (ksPlus > 0 && eddyViscosityModel().compare("baldwinLomax") == 0)
                {
                    DUNE_THROW(Dune::NotImplemented, "Roughness is not implemented for the Baldwin-Lomax model.");
                }
                if (ksPlus > 2000.)
                {
                    std::cout << "info: equivalent sand grain roughness ks+=" << ksPlus << " at " << asImp_().cellCenter(asImp_().wallElementIndex(elementIdx))
                            << " is not in the valid range (ksPlus < 2000),"
                            << " for high ksPlus values the roughness function reaches a turning point."<< std::endl;
                    DUNE_THROW(Dune::InvalidStateException, "Unphysical roughness behavior.");
                }
                else if (ksPlus > 0.0 && ksPlus < 4.535 && !printedRangeWarning)
                {
                    Dune::dinfo << "info: equivalent sand grain roughness ks+=" << ksPlus << " at " << asImp_().cellCenter(asImp_().wallElementIndex(elementIdx))
                                << " is not in the valid range (ksPlus > 4.535) and now set to 0.0"<< std::endl;
                    ksPlus = 0.0;
                    printedRangeWarning = true;
                }
                additionalRoughnessLength_[elementIdx] = 0.9 / (volVars.uStar() / volVars.kinematicViscosity())
                                                        * (sqrt(ksPlus) - ksPlus * exp(-ksPlus / 6.0));
            }
        }
    }

    std::vector<Scalar> additionalRoughnessLength_;
    std::vector<Scalar> kinematicEddyViscosity_;

    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }
};

} // end namespace Dumux

#endif
