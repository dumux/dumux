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
 * \ingroup ZeroEqModel
 * \brief Zero-equation turbulence problem base class.
 *
 * This implements some base functionality for zero-equation models.
 * \todo please doc me
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
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);

    enum {
        dim = Grid::dimension,
      };
    using GlobalPosition = Dune::FieldVector<Scalar, dim>;
    using DimVector = Dune::FieldVector<Scalar, dim>;
    using DimMatrix = Dune::FieldMatrix<Scalar, dim, dim>;

    enum {
        massBalanceIdx = Indices::massBalanceIdx,
        momentumBalanceIdx = Indices::momentumBalanceIdx
    };

    using DofTypeIndices = typename GET_PROP(TypeTag, DofTypeIndices);
    typename DofTypeIndices::CellCenterIdx cellCenterIdx;
    typename DofTypeIndices::FaceIdx faceIdx;

public:
    //! The constructor sets the gravity, if desired by the user.
    ZeroEqProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
        if (getParamFromGroup<bool>(GET_PROP_VALUE(TypeTag, ModelParameterGroup), "Problem.EnableGravity"))
            gravity_[dim-1]  = -9.81;

        updateStaticWallProperties();
    }

    /*!
     * \brief Update the static (solution independent) relations to the walls
     *
     * This function determines all element with a wall intersection,
     * the wall distances and the relation to the neighboring elements.
     */
    void updateStaticWallProperties() const
    {
        std::cout << "Update static wall properties2. " << this->fvGridGeometry().elementMapper().size();
        // update size and initial values of the global vectors
        kinematicEddyViscosity_.resize(this->fvGridGeometry().elementMapper().size(), 0.0);
        std::cout << " " << kinematicEddyViscosity_.size();
    }

    /*!
     * \brief Update the dynamic (solution dependent) relations to the walls
     *
     * \todo please doc me
     *
     * \param curSol The solution vector.
     */
    void updateDynamicWallProperties(const SolutionVector& curSol) const
    {
        ParentType::updateDynamicWallProperties(curSol);

        std::cout << "Update dynamic wall properties2. " << kinematicEddyViscosity_.size() << std::endl << std::flush;

//         kinematicEddyViscosity_.resize(this->fvGridGeometry().elementMapper().size(), 0.0);

        std::vector<Scalar> kinematicEddyViscosityInner(this->fvGridGeometry().elementMapper().size(), 0.0);
        std::vector<Scalar> kinematicEddyViscosityOuter(this->fvGridGeometry().elementMapper().size(), 0.0);
        std::vector<Scalar> kinematicEddyViscosityDifference(this->fvGridGeometry().elementMapper().size(), 0.0);
        std::vector<Scalar> switchingPosition(this->fvGridGeometry().elementMapper().size(), std::numeric_limits<Scalar>::max());

        static const int eddyViscosityModel
            = getParamFromGroup<int>(GET_PROP_VALUE(TypeTag, ModelParameterGroup), "RANS.EddyViscosityModel");
        if (eddyViscosityModel != Indices::baldwinLomax)
            return;

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

        static const Scalar karmanConstant
            = getParamFromGroup<Scalar>(GET_PROP_VALUE(TypeTag, ModelParameterGroup), "RANS.KarmanConstant");
        static const int flowNormalAxis
            = getParamFromGroup<int>(GET_PROP_VALUE(TypeTag, ModelParameterGroup), "RANS.FlowNormalAxis");
        static const int wallNormalAxis
            = getParamFromGroup<int>(GET_PROP_VALUE(TypeTag, ModelParameterGroup), "RANS.WallNormalAxis");

        std::vector<Scalar> storedFMax;
        std::vector<Scalar> storedYFMax;
        storedFMax.resize(this->fvGridGeometry().elementMapper().size(), 0.0);
        storedYFMax.resize(this->fvGridGeometry().elementMapper().size(), 0.0);

        // (1) calculate inner viscosity and Klebanoff function
        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            unsigned int elementID = this->fvGridGeometry().elementMapper().index(element);
            unsigned int wallElementID = this->wallElementIDs_[elementID];
            Scalar wallDistance = this->wallDistances_[elementID];

            Scalar omegaAbs = abs(this->velocityGradients_[elementID][flowNormalAxis][wallNormalAxis]
                                  - this->velocityGradients_[elementID][wallNormalAxis][flowNormalAxis]);
            Scalar uStar = sqrt(this->kinematicViscosity_[wallElementID]
                                * abs(this->velocityGradients_[wallElementID][flowNormalAxis][wallNormalAxis]));
            Scalar yPlus = wallDistance * uStar / this->kinematicViscosity_[elementID];
            Scalar mixingLength = karmanConstant * wallDistance * (1.0 - exp(-yPlus / aPlus));
            kinematicEddyViscosityInner[elementID] = mixingLength * mixingLength * omegaAbs;

            Scalar f = wallDistance * omegaAbs * (1.0 - exp(-yPlus / aPlus));
            if (f > storedFMax[wallElementID])
            {
                storedFMax[wallElementID] = f;
                storedYFMax[wallElementID] = wallDistance;
            }
        }

        // (2) calculate outer viscosity
        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            unsigned int elementID = this->fvGridGeometry().elementMapper().index(element);
            unsigned int wallElementID = this->wallElementIDs_[elementID];
            Scalar wallDistance = this->wallDistances_[elementID];

            Scalar maxVelocityNorm = 0.0;
            Scalar minVelocityNorm = 0.0;
            for (unsigned dimIdx = 0; dimIdx < dim; ++dimIdx)
            {
                maxVelocityNorm += this->velocityMaximum_[wallElementID][dimIdx]
                                   * this->velocityMaximum_[wallElementID][dimIdx];
                minVelocityNorm += this->velocityMinimum_[wallElementID][dimIdx]
                                   * this->velocityMinimum_[wallElementID][dimIdx];
            }
            Scalar deltaU = sqrt(maxVelocityNorm) - sqrt(minVelocityNorm);
            Scalar yFMax = storedYFMax[wallElementID];
            Scalar fMax = storedFMax[wallElementID];
            Scalar fWake = min(yFMax * fMax, cWake * yFMax * deltaU * deltaU / fMax);
            Scalar fKleb = 1.0 / (1.0 + 5.5 * pow(cKleb * wallDistance / yFMax, 6.0));
            kinematicEddyViscosityOuter[elementID] = k * cCP * fWake * fKleb;

            kinematicEddyViscosityDifference[elementID]
              = kinematicEddyViscosityInner[elementID] - kinematicEddyViscosityOuter[elementID];
        }

        // (3) switching point
        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            unsigned int elementID = this->fvGridGeometry().elementMapper().index(element);
            unsigned int wallElementID = this->wallElementIDs_[elementID];
            Scalar wallDistance = this->wallDistances_[elementID];

            // checks if sign switches, by multiplication
            Scalar check = kinematicEddyViscosityDifference[wallElementID] * kinematicEddyViscosityDifference[elementID];
            if (check < 0 // means sign has switched
                && switchingPosition[wallElementID] > wallDistance)
            {
                switchingPosition[wallElementID] = wallDistance;
            }
        }

        // (4) finally determine eddy viscosity
        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            unsigned int elementID = this->fvGridGeometry().elementMapper().index(element);
            unsigned int wallElementID = this->wallElementIDs_[elementID];
            Scalar wallDistance = this->wallDistances_[elementID];

            kinematicEddyViscosity_[elementID] = (wallDistance >= switchingPosition[wallElementID])
                                                 ? kinematicEddyViscosityOuter[elementID]
                                                 : kinematicEddyViscosityInner[elementID];
        }
    }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ at a given global position.
     *
     * This is not specific to the discretization. By default it just
     * calls temperature().
     *
     * \param globalPos The position in global coordinates where the temperature should be specified.
     */
    Scalar temperatureAtPos(const GlobalPosition &globalPos) const
    { return asImp_().temperature(); }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * This method MUST be overwritten by the actual problem.
     */
    Scalar temperature() const
    { DUNE_THROW(Dune::NotImplemented, "temperature() method not implemented by the actual problem"); }

    /*!
     * \brief Returns the acceleration due to gravity.
     *
     * If the <tt>Problem.EnableGravity</tt> parameter is true, this means
     * \f$\boldsymbol{g} = ( 0,\dots,\ -9.81)^T \f$, else \f$\boldsymbol{g} = ( 0,\dots, 0)^T \f$
     */
    const GlobalPosition &gravity() const
    { return gravity_; }

    //! Applys the initial face solution (velocities on the faces). Specialization for staggered grid discretization.
    template <class T = TypeTag>
    typename std::enable_if<GET_PROP_TYPE(T, FVGridGeometry)::discMethod == DiscretizationMethod::staggered, void>::type
    applyInititalFaceSolution(SolutionVector& sol,
                              const SubControlVolumeFace& scvf,
                              const PrimaryVariables& initSol) const
    {
        typename GET_PROP(TypeTag, DofTypeIndices)::FaceIdx faceIdx;
        const auto numEqCellCenter = GET_PROP_VALUE(TypeTag, NumEqCellCenter);
        sol[faceIdx][scvf.dofIndex()][numEqCellCenter] = initSol[Indices::velocity(scvf.directionIndex())];
    }

public:
    mutable std::vector<Scalar> kinematicEddyViscosity_;

private:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

    GlobalPosition gravity_;
};

}

#endif
