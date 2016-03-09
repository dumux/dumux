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
 * \brief This file contains the data which is required to calculate
 *        the fluxes of the ZeroEq model over a face of a finite volume.
 *
 * This means the methods to calculate the eddy viscosity.
 */
#ifndef DUMUX_ZEROEQ_FLUX_VARIABLES_HH
#define DUMUX_ZEROEQ_FLUX_VARIABLES_HH

#include <dumux/common/math.hh>
#include <dumux/common/valgrind.hh>

#include <dumux/freeflow/stokes/fluxvariables.hh>
#include <dumux/freeflow/zeroeq/indices.hh>
#include <dumux/freeflow/zeroeq/properties.hh>

namespace Dumux
{

/*!
 * \ingroup BoxZeroEqModel
 * \ingroup ImplicitFluxVariables
 * \brief This template class contains data which is required to
 *        calculate the component fluxes over a face of a finite
 *        volume for a ZeroEq model.
 *
 * This means the methods to calculate the eddy viscosity.
 */
template <class TypeTag>
class ZeroEqFluxVariables : public GET_PROP_TYPE(TypeTag, BaseStokesFluxVariables)
{
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, BaseStokesFluxVariables) ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;

    enum { dim = GridView::dimension };

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef Dune::FieldVector<Scalar, dim> DimVector;

public:
    //! \brief The old constructor
    DUNE_DEPRECATED_MSG("FluxVariables now have to be default constructed and updated.")
    ZeroEqFluxVariables(const Problem &problem,
                        const Element &element,
                        const FVElementGeometry &fvGeometry,
                        const int fIdx,
                        const ElementVolumeVariables &elemVolVars,
                        const bool onBoundary = false)
        : ParentType(problem, element, fvGeometry, fIdx, elemVolVars, onBoundary) {}

    /*!
     * \brief Default constructor
     * \note This can be removed when the deprecated constructor is removed.
     */
    ZeroEqFluxVariables() = default;

    /*!
     * \brief Compute / update the flux variables
     *
     * \param problem The problem
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param fIdx The local index of the SCV (sub-control-volume) face
     * \param elemVolVars The volume variables of the current element
     * \param onBoundary A boolean variable to specify whether the flux variables
     * are calculated for interior SCV faces or boundary faces, default=false
     */
    void update(const Problem &problem,
                const Element &element,
                const FVElementGeometry &fvGeometry,
                const int fIdx,
                const ElementVolumeVariables &elemVolVars,
                const bool onBoundary = false)
    {
        ParentType::update(problem, element, fvGeometry, fIdx, elemVolVars, onBoundary);

        flowNormal_ = GET_PARAM_FROM_GROUP(TypeTag, int, ZeroEq, FlowNormal);
        wallNormal_ = GET_PARAM_FROM_GROUP(TypeTag, int, ZeroEq, WallNormal);
        eddyViscosityModel_ = GET_PARAM_FROM_GROUP(TypeTag, int, ZeroEq, EddyViscosityModel);
        karmanConstant_= GET_PROP_VALUE(TypeTag, KarmanConstant);

        dynamicEddyViscosity_ = 0.0;
        mixingLength_ = 0.0;
        dynamicEddyViscosityInner_ = 0.0;
        dynamicEddyViscosityOuter_ = 0.0;
        fz_ = 0.0;
        posIdx_ = problem.model().getPosIdx(globalPos());
        wallIdx_ = problem.model().getWallIdx(globalPos(), posIdx_);
        distanceToWallRough_ = std::abs(problem.model().distanceToWallRough(globalPos(), wallIdx_, posIdx_));
        distanceToWallReal_ = std::abs(problem.model().distanceToWallReal(globalPos(), wallIdx_, posIdx_));
        for (int dimIdx = 0; dimIdx < dim; ++dimIdx)
            maxVelocity_[dimIdx] = problem.model().wall[wallIdx_].maxVelocity[posIdx_][dimIdx];
        for (int dimIdx = 0; dimIdx < dim; ++dimIdx)
            minVelocity_[dimIdx] = problem.model().wall[wallIdx_].minVelocity[posIdx_][dimIdx];
        velGrad_ = this->velocityGrad_[flowNormal_][wallNormal_];
        frictionVelocityWall_ = sqrt(problem.model().wall[wallIdx_].wallShearStress[posIdx_]
                                     / problem.model().wall[wallIdx_].wallDensity[posIdx_]);
        yPlusRough_ = distanceToWallRough_ * frictionVelocityWall_ / problem.model().wall[wallIdx_].wallKinematicViscosity[posIdx_];

        // calculation of an eddy viscosity only makes sense with Navier-Stokes equation
        if (GET_PROP_VALUE(TypeTag, EnableNavierStokes))
            asImp_().calculateEddyViscosity_(problem, element, elemVolVars);
    }

protected:
    //! Returns the implementation of the flux variables (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

    /*!
     * \brief This function calculates the dynamic viscosity.
     *
     * The eddy viscosity is added to the viscosity in stokeslocalresidual.hh at each scv
     * face.
     */
    void calculateEddyViscosity_(const Problem &problem,
                                 const Element &element,
                                 const ElementVolumeVariables &elemVolVars)
    {
        // no turbulence model
        if (eddyViscosityModel_ == EddyViscosityIndices::noEddyViscosityModel)
          return;

        // Prandtl mixing length
        // e.g. Wilcox, D. C., Turbulence Modeling for CFD, 2006
        else if (eddyViscosityModel_ == EddyViscosityIndices::prandtl)
        {
            mixingLength_ = distanceToWallRough_ * karmanConstant_;
            dynamicEddyViscosity_ = this->density() * mixingLength() * mixingLength() * std::abs(velGrad_);
        }

        // modified Van-Driest
        // e.g. Bird, Stewart, and Lightfoot, E. N. Transport phenomena, 2007
        else if (eddyViscosityModel_ == EddyViscosityIndices::modifiedVanDriest)
        {
            Scalar aPlus = 26.0;
            Scalar bPlus = 0.26;
            // eddy viscosity can only be calculated correctly for non-zero distance to walls
            mixingLength_ = 0.0;
            if (distanceToWallRough_ > 0.0 && yPlusRough_ > 0.0)
                mixingLength_= karmanConstant_ * distanceToWallRough_
                               * (1.0 - std::exp(-yPlusRough_ / aPlus ))
                               / std::sqrt(1.0 - std::exp(-bPlus * yPlusRough_));

            dynamicEddyViscosity_ = this->density() * mixingLength() * mixingLength() * std::abs(velGrad_);
        }

        // Baldwin and Lomax
        // Baldwin, B. S. & Lomax, H. "Thin Layer Approximation and Algebraic Model for Seperated Turbulent Flows"
        //   AIAA Journal, 1978, 78--257, 1-9
        else if (eddyViscosityModel_ == EddyViscosityIndices::baldwinLomax)
        {
            // LAW CONSTANTS
            const Scalar aPlus = 26.0;
            const Scalar cCP = 1.6;
            const Scalar cKleb = 0.30;
            const Scalar cWK = 0.25;
            const Scalar kUpper = 0.0168;

            // Calculate muInner
            mixingLength_ = 0.0;
            if (distanceToWallRough_ > 0.0 && yPlusRough_ > 0.0)
                mixingLength_ = karmanConstant_ * distanceToWallRough_
                                * (1.0 - std::exp(-yPlusRough_ / aPlus ));
            Scalar omega1 = this->velocityGrad_[0][1] - this->velocityGrad_[1][0];
            Scalar omega2 = 0.0;
            Scalar omega3 = 0.0;
            if (dim == 3)
            {
                omega2 = this->velocityGrad_[1][2] - this->velocityGrad_[2][1];
                omega3 = this->velocityGrad_[2][0] - this->velocityGrad_[0][2];
            }
            Scalar omega = sqrt(omega1 * omega1 + omega2 * omega2 + omega3 * omega3);
            dynamicEddyViscosityInner_ = this->density() * mixingLength()  * mixingLength() * omega;

            // Calculate muOuter
            fz_ = 0.0;
            if (distanceToWallRough_ > 0.0 && yPlusRough_ > 0.0)
                fz_ = distanceToWallRough_ * omega * (1.0 - std::exp(-yPlusRough_ / aPlus ));
            Scalar fMax = problem.model().wall[wallIdx_].fMax[posIdx_];
            Scalar yMax = std::abs(problem.model().wall[wallIdx_].yMax[posIdx_]);
            Scalar uDiff = maxVelocity_.two_norm() - minVelocity_.two_norm();

            Scalar f1 = yMax * fMax;
            Scalar f2 = cWK * yMax * uDiff * uDiff / fMax;
            Scalar fWake = fmin(f1, f2);
            Scalar fKleb = 1 / (1 + 5.5 * std::pow(cKleb * distanceToWallRough_ / yMax, 6.0));
            dynamicEddyViscosityOuter_ = this->density() * kUpper * cCP * fWake * fKleb;

            bool inner = problem.model().useViscosityInner(this->face().ipGlobal, posIdx_);

            // muOuter can only be calculated correctly for non-zero fmax
            if (fMax == 0)
                inner = true;

            if (inner)
                dynamicEddyViscosity_ = dynamicEddyViscosityInner_;
            else
                dynamicEddyViscosity_ = dynamicEddyViscosityOuter_;
        }

        else
        {
            DUNE_THROW(Dune::NotImplemented, "This eddy viscosity model is not implemented.");
        }

        Valgrind::CheckDefined(dynamicEddyViscosity_);
        Valgrind::CheckDefined(dynamicEddyViscosityInner_);
        Valgrind::CheckDefined(dynamicEddyViscosityOuter_);
    }

public:
    /*!
     * \brief Returns the global coordinates of the integration point.
     */
    DimVector globalPos() const
    { return this->face().ipGlobal; }

    /*!
     * \brief Returns the mixing length \f$\mathrm{[m]}\f$.
     */
    Scalar mixingLength() const
    { return mixingLength_; }

    /*!
     * \brief Returns the dynamic eddy viscosity
     *        \f$\mathrm{[Pa \cdot s]} = \mathrm{[N \cdot s/m^2]}\f$.
     */
    const Scalar dynamicEddyViscosity() const
    { return dynamicEddyViscosity_; }

    /*!
     * \brief Returns the kinematic eddy viscosity
     *        \f$\mathrm{[m^2/s]}\f$ (if implemented).
     */
    const Scalar kinematicEddyViscosity() const
    { return dynamicEddyViscosity() / this->density(); }

    /*!
     * \brief Returns the inner dynamic eddy Viscosity (only Baldwin-Lomax model)
     *        \f$\mu_\textrm{inner}\f$ in \f$\mathrm{[N\cdot s/m^2]}\f$.
     */
    Scalar dynamicEddyViscosityInner() const
    { return dynamicEddyViscosityInner_; }

    /*!
     * \brief Returns the outer dynamic eddy Viscosity (only Baldwin-Lomax model).
     *        \f$\mu_\textrm{outer}\f$ in \f$\mathrm{[N\cdot s/m^2]}\f$.
     */
    Scalar dynamicEddyViscosityOuter() const
    { return dynamicEddyViscosityOuter_; }

    /*!
     * \brief Returns the value of the f-function.
     *
     * \f$\mathrm{f = y \omega \left( 1 - exp \left[ -y^+ / A^+ \right] \right)}\f$.<br>
     * y = distanceToWall
     */
    Scalar fz() const
    { return fz_; }

    /*!
     * \brief Returns the friction velocity at the wall \f$\mathrm{[kg/(m*s^2)]}\f$.
     */
    Scalar frictionVelocityWall() const
    { return frictionVelocityWall_; }

    /*!
     * \brief Returns the distance to the corresponding Wall, including surface roughness \f$\mathrm{[m]}\f$.
     */
    Scalar distanceToWallRough() const
    { return distanceToWallRough_; }

    /*!
     * \brief Returns the real distance to the corresponding Wall \f$\mathrm{[m]}\f$.
     */
    Scalar distanceToWallReal() const
    { return distanceToWallReal_; }

    /*!
     * \brief Returns dimensionless wall distance, including surface roughness \f$\mathrm{[-]}\f$.
     */
    Scalar yPlusRough() const
    { return yPlusRough_; }

    /*!
     * \brief Returns a dimensionless velocity \f$\mathrm{[-]}\f$.
     */
    Scalar uPlus() const
    { return this->velocity()[flowNormal_] / frictionVelocityWall(); }

    /*!
     * \brief Returns the Karman constant \f$\mathrm{[-]}\f$.
     */
    Scalar karmanConstant() const
    { return karmanConstant_; }

private:
    int flowNormal_;
    int wallNormal_;
    int eddyViscosityModel_;
    Scalar karmanConstant_;

    int wallIdx_;
    int posIdx_;
    Scalar yPlusRough_;
    Scalar distanceToWallReal_;
    Scalar distanceToWallRough_;
    Scalar sandGrainRoughness_;
    Scalar sandGrainRoughnessDimensionless_;

    Scalar velGrad_;
    Scalar frictionVelocityWall_;
    DimVector maxVelocity_;
    DimVector minVelocity_;

    Scalar mixingLength_;
    Scalar dynamicEddyViscosity_;
    Scalar dynamicEddyViscosityInner_;
    Scalar dynamicEddyViscosityOuter_;
    Scalar fz_;

    Scalar eps_;
};

} // end namespace

#endif // DUMUX_ZEROEQ_FLUX_VARIABLES_HH
