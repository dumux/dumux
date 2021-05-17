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
 * \ingroup NavierStokesModel
 * \copydoc Dumux::NavierStokesProblem
 */
#ifndef DUMUX_NAVIERSTOKES_PROBLEM_HH
#define DUMUX_NAVIERSTOKES_PROBLEM_HH

#include <dune/common/exceptions.hh>
#include <dune/common/typetraits.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/staggeredfvproblem.hh>
#include <dumux/discretization/method.hh>

namespace Dumux {

//! The implementation is specialized for the different discretizations
template<class TypeTag, DiscretizationMethod discMethod> struct NavierStokesParentProblemImpl;

template<class TypeTag>
struct NavierStokesParentProblemImpl<TypeTag, DiscretizationMethod::staggered>
{
    using type = StaggeredFVProblem<TypeTag>;
};

//! The actual NavierStokesParentProblem
template<class TypeTag>
using NavierStokesParentProblem =
      typename NavierStokesParentProblemImpl<TypeTag,
      GetPropType<TypeTag, Properties::GridGeometry>::discMethod>::type;

/*!
 * \ingroup NavierStokesModel
 * \brief Navier-Stokes problem base class.
 *
 * This implements gravity (if desired) and a function returning the temperature.
 * Includes a specialized method used only by the staggered grid discretization.
 *
 */
template<class TypeTag>
class NavierStokesProblem : public NavierStokesParentProblem<TypeTag>
{
    using ParentType = NavierStokesParentProblem<TypeTag>;
    using Implementation = GetPropType<TypeTag, Properties::Problem>;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using GridFaceVariables = typename GridVariables::GridFaceVariables;
    using ElementFaceVariables = typename GridFaceVariables::LocalView;
    using FaceVariables = typename GridFaceVariables::FaceVariables;
    using GridVolumeVariables = typename GridVariables::GridVolumeVariables;
    using ElementVolumeVariables = typename GridVolumeVariables::LocalView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
      };

    using GlobalPosition = typename SubControlVolumeFace::GlobalPosition;
    using VelocityVector = Dune::FieldVector<Scalar, dimWorld>;
    using GravityVector = Dune::FieldVector<Scalar, dimWorld>;

public:
    /*!
     * \brief The constructor
     * \param gridGeometry The finite volume grid geometry
     * \param paramGroup The parameter group in which to look for runtime parameters first (default is "")
     */
    NavierStokesProblem(std::shared_ptr<const GridGeometry> gridGeometry, const std::string& paramGroup = "")
    : ParentType(gridGeometry, paramGroup)
    , gravity_(0.0)
    {
        if (getParamFromGroup<bool>(paramGroup, "Problem.EnableGravity"))
            gravity_[dim-1]  = -9.81;

        enableInertiaTerms_ = getParamFromGroup<bool>(paramGroup, "Problem.EnableInertiaTerms");
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

    //needed for residualcalc
    auto instationaryAnalyticalSolutionAtPos(const GlobalPosition& globalPos, Scalar t) const
    {
        return asImp_().analyticalSolutionAtPos(globalPos);
    }

    /*!
     * \brief Return the analytical solution of the problem at a given position
     *
     * \param globalPos The global position
     */
    PrimaryVariables analyticalSolution(const GlobalPosition& globalPos) const
    {
        if (!getParam<bool>("Problem.AveragedAnalyticalSolution", true))
            return asImp_().analyticalSolutionAtPos(globalPos);

        return meanAnalyticalSolution_(globalPos, [&](Scalar x) { return asImp_().xFactorAnalyticalSolutionAntiderivativeAtPos(x); },
                                                 [&](Scalar y) { return asImp_().yFactorAnalyticalSolutionAntiderivativeAtPos(y); },
                                                 [&](Scalar x) { return asImp_().xFactorAnalyticalSolutionAtPos(x); },
                                                 [&](Scalar y) { return asImp_().yFactorAnalyticalSolutionAtPos(y); });
    }

    /*!
     * \brief Return the analytical solution of the problem at a given position
     *
     * \param globalPos The global position
     */
    PrimaryVariables instationaryAnalyticalSolution(const GlobalPosition& globalPos, Scalar t) const
    {
        if (!getParam<bool>("Problem.AveragedAnalyticalSolution", true))
            return asImp_().instationaryAnalyticalSolutionAtPos(globalPos,t);

        return meanAnalyticalSolution_(globalPos,
                                      [&](Scalar x) { return asImp_().xFactorAnalyticalSolutionAntiderivativeAtPos(x,t); },
                                      [&](Scalar y) { return asImp_().yFactorAnalyticalSolutionAntiderivativeAtPos(y,t); },
                                      [&](Scalar x) { return asImp_().xFactorAnalyticalSolutionAtPos(x,t); },
                                      [&](Scalar y) { return asImp_().yFactorAnalyticalSolutionAtPos(y,t); });
    }

    /*!
     * \brief Returns the acceleration due to gravity.
     *
     * If the <tt>Problem.EnableGravity</tt> parameter is true, this means
     * \f$\boldsymbol{g} = ( 0,\dots,\ -9.81)^T \f$, else \f$\boldsymbol{g} = ( 0,\dots, 0)^T \f$
     */
    const GravityVector& gravity() const
    { return gravity_; }

    /*!
     * \brief Returns whether interia terms should be considered.
     */
    bool enableInertiaTerms() const
    { return enableInertiaTerms_; }

    //! Applys the initial face solution (velocities on the faces). Specialization for staggered grid discretization.
    template <class SolutionVector, class G = GridGeometry>
    typename std::enable_if<G::discMethod == DiscretizationMethod::staggered, void>::type
    applyInitialFaceSolution(SolutionVector& sol,
                             const SubControlVolumeFace& scvf,
                             const PrimaryVariables& initSol) const
    {
        sol[GridGeometry::faceIdx()][scvf.dofIndex()][0] = initSol[Indices::velocity(scvf.directionIndex())];
    }

    /*!
     * \brief An additional drag term can be included as source term for the momentum balance
     *        to mimic 3D flow behavior in 2D:
     *  \f[
     *        f_{drag} = -(8 \mu / h^2)v
     *  \f]
     *  Here, \f$h\f$ corresponds to the extruded height that is
     *  bounded by the imaginary walls. See Flekkoy et al. (1995) \cite flekkoy1995a<BR>
     *  A value of 8.0 is used as a default factor, corresponding
     *  to the velocity profile at  the center plane
     *  of the virtual height (maximum velocity). Setting this value to 12.0 corresponds
     *  to an depth-averaged velocity (Venturoli and Boek, 2006) \cite venturoli2006a.
     */
    Scalar pseudo3DWallFriction(const Scalar velocity,
                                const Scalar viscosity,
                                const Scalar height,
                                const Scalar factor = 8.0) const
    {
        static_assert(dim == 2, "Pseudo 3D wall friction may only be used in 2D");
        return -factor * velocity * viscosity / (height*height);
    }

    //! Convenience function for staggered grid implementation.
    template <class ElementVolumeVariables, class ElementFaceVariables, class G = GridGeometry>
    typename std::enable_if<G::discMethod == DiscretizationMethod::staggered, Scalar>::type
    pseudo3DWallFriction(const SubControlVolumeFace& scvf,
                         const ElementVolumeVariables& elemVolVars,
                         const ElementFaceVariables& elemFaceVars,
                         const Scalar height,
                         const Scalar factor = 8.0) const
    {
        const Scalar velocity = elemFaceVars[scvf].velocitySelf();
        const Scalar viscosity = elemVolVars[scvf.insideScvIdx()].effectiveViscosity();
        return pseudo3DWallFriction(velocity, viscosity, height, factor);
    }

    /*!
     * \brief Returns the intrinsic permeability of required as input parameter for the Beavers-Joseph-Saffman boundary condition
     *
     * This member function must be overloaded in the problem implementation, if the BJS boundary condition is used.
     */
    Scalar permeability(const Element& element, const SubControlVolumeFace& scvf) const
    {
        DUNE_THROW(Dune::NotImplemented, "When using the Beavers-Joseph-Saffman boundary condition, the permeability must be returned in the acutal problem");
    }

    /*!
     * \brief Returns the alpha value required as input parameter for the Beavers-Joseph-Saffman boundary condition
     *
     * This member function must be overloaded in the problem implementation, if the BJS boundary condition is used.
     */
    Scalar alphaBJ(const SubControlVolumeFace& scvf) const
    {
        DUNE_THROW(Dune::NotImplemented, "When using the Beavers-Joseph-Saffman boundary condition, the alpha value must be returned in the acutal problem");
    }

    /*!
     * \brief Returns the beta value which is the alpha value divided by the square root of the (scalar-valued) interface permeability.
     */
    Scalar betaBJ(const Element& element, const SubControlVolumeFace& scvf, const GlobalPosition& tangentialVector) const
    {
        const Scalar interfacePermeability = interfacePermeability_(element, scvf, tangentialVector);
        using std::sqrt;
        return asImp_().alphaBJ(scvf) / sqrt(interfacePermeability);
    }

    /*!
     * \brief Returns the velocity in the porous medium (which is 0 by default according to Saffmann).
     */
    VelocityVector porousMediumVelocity(const Element& element, const SubControlVolumeFace& scvf) const
    {
        return VelocityVector(0.0);
    }

    /*!
     * \brief Returns the slip velocity at a porous boundary based on the Beavers-Joseph(-Saffman) condition.
     */
    const Scalar beaversJosephVelocity(const Element& element,
                                       const SubControlVolume& scv,
                                       const SubControlVolumeFace& ownScvf,
                                       const SubControlVolumeFace& faceOnPorousBoundary,
                                       const Scalar velocitySelf,
                                       const Scalar tangentialVelocityGradient) const
    {
        // create a unit normal vector oriented in positive coordinate direction
        GlobalPosition orientation = ownScvf.unitOuterNormal();
        orientation[ownScvf.directionIndex()] = 1.0;

        // du/dy + dv/dx = alpha/sqrt(K) * (u_boundary-uPM)
        // beta = alpha/sqrt(K)
        const Scalar betaBJ = asImp_().betaBJ(element, faceOnPorousBoundary, orientation);
        const Scalar distanceNormalToBoundary = (faceOnPorousBoundary.center() - scv.center()).two_norm();

        return (tangentialVelocityGradient*distanceNormalToBoundary
              + asImp_().porousMediumVelocity(element, faceOnPorousBoundary) * orientation * betaBJ * distanceNormalToBoundary
              + velocitySelf) / (betaBJ*distanceNormalToBoundary + 1.0);
    }

private:
    /*!
     * \brief Return the analytical solution of the problem at a given position
     *
     * \param globalPos The global position
     */
    template <class LambdaA, class LambdaB, class LambdaC, class LambdaD>
    PrimaryVariables meanAnalyticalSolution_(const GlobalPosition& globalPos,
                                        const LambdaA& xFactorAnalyticalSolutionAntiderivativeAtPosLambdaFunction,
                                        const LambdaB& yFactorAnalyticalSolutionAntiderivativeAtPosLambdaFunction,
                                        const LambdaC& xFactorAnalyticalSolutionAtPosLambdaFunction,
                                        const LambdaD& yFactorAnalyticalSolutionAtPosLambdaFunction) const
    {
        Scalar xIntegralLeft = 0.;
        Scalar xIntegralRight = 0.;
        Scalar yIntegralDown = 0.;
        Scalar yIntegralUp = 0.;

        bool integrateX = false;
        bool integrateY = false;

        setIntegralRanges_(globalPos, integrateX, integrateY, xIntegralLeft, xIntegralRight, yIntegralDown, yIntegralUp);

        std::array<PrimaryVariables,2> xMean;
        std::array<PrimaryVariables,2> yMean;

        for (unsigned int j=0; j < 2; ++j) //up to two summands in the analytical solution
        {
            for (unsigned int i = 0; i < xMean[j].size(); ++i)
            {
                xMean[j][i] = (xFactorAnalyticalSolutionAntiderivativeAtPosLambdaFunction(xIntegralRight)[j][i] - xFactorAnalyticalSolutionAntiderivativeAtPosLambdaFunction(xIntegralLeft)[j][i])/(xIntegralRight - xIntegralLeft);
                yMean[j][i] = (yFactorAnalyticalSolutionAntiderivativeAtPosLambdaFunction(yIntegralUp)[j][i]    - yFactorAnalyticalSolutionAntiderivativeAtPosLambdaFunction(yIntegralDown)[j][i])/(yIntegralUp    - yIntegralDown);
            }
        }

        const std::array<PrimaryVariables,2> xPointValue = xFactorAnalyticalSolutionAtPosLambdaFunction(globalPos[0]);
        const std::array<PrimaryVariables,2> yPointValue = yFactorAnalyticalSolutionAtPosLambdaFunction(globalPos[1]);

        if (integrateX && !integrateY)
        {
            PrimaryVariables values;

            for (unsigned int j=0; j < 2; ++j)
            {
                for (unsigned int i = 0; i < values.size(); ++i)
                {
                    values[i] += xMean[j][i] * yPointValue[j][i];
                }
            }

            return values;
        }

        if (!integrateX && integrateY)
        {
            PrimaryVariables values;

            for (unsigned int j=0; j < 2; ++j)
            {
                for (unsigned int i = 0; i < values.size(); ++i)
                {
                    values[i] += xPointValue[j][i] * yMean[j][i];
                }
            }

            return values;
        }

        if (integrateX && integrateY)
        {
            PrimaryVariables values;

            for (unsigned int j=0; j < 2; ++j)
            {
                for (unsigned int i = 0; i < values.size(); ++i)
                {
                    values[i] += xMean[j][i] * yMean[j][i];
                }
            }

            return values;
        }
    }

    template<class SomeType>
    bool containerCmp_(const SomeType& a, const SomeType& b) const
    {
        double eps = 1e-10;

        if (a.size() != b.size())
        {
            std::cout << "containerCmp_ with different sizes" << std::endl;
            return false;
        }
        else
        {
            bool retVal = true;
            for (unsigned int i = 0; i < a.size(); ++i)
            {
                if (!((a[i] > b[i] - eps) && (a[i] < b[i] + eps)))
                {
                    retVal = false;
                    break;
                }
            }

            return retVal;
        }
    }

    template<class SomeType>
    bool scalarCmp_(const SomeType& a, const SomeType& b) const
    {
        double eps = 1e-10;

        return a > b - eps
        && a < b + eps;
    }

    void setVolumeOrLineForAveragingCorners_(const GlobalPosition& globalPos, std::array<GlobalPosition,4>& volumeForAveragingCorners, std::array<GlobalPosition,2>& lineForAveragingCorners, bool& setVolumeForAveragingCorners, bool& setLineForAveragingCorners) const
    {
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            auto fvGeometry = localView(this->gridGeometry());
            fvGeometry.bindElement(element);
            for (auto&& scv : scvs(fvGeometry))
            {
                if (containerCmp_(globalPos,scv.center()))
                {
                    if(scv.geometry().corners()!=4)
                    {
                        DUNE_THROW(Dune::InvalidStateException, "");
                    }

                    for (unsigned int i = 0; i < scv.geometry().corners(); ++i)
                    {
                        volumeForAveragingCorners[i] = scv.geometry().corner(i);
                    }

                    setVolumeForAveragingCorners = true;
                }

                for (auto&& scvf : scvfs(fvGeometry))
                {
                    if (!scvf.boundary())
                    {
                        if (containerCmp_(globalPos,scvf.center()))
                        {
                            lineForAveragingCorners[0] = scvf.corner(0);
                            lineForAveragingCorners[1] = scvf.corner(1);
                            setLineForAveragingCorners = true;
                        }
                    }
                    else
                    {
                        if (containerCmp_(globalPos,scvf.center()))
                        {
                            lineForAveragingCorners[0] = scvf.corner(0);
                            lineForAveragingCorners[1] = scvf.corner(1);
                            setLineForAveragingCorners = true;
                        }

                        auto treatVirtualParallelFaceDofPos = [&](unsigned int thisCornerIdx)
                        {
                            const GlobalPosition& thisCorner = scvf.corner(thisCornerIdx);

                            if (containerCmp_(globalPos, thisCorner))
                            {
                                lineForAveragingCorners[0] = scvf.center();

                                unsigned int dirIdx = scvf.directionIndex();
                                unsigned int nonDirIdx = (dirIdx == 0) ? 1 : 0;

                                unsigned int otherCornerIdx = (thisCornerIdx == 0) ? 1 : 0;
                                const GlobalPosition& otherCorner = scvf.corner(otherCornerIdx);

                                const SubControlVolumeFace& normalFacePair0 = this->gridGeometry().scvf(scvf.insideScvIdx(), scvf.pairData(0).localLateralFaceIdx);
                                const SubControlVolumeFace& normalFacePair1 = this->gridGeometry().scvf(scvf.insideScvIdx(), scvf.pairData(1).localLateralFaceIdx);

                                unsigned int localSubFaceIdx;

                                if (scalarCmp_(thisCorner[nonDirIdx], normalFacePair0.center()[nonDirIdx]))
                                {
                                    localSubFaceIdx = thisCornerIdx;
                                }
                                else if (scalarCmp_(thisCorner[nonDirIdx], normalFacePair1.center()[nonDirIdx]))
                                {
                                    localSubFaceIdx = otherCornerIdx;
                                }
                                else
                                {
                                    DUNE_THROW(Dune::InvalidStateException, "");
                                }

                                if (thisCorner[nonDirIdx] > otherCorner[nonDirIdx] )
                                {
                                    lineForAveragingCorners[1] = scvf.center();
                                    lineForAveragingCorners[1][nonDirIdx] += scvf.parallelDofsDistance(localSubFaceIdx, 0);
                                }
                                else
                                {
                                    lineForAveragingCorners[1] = scvf.center();
                                    lineForAveragingCorners[1][nonDirIdx] -= scvf.parallelDofsDistance(localSubFaceIdx, 0);
                                }

                                setLineForAveragingCorners = true;
                            }
                        };

                        treatVirtualParallelFaceDofPos(0);
                        treatVirtualParallelFaceDofPos(1);
                    }
                }
            }
        }

        if (setLineForAveragingCorners == setVolumeForAveragingCorners)
        {
            DUNE_THROW(Dune::InvalidStateException, "either they both stayed false or they were both set to true, none of which should happen");
        }
    }

    void setIntegralRanges_ (const GlobalPosition& globalPos, bool& integrateX, bool& integrateY, Scalar& xIntegralLeft, Scalar& xIntegralRight, Scalar& yIntegralDown, Scalar& yIntegralUp) const
    {
        std::array<GlobalPosition,4> volumeForAveragingCorners = {};
        std::array<GlobalPosition,2> lineForAveragingCorners = {};

        bool setVolumeForAveragingCorners = false;
        bool setLineForAveragingCorners = false;

        setVolumeOrLineForAveragingCorners_(globalPos, volumeForAveragingCorners, lineForAveragingCorners, setVolumeForAveragingCorners, setLineForAveragingCorners);

        if (setLineForAveragingCorners)
        {
            if (scalarCmp_(lineForAveragingCorners[0][1], lineForAveragingCorners[1][1]))
            {
                integrateX = true;
                if (lineForAveragingCorners[0][0] < lineForAveragingCorners[1][0])
                {
                    xIntegralLeft = lineForAveragingCorners[0][0];
                    xIntegralRight = lineForAveragingCorners[1][0];
                }
                else if (lineForAveragingCorners[0][0] > lineForAveragingCorners[1][0])
                {
                    xIntegralLeft = lineForAveragingCorners[1][0];
                    xIntegralRight = lineForAveragingCorners[0][0];
                }
                else
                {
                    DUNE_THROW(Dune::InvalidStateException, "");
                }
            }
            else if (scalarCmp_(lineForAveragingCorners[0][0], lineForAveragingCorners[1][0]))
            {
                integrateY = true;
                if (lineForAveragingCorners[0][1] < lineForAveragingCorners[1][1])
                {
                    yIntegralDown = lineForAveragingCorners[0][1];
                    yIntegralUp = lineForAveragingCorners[1][1];
                }
                else if (lineForAveragingCorners[0][1] > lineForAveragingCorners[1][1])
                {
                    yIntegralDown = lineForAveragingCorners[1][1];
                    yIntegralUp = lineForAveragingCorners[0][1];
                }
                else
                {
                    DUNE_THROW(Dune::InvalidStateException, "");
                }
            }
            else
            {
                DUNE_THROW(Dune::InvalidStateException, "");
            }
        }
        else //setVolumeForAveragingCorners
        {
            integrateX = true;
            integrateY = true;

            xIntegralLeft = std::numeric_limits<Scalar>::max();
            xIntegralRight = std::numeric_limits<Scalar>::min();
            yIntegralDown = std::numeric_limits<Scalar>::max();
            yIntegralUp = std::numeric_limits<Scalar>::min();

            for (unsigned int i = 0; i < volumeForAveragingCorners.size(); ++i)
            {
                xIntegralLeft  = std::min(xIntegralLeft,  volumeForAveragingCorners[i][0]);
                xIntegralRight = std::max(xIntegralRight, volumeForAveragingCorners[i][0]);
                yIntegralDown  = std::min(yIntegralDown,  volumeForAveragingCorners[i][1]);
                yIntegralUp    = std::max(yIntegralUp,    volumeForAveragingCorners[i][1]);
            }
        }
    }

    //! Returns a scalar permeability value at the coupling interface
    Scalar interfacePermeability_(const Element& element, const SubControlVolumeFace& scvf, const GlobalPosition& tangentialVector) const
    {
        const auto& K = asImp_().permeability(element, scvf);

        // use t*K*t for permeability tensors
        if constexpr (Dune::IsNumber<std::decay_t<decltype(K)>>::value)
            return K;
        else
            return vtmv(tangentialVector, K, tangentialVector);
    }

    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

    GravityVector gravity_;
    bool enableInertiaTerms_;
};

} // end namespace Dumux

#endif
