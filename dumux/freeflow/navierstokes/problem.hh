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
     * \brief Returns the beta value which is the alpha value divided by the square root of the (scalar-valued) interface permeability.
     */
    [[deprecated("Use betaBJ with tangential vector instead. Will be removed after 3.3")]]
    Scalar betaBJ(const Element& element, const SubControlVolumeFace& scvf) const
    {
        using std::sqrt;
        return asImp_().alphaBJ(scvf) / sqrt(asImp_().permeability(element, scvf));
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
