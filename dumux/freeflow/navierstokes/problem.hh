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
#include <dumux/common/numeqvector.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/fvproblemwithspatialparams.hh>
#include <dumux/common/staggeredfvproblem.hh>
#include <dumux/discretization/method.hh>
#include <dumux/freeflow/navierstokes/momentum/boundarytypes.hh>

namespace Dumux {

//! The implementation is specialized for the different discretizations
template<class TypeTag, class DiscretizationMethod> struct NavierStokesParentProblemImpl;

// compatibility with old-style Navier-Stokes models
template<class TypeTag>
struct NavierStokesParentProblemImpl<TypeTag, DiscretizationMethods::Staggered>
{
    using type = StaggeredFVProblem<TypeTag>;
};

//! The actual NavierStokesParentProblem
template<class TypeTag>
using NavierStokesParentProblem = typename NavierStokesParentProblemImpl<
    TypeTag, typename GetPropType<TypeTag, Properties::GridGeometry>::DiscretizationMethod
>::type;

template<class TypeTag, class DiscretizationMethod>
class NavierStokesProblemImpl;

template<class TypeTag>
class NavierStokesProblemImpl<TypeTag, DiscretizationMethods::FCStaggered>
: public FVProblemWithSpatialParams<TypeTag>
{
    using ParentType = FVProblemWithSpatialParams<TypeTag>;
    using Implementation = GetPropType<TypeTag, Properties::Problem>;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using GridVolumeVariables = typename GridVariables::GridVolumeVariables;
    using ElementVolumeVariables = typename GridVolumeVariables::LocalView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
      };

    using GlobalPosition = typename SubControlVolumeFace::GlobalPosition;
    using VelocityVector = Dune::FieldVector<Scalar, dimWorld>;
    using GravityVector = Dune::FieldVector<Scalar, dimWorld>;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;

    static constexpr bool isCoupled_ = !std::is_empty_v<CouplingManager>;


public:
    //! These types are used in place of the typical NumEqVector type.
    //! In the numeqvector assembly type, only one equation per DOF (face) is considered
    //! while the type here provides one entry for each world dimension.
    using InitialValues = Dune::FieldVector<Scalar, dimWorld>;
    using Sources = Dune::FieldVector<Scalar, dimWorld>;
    using DirichletValues = Dune::FieldVector<Scalar, dimWorld>;
    using BoundaryFluxes = Dune::FieldVector<Scalar, dimWorld>;

    using MomentumFluxType = Dune::FieldVector<Scalar, dimWorld>;

    //! Export the boundary types.
    using BoundaryTypes = NavierStokesMomentumBoundaryTypes<ModelTraits::dim()>;

    //! This problem is used for the momentum balance model.
    static constexpr bool isMomentumProblem() { return true; }

    /*!
     * \brief The constructor
     * \param gridGeometry The finite volume grid geometry
     * \param couplingManager The coupling manager (couples mass and momentum equations)
     * \param paramGroup The parameter group in which to look for runtime parameters first (default is "")
     */
    NavierStokesProblemImpl(std::shared_ptr<const GridGeometry> gridGeometry,
                            std::shared_ptr<CouplingManager> couplingManager,
                            const std::string& paramGroup = "")
    : ParentType(gridGeometry, paramGroup)
    , gravity_(0.0)
    , couplingManager_(couplingManager)
    {
        if (getParamFromGroup<bool>(paramGroup, "Problem.EnableGravity"))
            gravity_[dim-1]  = -9.81;

        enableInertiaTerms_ = getParamFromGroup<bool>(paramGroup, "Problem.EnableInertiaTerms");
    }

    /*!
     * \brief The constructor for usage without a coupling manager
     * \param gridGeometry The finite volume grid geometry
     * \param paramGroup The parameter group in which to look for runtime parameters first (default is "")
     */
    NavierStokesProblemImpl(std::shared_ptr<const GridGeometry> gridGeometry,
                            const std::string& paramGroup = "")
    : NavierStokesProblemImpl(gridGeometry, {}, paramGroup)
    {}

    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     *
     * This is the method for the case where the source term is
     * potentially solution dependent and requires some quantities that
     * are specific to the fully-implicit method.
     *
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars All volume variables for the element
     * \param scv The sub control volume
     *
     * For this method, the return parameter stores the conserved quantity rate
     * generated or annihilate per volume unit. Positive values mean
     * that the conserved quantity is created, negative ones mean that it vanishes.
     * E.g. for the mass balance that would be a mass rate in \f$ [ kg / (m^3 \cdot s)] \f$.
     */
    Sources source(const Element& element,
                   const FVElementGeometry& fvGeometry,
                   const ElementVolumeVariables& elemVolVars,
                   const SubControlVolume& scv) const
    {
        // forward to solution independent, fully-implicit specific interface
        return asImp_().sourceAtPos(scv.dofPosition());
    }

    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     *
     * \param globalPos The position of the center of the finite volume
     *            for which the source term ought to be
     *            specified in global coordinates
     *
     * For this method, the values parameter stores the conserved quantity rate
     * generated or annihilate per volume unit. Positive values mean
     * that the conserved quantity is created, negative ones mean that it vanishes.
     * E.g. for the mass balance that would be a mass rate in \f$ [ kg / (m^3 \cdot s)] \f$.
     */
    Sources sourceAtPos(const GlobalPosition& globalPos) const
    {
        //! As a default, i.e. if the user's problem does not overload any source method
        //! return 0.0 (no source terms)
        return Sources(0.0);
    }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param element The finite element
     * \param scvf The sub control volume face
     */
    auto boundaryTypes(const Element& element,
                       const SubControlVolumeFace& scvf) const
    {
        // Forward it to the method which only takes the global coordinate.
        // We evaluate the boundary type at the center of the sub control volume face
        // in order to avoid ambiguities at domain corners.
        return asImp_().boundaryTypesAtPos(scvf.center());
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume face.
     *
     * \param element The finite element
     * \param scvf the sub control volume face
     * \note used for cell-centered discretization schemes
     */
    DirichletValues dirichlet(const Element& element, const SubControlVolumeFace& scvf) const
    {
        return asImp_().dirichletAtPos(scvf.ipGlobal());
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann control volume.
     *
     * \param element The element for which the Neumann boundary condition is set
     * \param fvGeometry The fvGeometry
     * \param elemVolVars The element volume variables
     * \param elemFluxVarsCache The element flux variables cache
     * \param scvf The boundary sub control volume face
     */
    template<class ElementFluxVariablesCache>
    BoundaryFluxes neumann(const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const ElementVolumeVariables& elemVolVars,
                           const ElementFluxVariablesCache& elemFluxVarsCache,
                           const SubControlVolumeFace& scvf) const
    { return asImp_().neumannAtPos(scvf.ipGlobal()); }

    /*!
     * \brief Returns the neumann flux at a given position.
     */
    BoundaryFluxes neumannAtPos(const GlobalPosition& globalPos) const
    { return BoundaryFluxes(0.0); } //! A default, i.e. if the user's does not overload any neumann method

    /*!
     * \brief Returns the acceleration due to gravity.
     *
     * If the <tt>Problem.EnableGravity</tt> parameter is true, this means
     * \f$\boldsymbol{g} = ( 0,\dots,\ -9.81)^T \f$, else \f$\boldsymbol{g} = ( 0,\dots, 0)^T \f$
     */
    const GravityVector& gravity() const
    { return gravity_; }

    /*!
     * \brief Returns whether inertia terms should be considered.
     */
    bool enableInertiaTerms() const
    { return enableInertiaTerms_; }

    /*!
     * \brief Returns the pressure at a given sub control volume face.
     * \note  Overload this if a fixed pressure shall be prescribed (e.g., given by an analytical solution).
     */
    Scalar pressure(const Element& element,
                    const FVElementGeometry& fvGeometry,
                    const SubControlVolumeFace& scvf) const
    {
        if constexpr (isCoupled_)
            return couplingManager_->pressure(element, fvGeometry, scvf);
        else
            return asImp_().pressureAtPos(scvf.ipGlobal());
    }

    /*!
     * \brief Returns the pressure at a given position.
     */
    Scalar pressureAtPos(const GlobalPosition&) const
    {
        DUNE_THROW(Dune::NotImplemented, "pressureAtPos not implemented");
    }

    /*!
     * \brief Returns a reference pressure at a given sub control volume face.
     *        This pressure is subtracted from the actual pressure for the momentum balance
     *        which potentially helps to improve numerical accuracy by avoiding issues related do floating point arithmetic.
     * \note  Overload this for reference pressures other than zero.
     */
    Scalar referencePressure(const Element& element,
                             const FVElementGeometry& fvGeometry,
                             const SubControlVolumeFace& scvf) const
    { return 0.0; }

    /*!
     * \brief Returns the density at a given sub control volume face.
     * \note  Overload this if a fixed density shall be prescribed.
     */
    Scalar density(const Element& element,
                   const FVElementGeometry& fvGeometry,
                   const SubControlVolumeFace& scvf) const
    {
        if constexpr (isCoupled_)
            return couplingManager_->density(element, fvGeometry, scvf);
        else
            return asImp_().densityAtPos(scvf.ipGlobal());
    }

    /*!
     * \brief Returns the density at a given sub control volume.
     * \note  Overload this if a fixed density shall be prescribed.
     */
    Scalar density(const Element& element,
                   const SubControlVolume& scv,
                   const bool isPreviousTimeStep = false) const
    {
        if constexpr (isCoupled_)
            return couplingManager_->density(element, scv, isPreviousTimeStep);
        else
            return asImp_().densityAtPos(scv.dofPosition());
    }

    auto insideAndOutsideDensity(const Element& element,
                                 const FVElementGeometry& fvGeometry,
                                 const SubControlVolumeFace& scvf,
                                 const bool isPreviousTimeStep = false) const
    {
        if constexpr (isCoupled_)
            return couplingManager_->insideAndOutsideDensity(element, fvGeometry, scvf, isPreviousTimeStep);
        else
        {
            const auto rho = asImp_().densityAtPos(scvf.ipGlobal());
            return std::make_pair(rho, rho);
        }
    }

    /*!
     * \brief Returns the density at a given position.
     */
    Scalar densityAtPos(const GlobalPosition&) const
    {
        DUNE_THROW(Dune::NotImplemented, "densityAtPos not implemented");
    }

    /*!
     * \brief Returns the effective dynamic viscosity at a given sub control volume face.
     * \note  Overload this if a fixed viscosity shall be prescribed.
     */
    Scalar effectiveViscosity(const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const SubControlVolumeFace& scvf) const
    {
        if constexpr (isCoupled_)
            return couplingManager_->effectiveViscosity(element, fvGeometry, scvf);
        else
            return asImp_().effectiveViscosityAtPos(scvf.ipGlobal());
    }

    /*!
     * \brief Returns the effective dynamic viscosity at a given position.
     */
    Scalar effectiveViscosityAtPos(const GlobalPosition&) const
    {
        DUNE_THROW(Dune::NotImplemented, "effectiveViscosityAtPos not implemented");
    }

    /*!
     * \brief Applies the initial solution for all degrees of freedom of the grid.
     * \param sol the initial solution vector
     */
    template<class SolutionVector>
    void applyInitialSolution(SolutionVector& sol) const
    {
        sol.resize(this->gridGeometry().numDofs());
        std::vector<bool> dofHandled(this->gridGeometry().numDofs(), false);
        auto fvGeometry = localView(this->gridGeometry());
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            fvGeometry.bindElement(element);
            for (const auto& scv : scvs(fvGeometry))
            {
                const auto dofIdx = scv.dofIndex();
                if (!dofHandled[dofIdx])
                {
                    dofHandled[dofIdx] = true;
                    sol[dofIdx] = asImp_().initial(scv)[scv.dofAxis()];
                }
            }
        }
    }

    /*!
     * \brief Evaluate the initial value at an sub control volume
     */
    InitialValues initial(const SubControlVolume& scv) const
    {
        return asImp_().initialAtPos(scv.dofPosition());
    }

    //! Convenience function for staggered grid implementation.
    Scalar pseudo3DWallFriction(const Element& element,
                                const FVElementGeometry& fvGeometry,
                                const ElementVolumeVariables& elemVolVars,
                                const SubControlVolume& scv,
                                const Scalar height,
                                const Scalar factor = 8.0) const
    {
        const Scalar velocity = elemVolVars[scv].velocity();
        const auto scvf = scvfs(fvGeometry, scv).begin();
        const Scalar viscosity = effectiveViscosity(element, fvGeometry, *scvf);
        return pseudo3DWallFriction(velocity, viscosity, height, factor);
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

    /*!
     * \brief Returns true if the scvf is located on a boundary with a slip condition.
     * \note This member function must be overloaded in the problem implementation.
     *        Make sure to use scvf.center() for querying the spatial location.
     */
    bool onSlipBoundary(const FVElementGeometry& fvGeometry, const SubControlVolumeFace& scvf) const
    { return asImp_().onSlipBoundaryAtPos(scvf.center()); }

    /*!
     * \brief Returns true if the scvf is located on a boundary with a slip condition.
     * \note This member function must be overloaded in the problem implementation.
     */
    bool onSlipBoundaryAtPos(const GlobalPosition& pos) const
    { return false; }

    /*!
     * \brief Returns the intrinsic permeability of required as input parameter for the Beavers-Joseph-Saffman boundary condition
     * This member function must be overloaded in the problem implementation, if the BJS boundary condition is used.
     */
    Scalar permeability(const FVElementGeometry& fvGeometry, const SubControlVolumeFace& scvf) const
    {
        DUNE_THROW(Dune::NotImplemented,
            "When using the Beavers-Joseph-Saffman boundary condition, "
            "the permeability must be returned in the actual problem"
        );
    }

    /*!
     * \brief Returns the alpha value required as input parameter for the Beavers-Joseph-Saffman boundary condition
     * This member function must be overloaded in the problem implementation, if the BJS boundary condition is used.
     */
    Scalar alphaBJ(const FVElementGeometry& fvGeometry, const SubControlVolumeFace& scvf) const
    {
        DUNE_THROW(Dune::NotImplemented,
            "When using the Beavers-Joseph-Saffman boundary condition, "
            "the alpha value must be returned in the actual problem"
        );
    }

    /*!
     * \brief Returns the beta value which is the alpha value divided by the square root of the (scalar-valued) interface permeability.
     */
    Scalar betaBJ(const FVElementGeometry& fvGeometry, const SubControlVolumeFace& scvf, const GlobalPosition& tangentialVector) const
    {
        const Scalar interfacePermeability = interfacePermeability_(fvGeometry, scvf, tangentialVector);
        using std::sqrt;
        return asImp_().alphaBJ(fvGeometry, scvf) / sqrt(interfacePermeability);
    }

    /*!
     * \brief Returns the velocity in the porous medium (which is 0 by default according to Saffmann).
     */
    VelocityVector porousMediumVelocity(const FVElementGeometry& fvGeometry, const SubControlVolumeFace& scvf) const
    {
        return VelocityVector(0.0);
    }

    /*!
     * \brief Returns the slip velocity at a porous boundary based on the Beavers-Joseph(-Saffman) condition.
     * \note This only returns a vector filled with one component of the slip velocity (corresponding to the dof axis of the scv the svf belongs to)
     */
    const VelocityVector beaversJosephVelocity(const FVElementGeometry& fvGeometry,
                                               const SubControlVolumeFace& scvf,
                                               const ElementVolumeVariables& elemVolVars,
                                               Scalar tangentialVelocityGradient) const
    {
        assert(scvf.isLateral());
        assert(scvf.boundary());

        const auto& scv = fvGeometry.scv(scvf.insideScvIdx());

        // create a unit normal vector oriented in positive coordinate direction
        GlobalPosition orientation(0.0);
        orientation[scv.dofAxis()] = 1.0;

        // du/dy + dv/dx = beta * (u_boundary-uPM)
        // beta = alpha/sqrt(K)
        const Scalar betaBJ = asImp_().betaBJ(fvGeometry, scvf, orientation);
        const Scalar distanceNormalToBoundary = (scvf.ipGlobal() - scv.dofPosition()).two_norm();

        static const bool onlyNormalGradient = getParamFromGroup<bool>(this->paramGroup(), "FreeFlow.EnableUnsymmetrizedVelocityGradientForBeaversJoseph", false);
        if (onlyNormalGradient)
            tangentialVelocityGradient = 0.0;

        const Scalar scalarSlipVelocity = (tangentialVelocityGradient*distanceNormalToBoundary
            + asImp_().porousMediumVelocity(fvGeometry, scvf) * orientation * betaBJ * distanceNormalToBoundary
            + elemVolVars[scv].velocity()) / (betaBJ*distanceNormalToBoundary + 1.0);

        orientation *= scalarSlipVelocity;
        return orientation;
    }

    const CouplingManager& couplingManager() const
    {
        if constexpr (isCoupled_)
            return *couplingManager_;
        else
            DUNE_THROW(Dune::InvalidStateException,
                "Accessing coupling manager of an uncoupled problem is not possible."
            );
    }

private:
    //! Returns a scalar permeability value at the coupling interface
    template<class Scvf>
    Scalar interfacePermeability_(const FVElementGeometry& fvGeometry, const Scvf& scvf, const GlobalPosition& tangentialVector) const
    {
        const auto& K = asImp_().permeability(fvGeometry, scvf);

        // use t*K*t for permeability tensors
        if constexpr (Dune::IsNumber<std::decay_t<decltype(K)>>::value)
            return K;
        else
            return vtmv(tangentialVector, K, tangentialVector);
    }

    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation& asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation& asImp_() const
    { return *static_cast<const Implementation *>(this); }

    GravityVector gravity_;
    bool enableInertiaTerms_;
    std::shared_ptr<CouplingManager> couplingManager_;
};

/*!
 * \brief Navier-Stokes default problem implementation for FCDiamond discretizations
 */
template<class TypeTag>
class NavierStokesProblemImpl<TypeTag, DiscretizationMethods::FCDiamond>
: public FVProblemWithSpatialParams<TypeTag>
{
    using ParentType = FVProblemWithSpatialParams<TypeTag>;
    using Implementation = GetPropType<TypeTag, Properties::Problem>;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using GridVolumeVariables = typename GridVariables::GridVolumeVariables;
    using ElementVolumeVariables = typename GridVolumeVariables::LocalView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;

    using GlobalPosition = typename SubControlVolumeFace::GlobalPosition;
    using VelocityVector = Dune::FieldVector<Scalar, dimWorld>;
    using GravityVector = Dune::FieldVector<Scalar, dimWorld>;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;

public:
    //! These types are used in place of the typical NumEqVector type.
    //! In the numeqvector assembly type, only one equation per DOF (face) is considered
    //! while the type here provides one entry for each world dimension.
    using InitialValues = Dune::FieldVector<Scalar, dimWorld>;
    using Sources = Dune::FieldVector<Scalar, dimWorld>;
    using DirichletValues = Dune::FieldVector<Scalar, dimWorld>;
    using BoundaryFluxes = Dune::FieldVector<Scalar, dimWorld>;

    //! Export the boundary types.
    using BoundaryTypes = NavierStokesMomentumBoundaryTypes<ModelTraits::dim()>;

    //! This problem is used for the momentum balance model.
    static constexpr bool isMomentumProblem() { return true; }

    /*!
     * \brief The constructor
     * \param gridGeometry The finite volume grid geometry
     * \param couplingManager A coupling manager (in case this problem is a coupled "mass-momentum"/full Navier-Stokes problem)
     * \param paramGroup The parameter group in which to look for runtime parameters first (default is "")
     */
    NavierStokesProblemImpl(std::shared_ptr<const GridGeometry> gridGeometry,
                            std::shared_ptr<CouplingManager> couplingManager,
                            const std::string& paramGroup = "")
    : ParentType(gridGeometry, paramGroup)
    , gravity_(0.0)
    , couplingManager_(couplingManager)
    {
        if (getParamFromGroup<bool>(paramGroup, "Problem.EnableGravity"))
            gravity_[dim-1]  = -9.81;

        enableInertiaTerms_ = getParamFromGroup<bool>(paramGroup, "Problem.EnableInertiaTerms");
    }

    /*!
     * \brief The constructor for usage without a coupling manager
     * \param gridGeometry The finite volume grid geometry
     * \param paramGroup The parameter group in which to look for runtime parameters first (default is "")
     */
    NavierStokesProblemImpl(std::shared_ptr<const GridGeometry> gridGeometry,
                            const std::string& paramGroup = "")
    : NavierStokesProblemImpl(gridGeometry, {}, paramGroup)
    {}

    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     *
     * This is the method for the case where the source term is
     * potentially solution dependent and requires some quantities that
     * are specific to the fully-implicit method.
     *
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars All volume variables for the element
     * \param scv The sub control volume
     *
     * For this method, the return parameter stores the conserved quantity rate
     * generated or annihilate per volume unit. Positive values mean
     * that the conserved quantity is created, negative ones mean that it vanishes.
     * E.g. for the mass balance that would be a mass rate in \f$ [ kg / (m^3 \cdot s)] \f$.
     */
    template<class ElementVolumeVariables>
    Sources source(const Element &element,
                   const FVElementGeometry& fvGeometry,
                   const ElementVolumeVariables& elemVolVars,
                   const SubControlVolume &scv) const
    {
        // forward to solution independent, fully-implicit specific interface
        return asImp_().sourceAtPos(scv.center());
    }

    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     *
     * \param globalPos The position of the center of the finite volume
     *            for which the source term ought to be
     *            specified in global coordinates
     *
     * For this method, the values parameter stores the conserved quantity rate
     * generated or annihilate per volume unit. Positive values mean
     * that the conserved quantity is created, negative ones mean that it vanishes.
     * E.g. for the mass balance that would be a mass rate in \f$ [ kg / (m^3 \cdot s)] \f$.
     */
    Sources sourceAtPos(const GlobalPosition &globalPos) const
    {
        //! As a default, i.e. if the user's problem does not overload any source method
        //! return 0.0 (no source terms)
        return Sources(0.0);
    }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param element The finite element
     * \param scvf The sub control volume face
     */
    BoundaryTypes boundaryTypes(const Element &element,
                                const SubControlVolumeFace &scvf) const
    {
        // Forward it to the method which only takes the global coordinate.
        // We evaluate the boundary type at the center of the sub control volume face
        // in order to avoid ambiguities at domain corners.
        return asImp_().boundaryTypesAtPos(scvf.center());
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume face.
     *
     * \param element The finite element
     * \param scvf the sub control volume face
     * \note used for cell-centered discretization schemes
     */
    DirichletValues dirichlet(const Element &element, const SubControlVolumeFace &scvf) const
    {
        return asImp_().dirichletAtPos(scvf.ipGlobal());
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
     * \brief Returns whether inertia terms should be considered.
     */
    bool enableInertiaTerms() const
    { return enableInertiaTerms_; }

    /*!
     * \brief Returns a reference pressure
     *        This pressure is subtracted from the actual pressure for the momentum balance
     *        which potentially helps to improve numerical accuracy by avoiding issues related do floating point arithmetic.
     * \note  Overload this for reference pressures other than zero.
     */
    Scalar referencePressure() const
    { return 0.0; }

    /*!
     * \brief Returns the pressure at a given sub control volume face.
     * \note  Overload this if a fixed pressure shall be prescribed (e.g., given by an analytical solution).
     */
    Scalar pressure(const Element& element,
                    const FVElementGeometry& fvGeometry,
                    const SubControlVolumeFace& scvf) const
    {
        if constexpr (std::is_empty_v<CouplingManager>)
            return asImp_().pressureAtPos(scvf.ipGlobal());
        else
            return couplingManager_->pressure(element, fvGeometry, scvf);
    }

    /*!
     * \brief Returns the pressure at a given sub control volume.
     * \note  Overload this if a fixed pressure shall be prescribed (e.g., given by an analytical solution).
     */
    Scalar pressure(const Element& element,
                    const FVElementGeometry& fvGeometry,
                    const SubControlVolume& scv,
                    const bool isPreviousTimeStep = false) const
    {
        if constexpr (std::is_empty_v<CouplingManager>)
            return asImp_().pressureAtPos(scv.dofPosition());
        else
            return couplingManager_->pressure(element, fvGeometry, scv, isPreviousTimeStep);
    }

    /*!
     * \brief Returns the pressure at a given position.
     */
    Scalar pressureAtPos(const GlobalPosition&) const
    {
        DUNE_THROW(Dune::NotImplemented, "pressureAtPos not implemented");
    }

    /*!
     * \brief Returns the density at a given sub control volume face.
     * \note  Overload this if a fixed density shall be prescribed.
     */
    Scalar density(const Element& element,
                   const FVElementGeometry& fvGeometry,
                   const SubControlVolumeFace& scvf) const
    {
        if constexpr (std::is_empty_v<CouplingManager>)
            return asImp_().densityAtPos(scvf.ipGlobal());
        else
            return couplingManager_->density(element, fvGeometry, scvf);
    }

    /*!
     * \brief Returns the density at a given sub control volume.
     * \note  Overload this if a fixed density shall be prescribed.
     */
    Scalar density(const Element& element,
                   const FVElementGeometry& fvGeometry,
                   const SubControlVolume& scv,
                   const bool isPreviousTimeStep = false) const
    {
        if constexpr (std::is_empty_v<CouplingManager>)
            return asImp_().densityAtPos(scv.dofPosition());
        else
            return couplingManager_->density(element, fvGeometry, scv, isPreviousTimeStep);
    }


    /*!
     * \brief Returns the density at a given position.
     */
    Scalar densityAtPos(const GlobalPosition&) const
    {
        DUNE_THROW(Dune::NotImplemented, "densityAtPos not implemented");
    }

    /*!
     * \brief Returns the effective dynamic viscosity at a given sub control volume face.
     * \note  Overload this if a fixed viscosity shall be prescribed.
     */
    Scalar effectiveViscosity(const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const SubControlVolumeFace& scvf) const
    {
        if constexpr (std::is_empty_v<CouplingManager>)
            return asImp_().effectiveViscosityAtPos(scvf.ipGlobal());
        else
            return couplingManager_->effectiveViscosity(element, fvGeometry, scvf);
    }

    /*!
     * \brief Returns the effective dynamic viscosity at a given sub control volume.
     * \note  Overload this if a fixed viscosity shall be prescribed.
     */
    Scalar effectiveViscosity(const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const SubControlVolume& scv) const
    {
        if constexpr (std::is_empty_v<CouplingManager>)
            return asImp_().effectiveViscosityAtPos(scv.dofPosition());
        else
            return couplingManager_->effectiveViscosity(element, fvGeometry, scv);
    }

    /*!
     * \brief Returns the effective dynamic viscosity at a given position.
     */
    Scalar effectiveViscosityAtPos(const GlobalPosition&) const
    {
        DUNE_THROW(Dune::NotImplemented, "effectiveViscosityAtPos not implemented");
    }

    /*!
     * \brief Applies the initial solution for all degrees of freedom of the grid.
     * \param sol the initial solution vector
     */
    template<class SolutionVector>
    void applyInitialSolution(SolutionVector& sol) const
    {
        static_assert(GridGeometry::discMethod == DiscretizationMethods::fcdiamond);
        sol.resize(this->gridGeometry().numDofs());
        std::vector<bool> dofHandled(this->gridGeometry().numDofs(), false);
        auto fvGeometry = localView(this->gridGeometry());
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            fvGeometry.bindElement(element);
            for (const auto& scv : scvs(fvGeometry))
            {
                const auto dofIdx = scv.dofIndex();
                if (!dofHandled[dofIdx])
                {
                    dofHandled[dofIdx] = true;
                    sol[dofIdx] = asImp_().initial(scv);
                }
            }
        }
    }

    /*!
     * \brief Evaluate the initial value at an sub control volume
     */
    InitialValues initial(const SubControlVolume& scv) const
    {
        static_assert(GridGeometry::discMethod == DiscretizationMethods::fcdiamond);
        return asImp_().initialAtPos(scv.dofPosition());
    }

private:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

    GravityVector gravity_;
    bool enableInertiaTerms_;
    std::shared_ptr<CouplingManager> couplingManager_ = {};
};

template<class TypeTag>
class NavierStokesProblemImpl<TypeTag, DiscretizationMethods::CCTpfa>
: public FVProblemWithSpatialParams<TypeTag>
{
    using ParentType = FVProblemWithSpatialParams<TypeTag>;
    using Implementation = GetPropType<TypeTag, Properties::Problem>;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename SubControlVolumeFace::GlobalPosition;
    using VelocityVector = GlobalPosition;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;

    static constexpr bool isCoupled_ = !std::is_empty_v<CouplingManager>;

public:

    //! These types are used in place of the typical NumEqVector type.
    //! In the momentum problem, there is a descrepancy between
    //! the assembly NumEqVector type, and the type used in the residual.
    //! These aliases are used in both problems to differentiate.
    using InitialValues = Dune::FieldVector<Scalar, ModelTraits::numEq()>;
    using Sources = Dune::FieldVector<Scalar, ModelTraits::numEq()>;
    using DirichletValues = Dune::FieldVector<Scalar, ModelTraits::numEq()>;
    using BoundaryFluxes = Dune::FieldVector<Scalar, ModelTraits::numEq()>;

    //! Export the boundary types.
    using BoundaryTypes = Dumux::BoundaryTypes<ModelTraits::numEq()>;

    //! this problem is used for the mass balance model
    static constexpr bool isMomentumProblem() { return false; }

    /*!
     * \brief The constructor
     * \param gridGeometry The finite volume grid geometry
     * \param couplingManager The coupling manager (couples mass and momentum equations)
     * \param paramGroup The parameter group in which to look for runtime parameters first (default is "")
     */
    NavierStokesProblemImpl(std::shared_ptr<const GridGeometry> gridGeometry,
                            std::shared_ptr<CouplingManager> couplingManager,
                            const std::string& paramGroup = "")
    : ParentType(gridGeometry, paramGroup)
    , couplingManager_(couplingManager)
    {}

    /*!
     * \brief The constructor for usage without a coupling manager
     * \param gridGeometry The finite volume grid geometry
     * \param paramGroup The parameter group in which to look for runtime parameters first (default is "")
     */
    NavierStokesProblemImpl(std::shared_ptr<const GridGeometry> gridGeometry,
                            const std::string& paramGroup = "")
    : NavierStokesProblemImpl(gridGeometry, {}, paramGroup)
    {}

    /*!
     * \brief Returns the normal velocity at a given sub control volume face.
     */
    VelocityVector faceVelocity(const Element& element,
                                const FVElementGeometry& fvGeometry,
                                const SubControlVolumeFace& scvf) const
    {
        if constexpr (isCoupled_)
            return couplingManager_->faceVelocity(element, scvf);
        else
            return asImp_().velocityAtPos(scvf.ipGlobal());
    }

    /*!
     * \brief Returns the velocity at a given position.
     */
    VelocityVector velocityAtPos(const GlobalPosition&) const
    {
        DUNE_THROW(Dune::NotImplemented, "velocityAtPos not implemented");
    }

    const CouplingManager& couplingManager() const
    {
        if constexpr (isCoupled_)
            return *couplingManager_;
        else
            DUNE_THROW(Dune::InvalidStateException,
                "Accessing coupling manager of an uncoupled problem is not possible."
            );
    }

private:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

    std::shared_ptr<CouplingManager> couplingManager_;
};

/*!
 * \ingroup NavierStokesModel
 * \brief Navier-Stokes problem base class.
 *
 * This implements gravity (if desired) and a function returning the temperature.
 * Includes a specialized method used only by the staggered grid discretization.
 */
template<class TypeTag>
class NavierStokesProblemImpl<TypeTag, DiscretizationMethods::Staggered>
: public NavierStokesParentProblem<TypeTag>
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
    NavierStokesProblemImpl(std::shared_ptr<const GridGeometry> gridGeometry, const std::string& paramGroup = "")
    : ParentType(gridGeometry, paramGroup)
    , gravity_(0.0)
    {
        if (getParamFromGroup<bool>(paramGroup, "Problem.EnableGravity"))
            gravity_[dim-1]  = -9.81;

        enableInertiaTerms_ = getParamFromGroup<bool>(paramGroup, "Problem.EnableInertiaTerms");
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

    //! Applies the initial face solution (velocities on the faces). Specialization for staggered grid discretization.
    template <class SolutionVector, class G = GridGeometry>
    typename std::enable_if<G::discMethod == DiscretizationMethods::staggered, void>::type
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
    typename std::enable_if<G::discMethod == DiscretizationMethods::staggered, Scalar>::type
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
        DUNE_THROW(Dune::NotImplemented,
            "When using the Beavers-Joseph-Saffman boundary condition, "
            "the permeability must be returned in the actual problem"
        );
    }

    /*!
     * \brief Returns the alpha value required as input parameter for the Beavers-Joseph-Saffman boundary condition
     *
     * This member function must be overloaded in the problem implementation, if the BJS boundary condition is used.
     */
    Scalar alphaBJ(const SubControlVolumeFace& scvf) const
    {
        DUNE_THROW(Dune::NotImplemented,
            "When using the Beavers-Joseph-Saffman boundary condition, "
            "the alpha value must be returned in the actual problem"
        );
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


/*!
 * \ingroup NavierStokesModel
 * \brief Navier-Stokes problem class
 *
 * Inherit from this problem to implement Navier-Stokes problems
 */
template<class TypeTag>
using NavierStokesProblem = NavierStokesProblemImpl<
    TypeTag, typename GetPropType<TypeTag, Properties::GridGeometry>::DiscretizationMethod
>;

} // end namespace Dumux

#endif
