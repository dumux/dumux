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
 *
 * \brief A simple Darcy test problem (cell-centered finite volume method).
 */
#ifndef DUMUX_DARCY_SUBPROBLEM_HH
#define DUMUX_DARCY_SUBPROBLEM_HH

#include <dune/grid/yaspgrid.hh>
#include <dune/alugrid/grid.hh>

#include <dumux/discretization/cellcentered/mpfa/properties.hh>

#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/porousmediumflow/problem.hh>

#include "1pspatialparams.hh"

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include <dumux/multidomain/boundary/stokesdarcy/mpfa/upwindscheme.hh>

#ifndef DARCYGRIDTYPE
#define DARCYGRIDTYPE Dune::YaspGrid<2>
#endif

namespace Dumux
{
template <class TypeTag>
class DarcySubProblem;

namespace Properties
{
NEW_TYPE_TAG(DarcyOneP, INHERITS_FROM(CCMpfaModel, OneP));

// Set the problem property
SET_TYPE_PROP(DarcyOneP, Problem, Dumux::DarcySubProblem<TypeTag>);

// the fluid system
SET_PROP(DarcyOneP, FluidSystem)
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};

// Set the grid type
SET_TYPE_PROP(DarcyOneP, Grid, DARCYGRIDTYPE);

SET_PROP(DarcyOneP, SpatialParams)
{
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using type = OnePConvSpatialParams<FVGridGeometry, Scalar>;
};

//! The flux variables for models involving flow in porous media
SET_TYPE_PROP(DarcyOneP,
              FluxVariables,
              PorousMediumFluxVariables<TypeTag, CCMpfaStokesDarcyUpwindScheme<typename GET_PROP_TYPE(TypeTag, FVGridGeometry)> >);
}

template <class TypeTag>
class DarcySubProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);

    using Indices = typename GET_PROP_TYPE(TypeTag, ModelTraits)::Indices;

    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    using CouplingManager = typename GET_PROP_TYPE(TypeTag, CouplingManager);

    DarcySubProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry,
                   std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(fvGridGeometry, "Darcy"), eps_(1e-7), couplingManager_(couplingManager)
    {
        createAnalyticalSolution_();
    }

    /*!
     * \name Simulation steering
     */
    // \{

    /*!
     * \brief Returns true if a restart file should be written to
     *        disk.
     */
    bool shouldWriteRestartFile() const
    { return false; }

    /*!
     * \name Problem parameters
     */
    // \{

    bool shouldWriteOutput() const // define output
    { return true; }

    /*!
     * \brief Return the temperature within the domain in [K].
     *
     */
    Scalar temperature() const
    { return 273.15 + 10; } // 10°C
    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
      * \brief Specifies which kind of boundary condition should be
      *        used for which equation on a given boundary control volume.
      *
      * \param element The element
      * \param scvf The boundary sub control volume face
      */
    BoundaryTypes boundaryTypes(const Element &element, const SubControlVolumeFace &scvf) const
    {
        BoundaryTypes values;
        values.setAllDirichlet();

        if (couplingManager().isCoupledEntity(CouplingManager::darcyIdx, scvf))
            values.setAllCouplingNeumann();

        return values;
    }

        /*!
     * \brief Evaluate the boundary conditions for a Dirichlet control volume.
     *
     * \param element The element for which the Dirichlet boundary condition is set
     * \param scvf The boundary subcontrolvolumeface
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        return exactPressure(globalPos);
    }

    /*!
     * \brief Evaluate the boundary conditions for a Neumann control volume.
     *
     * \param element The element for which the Neumann boundary condition is set
     * \param fvGeomentry The fvGeometry
     * \param elemVolVars The element volume variables
     * \param scvf The boundary sub control volume face
     *
     * For this method, the \a values variable stores primary variables.
     */
    template<class ElementVolumeVariables>
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const SubControlVolumeFace& scvf) const
    {
        NumEqVector values(0.0);

        if (couplingManager().isCoupledEntity(CouplingManager::darcyIdx, scvf))
            values[Indices::conti0EqIdx] = couplingManager().couplingData().massCouplingCondition(fvGeometry, elemVolVars, scvf);

        return values;
    }

    /*!
     * \brief Function to add the coupling Neumann boundary condition
     *        into the local system of equations in the mpfa interaction volume.
     *
     * \param element The element for which the Neumann boundary condition is set
     * \param fvGeomentry The fvGeometry
     * \param elemVolVars The element volume variables
     * \param scvf The boundary sub control volume face
     */
    template<class ElementVolumeVariables>
    Scalar neumannTerm(const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolumeFace& scvf) const
    {
        NumEqVector values(0.0);

        if (couplingManager().isCoupledEntity(CouplingManager::darcyIdx, scvf))
            values[Indices::conti0EqIdx] = scvf.area()*couplingManager().couplingData().neumannCouplingCondition(element, fvGeometry, elemVolVars, scvf);

        return values;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
      * \brief Return the sources within the domain.
      *
      * \param globalPos The global position
      */
     NumEqVector sourceAtPos(const GlobalPosition &globalPos) const
     {
        NumEqVector source(0.0);
        const Scalar x = globalPos[0];
        const Scalar y = globalPos[1];

        using std::cos;
        using std::sin;
        using std::exp;

        source[Indices::conti0EqIdx] = (4*M_PI*M_PI*(exp(y+1) + 2 - exp(2)) - exp(y-1))*sin(2*M_PI*x);

        return source;
    }

    // \}

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param element The element
     *
     * For this method, the \a priVars parameter stores primary
     * variables.
     */
    PrimaryVariables initial(const Element &element) const
    {
        return PrimaryVariables(exactPressure(element.geometry().center()));
    }

    // \}

    //! Set the coupling manager
    void setCouplingManager(std::shared_ptr<CouplingManager> cm)
    { couplingManager_ = cm; }

    //! Get the coupling manager
    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

    /*!
      * \brief Returns the analytical solution for the pressure
      */
     auto& getAnalyticalPressureSolution() const
     {
         return analyticalPressure_;
     }

    /*!
      * \brief Returns the analytical solution for the velocity
      */
     auto& getAnalyticalVelocitySolution() const
     {
         return analyticalVelocity_;
     }


private:

    Scalar exactPressure(const GlobalPosition &globalPos) const
    {
        using std::exp;
        using std::sin;
        Scalar x = globalPos[0];
        Scalar y = globalPos[1];
        return (exp(y+1) + 2 - exp(2))*sin(2*M_PI*x);
    }

    GlobalPosition exactVelocity(const GlobalPosition &globalPos) const
    {
        using std::sin;
        using std::cos;
        using std::exp;
        Scalar x = globalPos[0];
        Scalar y = globalPos[1];
        GlobalPosition velocity(0.0);
        velocity[0] = -2*M_PI*(exp(y+1) + 2 - exp(2))*cos(2*M_PI*x);
        velocity[1] = -exp(y-1)*sin(2*M_PI*x);
        return velocity;
    }

    /*!
      * \brief Adds additional VTK output data to the VTKWriter. Function is called by the output module on every write.
      */
     void createAnalyticalSolution_()
     {
         analyticalPressure_.resize(this->fvGridGeometry().numDofs());
         analyticalVelocity_.resize(this->fvGridGeometry().numDofs());

         for (const auto& element : elements(this->fvGridGeometry().gridView()))
         {
             auto fvGeometry = localView(this->fvGridGeometry());
             fvGeometry.bindElement(element);
             for (auto&& scv : scvs(fvGeometry))
             {
                 auto dofIdx = scv.dofIndex();
                 auto dofPosition = scv.dofPosition();

                 analyticalPressure_[dofIdx] = exactPressure(dofPosition);
                 analyticalVelocity_[dofIdx] = exactVelocity(dofPosition);
             }
         }
      }

    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] < this->fvGridGeometry().bBoxMin()[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] > this->fvGridGeometry().bBoxMax()[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] < this->fvGridGeometry().bBoxMin()[1] + eps_; }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] > this->fvGridGeometry().bBoxMax()[1] - eps_; }

    Scalar eps_;
    std::shared_ptr<CouplingManager> couplingManager_;
    std::vector<Scalar> analyticalPressure_;
    std::vector<GlobalPosition> analyticalVelocity_;
};
} //end namespace

#endif //DUMUX_DARCY_SUBPROBLEM_HH
