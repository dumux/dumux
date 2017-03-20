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
/**
 * \file
 * \brief  Definition of a simple Navier-Stokes problem
 */
#ifndef DUMUX_TEST_FREEFLOW_NAVIERSTOKESTESTPROBLEM_HH
#define DUMUX_TEST_FREEFLOW_NAVIERSTOKESTESTPROBLEM_HH

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/gasphase.hh>

#include <dumux/freeflow/stokes/model.hh>

namespace Dumux
{
  template <class TypeTag>
  class NavierStokesTestProblem;

  // Specify the properties for the stokes problem
  namespace Properties
  {
    NEW_TYPE_TAG(NavierStokesTestProblem, INHERITS_FROM(BoxStokes));

    // Set the grid type
    SET_TYPE_PROP(NavierStokesTestProblem, Grid,
                  Dune::YaspGrid<2, Dune::TensorProductCoordinates<typename GET_PROP_TYPE(TypeTag, Scalar), 2> >);

    // Set the problem property
    SET_TYPE_PROP(NavierStokesTestProblem, Problem, NavierStokesTestProblem<TypeTag>);

    // Set calculation to Navier-Stokes, not Stokes
    SET_BOOL_PROP(NavierStokesTestProblem, EnableNavierStokes, true);

    // Use nitrogen with a constant viscosity as gas phase
    SET_TYPE_PROP(NavierStokesTestProblem, Fluid,
                  FluidSystems::GasPhase<typename GET_PROP_TYPE(TypeTag, Scalar),
                                                Constant<TypeTag, typename GET_PROP_TYPE(TypeTag, Scalar)> >);
  }

  /*!
   * \ingroup BoxStokesModel
   * \ingroup ImplicitTestProblems
   * \brief Stokes flow problem with modified nitrogen (N2) circulating in
   *        a cavity. (lid-driven cavity-flow)
   *
   * The example is taken from Ghia, Ghia, and Shin (1982), "High-Re solutions
   * for incompressible flow using the Navier-Stokes equations and a multigrid
   * method", Journal of Computational Physics, Vol. 48, pp. 387-411.
   *
   * The domain is two-dimensional and sized 1m times 1m. The boundary conditions
   * for the momentum balances are Neumann zero boundary conditions except for
   * the top, which is floating from left to right with 1 m/s. The mass balance
   * has outflow boundary conditions, which are replaced in the localresidual by
   * the sum of the two momentum balances. All vertices at the bottom receive
   * Dirichlet boundary conditions to set the pressure level.
   *
   * This problem uses the \ref BoxStokesModel with <code>EnableNavierStokes</code>
   * set to <code>true</code>.
   */
  template <class TypeTag>
  class NavierStokesTestProblem : public StokesProblem<TypeTag>
  {
    typedef StokesProblem<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum
    {
      // Number of equations and grid dimension
      dim = GridView::dimension,
      dimWorld = GridView::dimensionworld,

      // copy some indices for convenience
      massBalanceIdx = Indices::massBalanceIdx, //!< Index of the mass balance
      pressureIdx = Indices::pressureIdx, //!< Index of the pressure
      velocityXIdx = Indices::velocityXIdx, //!< Index of the x-component of the velocity
      velocityYIdx = Indices::velocityYIdx //!< Index of the y-component of the velocity
    };

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::ctype CoordScalar;
    typedef typename GridView::Intersection Intersection;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, Fluid) Fluid;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;

  public:
      NavierStokesTestProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
      {
        eps_ = 1e-6;
      }

      /**
       * \name Problem parameters
       */
      // \{

      /**
       * \brief The problem name.
       *
       * This is used as a prefix for files generated by the simulation.
       */
      std::string name() const
      {
        return "navierstokes";
      }

      /**
       * \brief Returns the temperature within the domain.
       *
       * This problem assumes a constant temperature of 10 degrees Celsius.
       */
      Scalar temperatureAtPos(const GlobalPosition &globalPos) const
      {
        return 273.15 + 10; // -> 10C
      }

      // \}

      /**
       * \name Boundary conditions
       */
      // \{

      /**
       * \brief Specifies which kind of boundary condition should be
       *        used for which equation on a given boundary segment.
       *
       * \param values The boundary types for the conservation equations
       * \param globalPos The globalPos on the boundary for which the
       *                  conditions needs to be specified
       */
      void boundaryTypesAtPos(BoundaryTypes &values,
                            const GlobalPosition &globalPos) const
      {
        values.setOutflow(massBalanceIdx);
        values.setDirichlet(velocityXIdx);
        values.setDirichlet(velocityYIdx);

        // set pressure for all vertices at the bottom
        if (onLowerBoundary_(globalPos))
        {
          values.setDirichlet(pressureIdx);
        }
      }

      /**
       * \brief Evaluate the boundary conditions for a dirichlet
       *        control volume.
       *
       * \param values The dirichlet values for the primary variables
       * \param globalPos The globalPos representing the "half volume on the boundary"
       *
       * For this method, the \a values parameter stores primary variables.
       */
      void dirichletAtPos(PrimaryVariables &values,
                        const GlobalPosition &globalPos) const
      {
        initial_(values, globalPos);

        // lid moves from left to right
        const Scalar lidVelocity = 1.0;
        if (onUpperBoundary_(globalPos))
        {
          values[Indices::velocityXIdx] = lidVelocity;
        }
      }

      /**
       * \brief Evaluate the boundary conditions for a neumann
       *        boundary segment.
       *
       * For this method, the \a values parameter stores the mass flux
       * in normal direction of each phase. Negative values mean influx.
       *
       * A Neumann condition for the Stokes equation corresponds to:
       * \f[ -\mu \nabla {\bf v} \cdot {\bf n} + p \cdot {\bf n} = q_N \f]
       */
      void neumann(PrimaryVariables &values,
                   const Element &element,
                   const FVElementGeometry &fvGeometry,
                   const Intersection &intersection,
                   const int scvIdx,
                   const int boundaryFaceIdx) const
      {
        values = 0.0;
      }
      // \}

      /**
       * \name Volume terms
       */
      // \{

      /**
       * \brief Evaluate the source term for all phases within a given
       *        sub-control-volume.
       *
       * For this method, the \a values parameter stores the rate mass
       * generated or annihilate per volume unit. Positive values mean
       * that mass is created, negative ones mean that it vanishes.
       */
      void source(PrimaryVariables &values,
                  const Element &element,
                  const FVElementGeometry &,
                  const int scvIdx) const
      {
        // ATTENTION: The source term of the mass balance has to be chosen as
        // div (q_momentum) in the problem file
        values = Scalar(0.0);
      }

      /**
       * \brief Evaluate the initial value for a control volume.
       *
       * For this method, the \a values parameter stores primary
       * variables.
       */
      void initialAtPos(PrimaryVariables &values,
                      const GlobalPosition &globalPos) const
      {
        initial_(values, globalPos);
      }
      // \}

  private:
      // internal method for the initial condition
      void initial_(PrimaryVariables &values,
                    const GlobalPosition &globalPos) const
      {
        // pressure
        values[pressureIdx] = 1e5;
        // velocity in x- and y-direction
        values[velocityXIdx] = 0.0;
        values[velocityYIdx] = 0.0;
      }

      bool onLeftBoundary_(const GlobalPosition &globalPos) const
      {
        return globalPos[0] < this->bBoxMin()[0] + eps_;
      }

      bool onRightBoundary_(const GlobalPosition &globalPos) const
      {
        return globalPos[0] > this->bBoxMax()[0] - eps_;
      }

      bool onLowerBoundary_(const GlobalPosition &globalPos) const
      {
        return globalPos[1] < this->bBoxMin()[1] + eps_;
      }

      bool onUpperBoundary_(const GlobalPosition &globalPos) const
      {
        return globalPos[1] > this->bBoxMax()[1] - eps_;
      }

      Scalar eps_;
  };
} //end namespace

#endif // DUMUX_TEST_FREEFLOW_NAVIERSTOKESTESTPROBLEM_HH
