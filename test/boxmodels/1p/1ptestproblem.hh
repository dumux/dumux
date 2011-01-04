// $Id$
/*****************************************************************************
 *   Copyright (C) 2009 by Onur Dogan                                        *
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief A test problem for the one-phase box model
 */
#ifndef DUMUX_1PTEST_PROBLEM_HH
#define DUMUX_1PTEST_PROBLEM_HH

#ifdef HAVE_UG
#include <dune/grid/io/file/dgfparser/dgfug.hh>
#endif
#ifdef HAVE_ALUGRID
#include <dune/grid/io/file/dgfparser/dgfalu.hh>
#endif
#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#include <dune/grid/io/file/gmshreader.hh>

#include <dumux/boxmodels/1p/1pmodel.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>

#include "1ptestspatialparameters.hh"

namespace Dumux
{
template <class TypeTag>
class OnePTestProblem;

namespace Properties
{
NEW_TYPE_TAG(OnePTestProblem, INHERITS_FROM(BoxOneP));

SET_PROP(OnePTestProblem, Fluid)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::SimpleH2O<Scalar> > type;
};

// Set the grid type
SET_PROP(OnePTestProblem, Grid)
{
    typedef Dune::SGrid<2, 2> type;
};

#if HAVE_DUNE_PDELAB
SET_PROP(OnePTestProblem, LocalFEMSpace)
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    enum{dim = GridView::dimension};

public:
    typedef Dune::PDELab::Q1LocalFiniteElementMap<Scalar,Scalar,dim> type; // for cubes
//    typedef Dune::PDELab::P1LocalFiniteElementMap<Scalar,Scalar,dim> type; // for simplices
};
#endif // HAVE_DUNE_PDELAB

SET_INT_PROP(OnePTestProblem, NewtonLinearSolverVerbosity, 0);

// Set the problem property
SET_PROP(OnePTestProblem, Problem)
{ typedef Dumux::OnePTestProblem<TTAG(OnePTestProblem)> type; };

// Set the spatial parameters
SET_PROP(OnePTestProblem, SpatialParameters)
{ typedef Dumux::OnePTestSpatialParameters<TypeTag> type; };

// Enable gravity
SET_BOOL_PROP(OnePTestProblem, EnableGravity, true);
}

/*!
 * \ingroup OnePBoxModel
 * \brief Air flow in porous media
 *
 * The domain is box shaped. All sides are closed (Neumann 0 boundary)
 * except the top and bottom boundaries (Dirichlet), where water is
 * flowing from bottom to top.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_1p grids/1p_2d.dgf 10 0.01</tt>
 * where start simulation time = 0.01 second, end simulation time = 10 seconds
 * The same file can be also used for 3d simulation but you need to change line
 * <tt>typedef Dune::SGrid<2,2> type;</tt> to
 * <tt>typedef Dune::SGrid<3,3> type;</tt> and use <tt>1p_3d.dgf</tt> grid.
 */
template <class TypeTag = TTAG(OnePTestProblem) >
class OnePTestProblem : public OnePBoxProblem<TypeTag>
{
    typedef OnePTestProblem<TypeTag> ThisType;
    typedef OnePBoxProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(OnePIndices)) Indices;
    enum {
        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,

        // indices of the primary variables
        pressureIdx = Indices::pressureIdx,
    };


    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(BoundaryTypes)) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TimeManager)) TimeManager;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::Intersection Intersection;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;

    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
    OnePTestProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
    {
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const char *name() const
    { return "1ptest"; }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature(const Element &element,
                       const FVElementGeometry &fvElemGeom,
                       int scvIdx) const
    {
        return 273.15 + 10; // 10°C
    };

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     */
    void boundaryTypes(BoundaryTypes &values, const Vertex &vertex) const
    {
        const GlobalPosition globalPos = vertex.geometry().center();

        double eps = 1.0e-3;
        if (globalPos[dim-1] < eps || globalPos[dim-1] > this->bboxMax()[dim-1] - eps)
            values.setAllDirichlet();
        else
            values.setAllNeumann();
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        boundary segment.
     *
     * For this method, the \a values parameter stores primary variables.
     */
    void dirichlet(PrimaryVariables &values, const Vertex &vertex) const
    {
        double eps = 1.0e-3;
        const GlobalPosition globalPos = vertex.geometry().center();

        if (globalPos[dim-1] < eps) {
            values[pressureIdx] = 2.0e+5;
        }
        else if (globalPos[dim-1] > this->bboxMax()[dim-1] - eps) {
            values[pressureIdx] = 1.0e+5;
        }
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each component. Negative values mean
     * influx.
     */
    void neumann(PrimaryVariables &values,
                 const Element &element,
                 const FVElementGeometry &fvElemGeom,
                 const Intersection &is,
                 int scvIdx,
                 int boundaryFaceIdx) const
    {
        //  const GlobalPosition &globalPos = fvElemGeom.boundaryFace[boundaryFaceIdx].ipGlobal;

        values[pressureIdx] = 0;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     *
     * For this method, the \a values parameter stores the rate mass
     * of a component is generated or annihilate per volume
     * unit. Positive values mean that mass is created, negative ones
     * mean that it vanishes.
     */
    void source(PrimaryVariables &values,
                const Element &element,
                const FVElementGeometry &fvElemGeom,
                int scvIdx) const
    {
        values = Scalar(0.0);
    }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    void initial(PrimaryVariables &values,
                 const Element &element,
                 const FVElementGeometry &fvElemGeom,
                 int scvIdx) const
    {
        //const GlobalPosition &globalPos = element.geometry().corner(scvIdx);
        values[pressureIdx] = 1.0e+5;// + 9.81*1.23*(20-globalPos[dim-1]);
    }

    // \}
};
} //end namespace

#endif
