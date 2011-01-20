// $Id$
/*****************************************************************************
 *   Copyright (C) 2009 by Karin Erbertseder                                 *
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
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
/**
 * \file
 * \brief Definition of a problem, where the distribution of a therapeutic agent
 * within pulmonary tissue is described
 * \author Karin Erbertseder, Bernd Flemisch
 */
#ifndef DUMUX_TISSUE_TUMOR_PROBLEM_HH
#define DUMUX_TISSUE_TUMOR_PROBLEM_HH

#ifdef HAVE_UG
#include <dune/grid/io/file/dgfparser/dgfug.hh>
#endif
#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

#include <dumux/material/fluidsystems/isfluid_trail_system.hh>
#include <dumux/boxmodels/1p2c/1p2cmodel.hh>

#include "tissue_tumor_spatialparameters.hh"

namespace Dumux
{

template <class TypeTag>
class TissueTumorProblem;

namespace Properties
{
NEW_TYPE_TAG(TissueTumorProblem, INHERITS_FROM(BoxOnePTwoC));

// Set the grid type
SET_PROP(TissueTumorProblem, Grid)
{
#if HAVE_UG
    typedef Dune::UGGrid<2> type;
#else
    typedef Dune::SGrid<2, 2> type;
    //typedef Dune::YaspGrid<2> type;
#endif
};

#if HAVE_DUNE_PDELAB
SET_PROP(TissueTumorProblem, LocalFEMSpace)
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    enum{dim = GridView::dimension};

public:
    typedef Dune::PDELab::Q1LocalFiniteElementMap<Scalar,Scalar,dim> type; // for cubes
//    typedef Dune::PDELab::P1LocalFiniteElementMap<Scalar,Scalar,dim> type; // for simplices
};
#endif // HAVE_DUNE_PDELAB


// Set the problem property
SET_PROP(TissueTumorProblem, Problem)
{
    typedef Dumux::TissueTumorProblem<TTAG(TissueTumorProblem)> type;
};

// Set fluid configuration
SET_PROP(TissueTumorProblem, FluidSystem)
{
    typedef Dumux::ISFluid_Trail_System<TypeTag> type;
};

// Set the spatial parameters
SET_TYPE_PROP(TissueTumorProblem,
              SpatialParameters,
              Dumux::TissueTumorSpatialParameters<TypeTag>);


// Disable gravity
SET_BOOL_PROP(TissueTumorProblem, EnableGravity, false);
}


/*!
 * \ingroup OnePTwoCBoxModel
 *
 * \brief Definition of a problem, where the distribution of a therapeutic agent
 *         within pulmonary tissue is described
 *
 * The model domain is 22 mm long in x-direction and in y-direction with a discretization length of 0.1
 * mm. The tumour area is located in the middle of the model domain. The diameter of the tumour is
 * assumed to be 2 mm.
 *
 * The intercapillary distance is in the range of 0.1 mm. So the distance between the grid nodes is
 * equal to the intercapillary distance. It is assumed that at each node within the model domain the
 * transition of the therapeutic agent from the blood capillary into the tissue can take place.
 * Based on this assumption, the initial conditions are adapted. The mole fraction of the dissolved
 * therapeutic agent x is set to 1.1249 e-8 within the normal pulmonary tissue. Within the tumour the
 * mole fraction of dissolved therapeutic agent is set to zero due to the assumption of a blood vessel
 * free tumour. As initial condition for the pressure p, an interstitial fluid pressure of -1067 Pa is
 * assumed. This value corresponds to the interstitial fluid pressure in healthy pulmonary tissue.
 *
 * All four sides of the model domain are described by Dirichlet boundary conditions. The primary
 * variable p is instantiated with the value -1067 Pa and the primary variable x with the value
 * 1.1249 e-8.
 *
 * The pressure field of the tumour is generated by an additional source term of 1.98 e-9 l/h,
 * that is set in the center of the tumour region.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_1p2c grids/test_1p2c.dgf 1 1</tt>
 */

template <class TypeTag = TTAG(TissueTumorProblem) >
class TissueTumorProblem : public OnePTwoCBoxProblem<TypeTag>
{
    typedef TissueTumorProblem<TypeTag> ThisType;
    typedef OnePTwoCBoxProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(PrimaryVariables)) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(BoundaryTypes)) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TimeManager)) TimeManager;

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(OnePTwoCIndices)) Indices;
    enum {
        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,

        // indices of the primary variables
        pressureIdx = Indices::pressureIdx,
        x1Idx = Indices::x1Idx,

        // indices of the equations
        contiEqIdx = Indices::contiEqIdx,
        transEqIdx = Indices::transEqIdx,
    };


    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::Intersection Intersection;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
    TissueTumorProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
    {
        // calculate the injection volume
        totalInjectionVolume_ = 0;
        FVElementGeometry fvGeom;
        ElementIterator elemIt = gridView.template begin<0>();
        const ElementIterator endIt = gridView.template end<0>();
        for (; elemIt != endIt; ++ elemIt) {
            fvGeom.update(gridView, *elemIt);
            for (int i = 0; i < fvGeom.numVertices; ++i) {
                const GlobalPosition &pos = fvGeom.subContVol[i].global;
                if (inInjectionVolume_(pos))
                    totalInjectionVolume_ += fvGeom.subContVol[i].volume;
            };
        }
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
    { return "tissue"; }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * This problem assumes a temperature of 36 degrees Celsius.
     */
    Scalar temperature(const Element &element,
                       const FVElementGeometry &fvElemGeom,
                       int scvIdx) const
    {
        return 273.15 + 36; // in [K]
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
        values.setAllDirichlet();
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        boundary segment.
     *
     * For this method, the \a values parameter stores primary variables.
     */
    void dirichlet(PrimaryVariables &values, const Vertex &vertex) const
    {
        const GlobalPosition globalPos = vertex.geometry().center();

        initial_(values, globalPos);
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
        const GlobalPosition &globalPos
            = element.geometry().corner(scvIdx);
        values = 0;

        //int globalIdx = this->model().vertexMapper().map(element, scvIdx, dim);

        //Scalar lambda = (globalPos[1])/height_;
        if (globalPos[0] < eps_ ) {
            values[contiEqIdx] = -3.8676e-2; // [kg/(m^2 * s)]

            //values[transEqIdx] = -4.35064e-4; // [mol/(m^2*s)
            //Robin-Boundary
            //values[transEqIdx] = (*this->model().curSolFunction())[globalIdx][transEqIdx];
        }
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
        const GlobalPosition &globalPos
            = element.geometry().corner(scvIdx);

        values = Scalar(0.0);
        if (inInjectionVolume_(globalPos)) {
            // total volumetric injection rate in ml/h
            Scalar injRateVol = 0.1;
            // convert to m^3/s
            injRateVol *= 1e-6/3600;
            // total mass injection rate. assume a density of 1030kg/m^3
            Scalar injRateMass = injRateVol*1030.0;

            // trail concentration in injected fluid in [mol/ml]
            Scalar trailInjRate = 1e-5;
            // convert to mol/m^3
            trailInjRate *= 1e6;
            // convert to mol/s
            trailInjRate *= injRateVol;

            // source term of the total mass
            values[contiEqIdx] = injRateMass / totalInjectionVolume_; // [kg/(s*m^3)]
            values[transEqIdx] = trailInjRate / totalInjectionVolume_; // [mol/(s*m^3)]
        }
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
        const GlobalPosition &globalPos
            = element.geometry().corner(scvIdx);

        initial_(values, globalPos);
    }

    // \}

private:
    bool inInjectionVolume_(const GlobalPosition &globalPos) const
    {
        return
            10e-3 < globalPos[0] && globalPos[0] < 12e-3 &&
            10e-3 < globalPos[1] && globalPos[1] < 12e-3;
    };
    // the internal method for the initial condition
    void initial_(PrimaryVariables &values,
                  const GlobalPosition &globalPos) const
    {
        values[pressureIdx] = 0; //initial condition for the pressure
        values[x1Idx] = 0; //initial condition for the trail molefraction
    }

    Scalar totalInjectionVolume_;
    static const Scalar eps_ = 1e-6;
}; //end namespace
}
#endif
