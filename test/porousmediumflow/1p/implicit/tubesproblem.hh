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
 * \ingroup OnePTests
 * \brief A test problem for the 1p model. A pipe system with circular cross-section
 *        and a branching point embedded in a three-dimensional world
 */
#ifndef DUMUX_ONEP_TUBES_TEST_PROBLEM_HH
#define DUMUX_ONEP_TUBES_TEST_PROBLEM_HH

#include <dune/localfunctions/lagrange/pqkfactory.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dumux/discretization/cellcentered/tpfa/properties.hh>
#include <dumux/discretization/box/properties.hh>
#include <dumux/discretization/methods.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>

#include "tubesspatialparams.hh"

namespace Dumux {

template <class TypeTag>
class TubesTestProblem;

namespace Properties {

NEW_TYPE_TAG(TubesTestTypeTag, INHERITS_FROM(OneP));
NEW_TYPE_TAG(TubesTestCCTpfaTypeTag, INHERITS_FROM(CCTpfaModel, TubesTestTypeTag));
NEW_TYPE_TAG(TubesTestBoxTypeTag, INHERITS_FROM(BoxModel, TubesTestTypeTag));

// Set the grid type
SET_TYPE_PROP(TubesTestTypeTag, Grid, Dune::FoamGrid<1, 3>);

// Set the problem property
SET_TYPE_PROP(TubesTestTypeTag, Problem, TubesTestProblem<TypeTag>);

// Set the spatial parameters
SET_TYPE_PROP(TubesTestTypeTag, SpatialParams, TubesTestSpatialParams<TypeTag>);

// the fluid system
SET_PROP(TubesTestTypeTag, FluidSystem)
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using type = FluidSystems::LiquidPhase<Scalar, Components::Constant<1, Scalar> >;
};
} // end namespace Properties

/*!
 * \ingroup OnePTests
 * \brief A test problem for the 1p model. A pipe system with circular cross-section
 *        and a branching point embedded in a three-dimensional world
 */
template <class TypeTag>
class TubesTestProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);

    // Grid and world dimension
    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;

    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    enum {
        // indices of the primary variables
        conti0EqIdx = Indices::conti0EqIdx,
        pressureIdx = Indices::pressureIdx
    };

    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using Element = typename GridView::template Codim<0>::Entity;
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

    enum { isBox = GET_PROP_TYPE(TypeTag, FVGridGeometry)::discMethod == DiscretizationMethod::box };

public:
    TubesTestProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
        name_ = getParam<std::string>("Problem.Name");

        //get hMax_ of the grid
        hMax_ = 0.0;
        for (const auto& element : elements(fvGridGeometry->gridView()))
           hMax_ = std::max(element.geometry().volume(), hMax_);
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
    const std::string& name() const
    {
        return name_;
    }

    /*!
     * \brief Return the temperature within the domain in [K].
     *
     */
    Scalar temperature() const
    { return 273.15 + 37.0; } // Body temperature

    /*!
     * \brief Return how much the domain is extruded at a given sub-control volume.
     *
     * This means the factor by which a lower-dimensional (1D or 2D)
     * entity needs to be expanded to get a full dimensional cell. The
     * default is 1.0 which means that 1D problems are actually
     * thought as pipes with a cross section of 1 m^2 and 2D problems
     * are assumed to extend 1 m to the back.
     */
    template<class ElementSolution>
    Scalar extrusionFactor(const Element &element,
                           const SubControlVolume &scv,
                           const ElementSolution& elemSol) const
    {
        const auto radius = this->spatialParams().radius(scv);
        return M_PI*radius*radius;
    }

    // \}
    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary control volume.
     *
     * \param globalPos The position of the center of the finite volume
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes bcTypes;
        bcTypes.setAllDirichlet();
        return bcTypes;
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     *
     * \param globalPos The center of the finite volume which ought to be set.
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        return PrimaryVariables(0.0);
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
     * This is the method for the case where the source term is
     * potentially solution dependent and requires some quantities that
     * are specific to the fully-implicit method.
     *
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars All volume variables for the element
     * \param scv The sub control volume
     *
     * For this method, the \a values parameter stores the conserved quantity rate
     * generated or annihilate per volume unit. Positive values mean
     * that the conserved quantity is created, negative ones mean that it vanishes.
     * E.g. for the mass balance that would be a mass rate in \f$ [ kg / (m^3 \cdot s)] \f$.
     */
    NumEqVector source(const Element &element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolume &scv) const
    {
        NumEqVector source(0.0);

        const auto& globalPos = scv.center();
        const auto& volVars = elemVolVars[scv];
        const auto K = this->spatialParams().permeability(element, scv, EmptyElementSolution{});

        using std::sin;
        if (globalPos[2] > 0.5 - eps_)
            source[conti0EqIdx] = K/volVars.viscosity()*volVars.density()
                                    *4.0*4.0*M_PI*M_PI*sin(4.0*M_PI*globalPos[2]);
        else
            // the "/3.0" stems from the coordindate transformations on the lower branches
            source[conti0EqIdx] = K/volVars.viscosity()*volVars.density()
                                    *4.0*4.0*M_PI*M_PI*sin(4.0*M_PI*globalPos[2])/3.0;

        return source;
    }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * For this method, the \a priVars parameter stores primary
     * variables.
     */
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        return PrimaryVariables(0.0);
    }

    // \}

    void outputL2Norm(const SolutionVector solution) const
    {
        // calculate the discrete L2-Norm
        Scalar lTwoNorm = 0.0;

        // get the Gaussian quadrature rule for intervals
        const auto& quad = Dune::QuadratureRules<Scalar, dim>::rule(Dune::GeometryType(1), 1);

        const auto& gg = this->fvGridGeometry();
        for (const auto& element : elements(gg.gridView()))
        {
            const auto eIdx = gg.elementMapper().index(element);
            const auto geometry = element.geometry();
            for(auto&& qp : quad)
            {
                const auto globalPos = geometry.global(qp.position());
                const Scalar pe = exactPressure_(globalPos);
                const Scalar integrationElement = geometry.integrationElement(qp.position());
                Scalar p = 0.0;
                if (!isBox)
                    p = solution[eIdx][pressureIdx];
                else
                {
                    // do interpolation with ansatz functions
                    std::vector<Dune::FieldVector<Scalar, 1> > shapeValues;
                    const auto& localFiniteElement = feCache_.get(geometry.type());
                    localFiniteElement.localBasis().evaluateFunction(qp.position(), shapeValues);
                    for (unsigned int i = 0; i < shapeValues.size(); ++i)
                        p += shapeValues[i]*solution[gg.dofMapper().subIndex(element, i, dim)][pressureIdx];
                }
                lTwoNorm += (p - pe)*(p - pe)*qp.weight()*integrationElement;
            }
        }
        lTwoNorm = std::sqrt(lTwoNorm);

        // write the norm into a log file
        std::ofstream logFile;
        logFile.open(this->name() + ".log", std::ios::app);
        logFile << "[ConvergenceTest] L2-norm(pressure) = " << lTwoNorm  << " hMax = " << hMax_ << std::endl;
        logFile.close();
    }

private:

    Scalar exactPressure_(const GlobalPosition& globalPos) const
    {
        using std::sin;
        return sin(4.0*M_PI*globalPos[2]);
    }

    static constexpr Scalar eps_ = 1e-8;
    std::string name_;

    Scalar hMax_;

    typename Dune::PQkLocalFiniteElementCache<Scalar, Scalar, dim, 1> feCache_;
};

} //end namespace Dumux

#endif
