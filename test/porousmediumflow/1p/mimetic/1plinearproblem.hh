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
 * \brief A test problem for the one-phase model:
 * water is flowing from bottom to top through and around a low permeable lens.
 */
#ifndef DUMUX_1PLINEAR_MIMETIC_PROBLEM_HH
#define DUMUX_1PLINEAR_MIMETIC_PROBLEM_HH

#if PROBLEM==1
#include <dumux/porousmediumflow/1p/mimetic/model.hh>
#include "resultevaluationmimetic.hh"
#include <dumux/discretization/staggered/mimetic/mimeticcpgeometryhelper.hh>
#else
#include <dumux/implicit/cellcentered/tpfa/properties.hh>
#include <dumux/implicit/cellcentered/mpfa/properties.hh>
#include <dumux/porousmediumflow/1p/implicit/model.hh>
#include "resultevaluationcc.hh"
#endif
#include <dumux/porousmediumflow/implicit/problem.hh>
#include <dumux/material/components/unit.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>

#include <dumux/io/cpgridcreator.hh>
#include <dumux/common/intersectionmapper.hh>

#include "1plinearspatialparams.hh"

namespace Dumux
{
template <class TypeTag>
class OnePLinearProblem;

namespace Capabilities
{
    template<class TypeTag>
    struct isStationary<OnePLinearProblem<TypeTag>>
    { static const bool value = true; };
}

namespace Properties
{
//NEW_TYPE_TAG(OnePLinearProblem, INHERITS_FROM(CCTpfaModel, OneP));
#if PROBLEM==1
NEW_TYPE_TAG(OnePLinearProblem, INHERITS_FROM(OnePMimetic));
#else
NEW_TYPE_TAG(OnePLinearProblem, INHERITS_FROM(CCMpfaModel, OneP));
#endif

SET_PROP(OnePLinearProblem, Fluid)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef FluidSystems::LiquidPhase<Scalar, Dumux::Unit<Scalar> > type;
};

// Set the grid type
SET_TYPE_PROP(OnePLinearProblem, Grid, Dune::CpGrid);
//SET_TYPE_PROP(OnePLinearProblem, Grid, Dune::PolyhedralGrid< 3, 3 >);
//SET_TYPE_PROP(OnePLinearProblem, Grid, Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming>);

// Set the grid creator
SET_TYPE_PROP(OnePLinearProblem, GridCreator, Dumux::CpGridCreator<TypeTag>);

// Set the problem property
SET_TYPE_PROP(OnePLinearProblem, Problem, Dumux::OnePLinearProblem<TypeTag> );

// Set the spatial parameters
SET_TYPE_PROP(OnePLinearProblem, SpatialParams, Dumux::OnePLinearSpatialParams<TypeTag> );


SET_BOOL_PROP(OnePLinearProblem, EnableGlobalFVGeometryCache, true);

SET_BOOL_PROP(OnePLinearProblem, EnableGlobalFluxVariablesCache, true);
SET_BOOL_PROP(OnePLinearProblem, EnableGlobalVolumeVariablesCache, true);

// Enable gravity
SET_BOOL_PROP(OnePLinearProblem, ProblemEnableGravity, false);

SET_TYPE_PROP(OnePLinearProblem, LinearSolver, SuperLUBackend<TypeTag> );

#if PROBLEM==1
SET_BOOL_PROP(OnePLinearProblem, VtkWriteFaceData, false);

// The geometry helper required for the stencils, etc.
SET_PROP(OnePLinearProblem, StaggeredGeometryHelper)
{
private:
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
public:
    using type = MimeticCPGeometryHelper<GridView>;
};

SET_TYPE_PROP(OnePLinearProblem, IntersectionMapper, Dumux::NonConformingGridIntersectionMapper<TypeTag>);
#endif

}

/*!
 * \ingroup OnePModel
 * \ingroup ImplicitTestProblems
 * \brief  Test problem for the one-phase model:
 * water is flowing from bottom to top through and around a low permeable lens.
 *
 * The domain is box shaped. All sides are closed (Neumann 0 boundary)
 * except the top and bottom boundaries (Dirichlet), where water is
 * flowing from bottom to top.
 *
 * In the middle of the domain, a lens with low permeability (\f$K=10e-12\f$)
 * compared to the surrounding material (\f$ K=10e-10\f$) is defined.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_box1p -parameterFile test_box1p.input</tt> or
 * <tt>./test_cc1p -parameterFile test_cc1p.input</tt>
 *
 * The same parameter file can be also used for 3d simulation but you need to change line
 * <tt>typedef Dune::YaspGrid<2> type;</tt> to
 * <tt>typedef Dune::YaspGrid<3> type;</tt> in the problem file
 * and use <tt>test_1p_3d.dgf</tt> in the parameter file.
 */


template <class TypeTag>
class OnePLinearProblem : public ImplicitPorousMediaProblem<TypeTag>
{
    using ParentType = ImplicitPorousMediaProblem<TypeTag>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);

    // Grid and world dimension
    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;

    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    enum {
     // indices of the primary variables
     conti0EqIdx = Indices::conti0EqIdx,
     pressureIdx = Indices::pressureIdx,
     //facePressureIdx = Indices::facePressureIdx
    };

    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using TimeManager = typename GET_PROP_TYPE(TypeTag, TimeManager);
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };

    typedef typename GridView::Intersection Intersection;

    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimWorldMatrix;

    typedef std::array<Scalar, 4> CoeffArray;

public:
    OnePLinearProblem(TimeManager &timeManager, const GridView &gridView)
    : ParentType(timeManager, gridView)
    {
        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag,
                                             std::string,
                                             Problem,
                                             Name);

        testCase_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, Problem, TestCase);

        a_ = GET_RUNTIME_PARAM_FROM_GROUP_CSTRING(TypeTag, CoeffArray, "Problem", a);
        b_ = GET_RUNTIME_PARAM_FROM_GROUP_CSTRING(TypeTag, CoeffArray, "Problem", b);

        pLow_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag,
                        Scalar,
                        Problem,
                        Plow);

        pHigh_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag,
                        Scalar,
                        Problem,
                        Phigh);
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
    std::string name() const
    {
        return name_;
    }

    /*!
     * \brief Return the temperature within the domain in [K].
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return 273.15 + 10; } // 10C

    PrimaryVariables source(const Element &element,
                const FVElementGeometry& fvGeometry,
                const ElementVolumeVariables& elemVolVars,
                const SubControlVolume &scv) const
    {
        PrimaryVariables values(0.0);
        return values;
    }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary control volume.
     *
     * \param values The boundary types for the conservation equations
     * \param globalPos The position of the center of the finite volume
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        values.setAllDirichlet();

        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     *
     * \param values The dirichlet values for the primary variables
     * \param globalPos The center of the finite volume which ought to be set.
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(exact(globalPos));

        return values;
    }


    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * For this method, the \a priVars parameter stores primary
     * variables.
     */
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables priVars(exact(globalPos));

        return priVars;
    }


    void resultEvaluation()
    {
        result.evaluate(*this);

        std::ofstream file;
        std::string outname = name_;
        outname.append("_rates.txt");
        file.open(outname, std::ios::out | std::ios::app);
        if (file.fail())
            throw std::ios_base::failure(std::strerror(errno));

        file    << result.absL2Error  << "\t\t " << result.relativeL2Error  << "\t\t "
                << result.absH1ErrorApproxMin << "\t\t " << result.absH1ErrorDiffMin << "\t\t "
                << this->newtonController().newtonNumSteps() << "\t\t "
                << result.hMax << "\t\t " << result.hMin << "\t\t "
                << result.uMax << "\t\t " << result.uMin << "\t\ "
                << this->gridView().size(0) << "\t " << std::endl;

        std::cout.setf(std::ios_base::scientific, std::ios_base::floatfield);
        std::cout.precision(2);

        std::cout
                << "\t absL2Error \t pRelErrorL2  \t absH1ErrorApproxMin \t absH1ErrorDiffMin \t hMax \t\t hMin \t\t numEle"
                << std::endl;
        std::cout << "\t " << result.absL2Error  << "\t " << result.relativeL2Error  << "\t "
                << result.absH1ErrorApproxMin << "\t\t " << result.absH1ErrorDiffMin << "\t\t "
                << result.hMax << "\t " << result.hMin << "\t " << this->gridView().size(0) << "\t " << std::endl;
    }

    /*!
     * \brief Adds additional VTK output data to the VTKWriter. Function is called by the output module on every write.
     */
    template<class VtkOutputModule>
    void addVtkOutputFields(VtkOutputModule& outputModule) const
    {
        auto& exactPressure = outputModule.createScalarField("exactPressure", 0);
        auto& volumes = outputModule.createScalarField("volumes", 0);

        Scalar pressureL2 = 0.0;

        for (const auto& element : elements(this->gridView()))
        {
            auto fvGeometry = localView(this->model().globalFvGeometry());
            fvGeometry.bindElement(element);

            for (auto&& scv : scvs(fvGeometry))
            {
                auto ccDofIdx = scv.dofIndex();
                auto ccDofPosition = scv.dofPosition();

                auto elemVolVars = localView(this->model().curGlobalVolVars());
                elemVolVars.bind(element, fvGeometry, this->model().curSol());

                Scalar exactP = exact(ccDofPosition);
                exactPressure[ccDofIdx] = exactP;
                volumes[ccDofIdx] = scv.volume();

                pressureL2 += (exactP-elemVolVars[scv].pressure())*(exactP-elemVolVars[scv].pressure())*scv.volume();
            }
        }

        std::cout << "L2 error of pressure: " << std::sqrt(pressureL2) << std::endl;
    }

    Scalar exact (const GlobalPosition& globalPos) const
    {
        Scalar x = globalPos[0];
        Scalar y = globalPos[1];
        Scalar z = 0.0;

        if(dimWorld == 3)
            z = globalPos[2];

        Scalar p_ex = 0.0;

        if(testCase_ == 1)
        {
            if(x < 0.6 +1.0e-6)
            {
                p_ex = a_[0]*x + a_[1]*y + a_[2]*z + a_[3];
            }
            else
            {
                p_ex = b_[0]*x + b_[1]*y + b_[2]*z + b_[3];
            }
        }
        else
        {
            Scalar x_max = this->bBoxMax()[0];
            Scalar x_min = this->bBoxMin()[0];
            Scalar y_max = this->bBoxMax()[1];
            Scalar y_min = this->bBoxMin()[1];

            Scalar z_max = 0.0;
            Scalar z_min = 0.0;

            if(dimWorld == 3)
            {
                z_max = this->bBoxMax()[2];
                z_min = this->bBoxMin()[2];
            }

            if(dimWorld == 3)
            {
            p_ex = 1.0/dimWorld*((x - x_min)/(x_max - x_min) + (y - y_min)/(y_max - y_min) + (z - z_min)/(z_max - z_min))*pLow_
                   - 1.0/dimWorld*((x - x_max)/(x_max - x_min) + (y - y_max)/(y_max - y_min) + (z - z_max)/(z_max - z_min))*pHigh_;
            }
            else
            {
            p_ex = 1.0/dimWorld*((x - x_min)/(x_max - x_min) + (y - y_min)/(y_max - y_min))*pLow_
                   - 1.0/dimWorld*((x - x_max)/(x_max - x_min) + (y - y_max)/(y_max - y_min))*pHigh_;
            }
        }

        return p_ex;
    }

    Dune::FieldVector<Scalar,dim> exactGrad(const GlobalPosition& evalPos, const GlobalPosition& regionPos) const
    {

        Scalar x = regionPos[0];

        Dune::FieldVector<Scalar,dim> grad(0);
        if(testCase_ == 1)
        {
            if(x < 0.6 +1.0e-6)
            {
                grad[0] = a_[0];
                grad[1] = a_[1];
                grad[2] = a_[2];
            }
            else
            {
                grad[0] = b_[0];
                grad[1] = b_[1];
                grad[2] = b_[2];
            }
        }
        else
        {
            Scalar x_max = this->bBoxMax()[0];
            Scalar x_min = this->bBoxMin()[0];
            Scalar y_max = this->bBoxMax()[1];
            Scalar y_min = this->bBoxMin()[1];

            Scalar z_max = 0.0;
            Scalar z_min = 0.0;

            if(dimWorld == 3)
            {
                z_max = this->bBoxMax()[2];
                z_min = this->bBoxMin()[2];
            }

            grad[0] = 1.0/dimWorld*(pLow_-pHigh_)/(x_max - x_min);
            grad[1] = 1.0/dimWorld*(pLow_-pHigh_)/(y_max - y_min);
            if(dimWorld == 3)
            {
                grad[2] = 1.0/dimWorld*(pLow_-pHigh_)/(z_max - z_min);
            }

        }

        return grad;
    }

    bool shouldWriteOutput() const
    {
        return
            this->timeManager().willBeFinished();
    }

    void postTimeStep()
    {
        ParentType::postTimeStep();

        if(this->timeManager().willBeFinished())
            this->resultEvaluation();
    }

    // \}

private:
    std::string name_;
    int testCase_;
    Dumux::ResultEvaluation<TypeTag> result;
    CoeffArray a_;
    CoeffArray b_;
    Scalar pLow_;
    Scalar pHigh_;
};
} //end namespace

#endif
