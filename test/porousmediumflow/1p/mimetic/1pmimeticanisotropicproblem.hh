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
#ifndef DUMUX_1PMIMETICANISOTROPIC_PROBLEM_HH
#define DUMUX_1PMIMETICANISOTROPIC_PROBLEM_HH

//#include <dumux/implicit/cellcentered/tpfa/properties.hh>
//#include <dumux/implicit/cellcentered/mpfa/properties.hh>
//#include <dumux/porousmediumflow/1p/implicit/model.hh>
#include <dumux/porousmediumflow/1p/mimetic/model.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>
#include <dumux/material/components/unit.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>

#include "1pmimeticanisotropicspatialparams.hh"
#include "resultevaluationmimetic.hh"

namespace Dumux
{
template <class TypeTag>
class OnePMimeticAnisotropicProblem;

namespace Capabilities
{
    template<class TypeTag>
    struct isStationary<OnePMimeticAnisotropicProblem<TypeTag>>
    { static const bool value = true; };
}

namespace Properties
{
//NEW_TYPE_TAG(OnePMimeticAnisotropicProblem, INHERITS_FROM(CCMpfaModel, OneP, OnePMimeticTestSpatialParams));
NEW_TYPE_TAG(OnePMimeticAnisotropicProblem, INHERITS_FROM(OnePMimetic, OnePMimeticTestSpatialParams));

SET_PROP(OnePMimeticAnisotropicProblem, Fluid)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef FluidSystems::LiquidPhase<Scalar, Dumux::Unit<Scalar> > type;
};

// Set the grid type
SET_TYPE_PROP(OnePMimeticAnisotropicProblem, Grid, Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming>);

// Set the problem property
SET_TYPE_PROP(OnePMimeticAnisotropicProblem, Problem, Dumux::OnePMimeticAnisotropicProblem<TypeTag> );

// Set the spatial parameters
SET_TYPE_PROP(OnePMimeticAnisotropicProblem, SpatialParams, Dumux::OnePMimeticTestSpatialParams<TypeTag> );


SET_BOOL_PROP(OnePMimeticAnisotropicProblem, EnableGlobalFVGeometryCache, true);

SET_BOOL_PROP(OnePMimeticAnisotropicProblem, EnableGlobalFluxVariablesCache, true);
SET_BOOL_PROP(OnePMimeticAnisotropicProblem, EnableGlobalVolumeVariablesCache, true);

// Enable gravity
SET_BOOL_PROP(OnePMimeticAnisotropicProblem, ProblemEnableGravity, false);

SET_TYPE_PROP(OnePMimeticAnisotropicProblem, LinearSolver, UMFPackBackend<TypeTag> );

SET_BOOL_PROP(OnePMimeticAnisotropicProblem, VtkWriteFaceData, false);
}

template <class TypeTag>
class OnePMimeticAnisotropicProblem : public ImplicitPorousMediaProblem<TypeTag>
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
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

    enum { isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox) };

    typedef typename GridView::Intersection Intersection;

    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimWorldMatrix;

    //quadrature Rule
    typedef Dune::QuadratureRule<Scalar, dim> Quad;
    typedef typename Quad::iterator QuadIterator;
    typedef typename Dune::ReferenceElements<Scalar, dim> ReferenceElements;
    typedef typename Dune::ReferenceElement<Scalar, dim> ReferenceElement;

    typedef Dune::QuadratureRule<Scalar, dim-1> QuadFace;
    typedef typename QuadFace::iterator QuadIteratorFace;
    typedef typename Dune::ReferenceElements<Scalar, dim-1> ReferenceElementsFace;
    typedef typename Dune::ReferenceElement<Scalar, dim-1> ReferenceElementFace;

    typedef typename GridView::ctype CoordScalar;

    typedef Dune::PQkLocalFiniteElementCache<CoordScalar, Scalar, dim, 1> LocalFiniteElementCache;
    typedef typename LocalFiniteElementCache::FiniteElementType LocalFiniteElement;

public:
    OnePMimeticAnisotropicProblem(TimeManager &timeManager, const GridView &gridView)
    : ParentType(timeManager, gridView)
    {
        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag,
                                             std::string,
                                             Problem,
                                             Name);

        testCase_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag,
                                                int,
                                                Problem,
                                                TestCase);

        pi_ = 4.0*atan(1.0);
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
        Scalar source = 0.0;

        if(testCase_ == 1)
        {
            auto K = this->spatialParams().permeability(element, scv,
                    this->model().elementSolution(element, this->model().curSol()));

            //get the Gaussian quadrature rule for intervals
            const ReferenceElement& referenceElement = ReferenceElements::general(element.geometry().type());
            const Quad &quad = Dune::QuadratureRules<Scalar, dim>::rule(referenceElement.type(), 5);
            //std::cout << "Using quadrature rule of order: " << quad.order() << std::endl;
            const QuadIterator qend = quad.end();

            QuadIterator qp = quad.begin();
            for(; qp != qend; ++qp)
            {
                GlobalPosition globalPos = element.geometry().global(qp->position());

                Scalar integrationElement = element.geometry().integrationElement(qp->position());
                source += (pi_*pi_*((K[0][0]+K[1][1])*sin(pi_*globalPos[0])*sin(pi_*globalPos[1]) - (K[0][1]+K[1][0])*cos(pi_*globalPos[0])*cos(pi_*globalPos[1])))*qp->weight()*integrationElement;
            }

            source /=element.geometry().volume();
        }
        else if(testCase_ == 2)
        {
            //get the Gaussian quadrature rule for intervals
            const ReferenceElement& referenceElement = ReferenceElements::general(element.geometry().type());
            const Quad &quad = Dune::QuadratureRules<Scalar, dim>::rule(referenceElement.type(), 5);
            //std::cout << "Using quadrature rule of order: " << quad.order() << std::endl;
            const QuadIterator qend = quad.end();

            Scalar alpha = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag,
                    Scalar,
                    Problem,
                    Alpha);

            QuadIterator qp = quad.begin();
            for(; qp != qend; ++qp)
            {
                GlobalPosition globalPos = element.geometry().global(qp->position());
                Scalar x = globalPos[0];
                Scalar y = globalPos[1];
                Scalar integrationElement = element.geometry().integrationElement(qp->position());
                Scalar val = (pi_*(x*cos(pi_*x)*sin(pi_*y) + y*cos(pi_*y)*sin(pi_*x) -
                              alpha*x*cos(pi_*x)*sin(pi_*y) - alpha*y*cos(pi_*y)*sin(pi_*x) +
                              x*x*pi_*sin(pi_*x)*sin(pi_*y) + y*y*pi_*sin(pi_*x)*sin(pi_*y) +
                              alpha*x*x*pi_*sin(pi_*x)*sin(pi_*y) + alpha*y*y*pi_*sin(pi_*x)*sin(pi_*y) +
                              2.0*x*y*pi_*cos(pi_*x)*cos(pi_*y) - 2.0*alpha*x*y*pi_*cos(pi_*x)*cos(pi_*y)))/(x*x + y*y);

                source += val*qp->weight()*integrationElement;
            }
            source /=element.geometry().volume();
        }
        else if(testCase_ == 3)
        {
            const GlobalPosition& globalPos = element.geometry().center();

            if(globalPos[0] > 0.5 + 1.0e-4 && globalPos[1] > 0.5 + 1.0e-4)
                source = 10.0;
        }

        PrimaryVariables values(source);

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

//        if(onUpperBoundary_(globalPos) || onLowerBoundary_(globalPos))
//            values.setAllNeumann();
//        else
//            values.setAllDirichlet();

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
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * For this method, the \a priVars parameter stores the mass flux
     * in normal direction of each component. Negative values mean
     * influx.
     */
    PrimaryVariables neumann(const Element &element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolumeFace &scvf) const
    {
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        auto gradu = this->exactGrad(scvf.center(), insideScv.center());

        auto K = this->spatialParams().permeability(element, insideScv,
                this->model().elementSolution(element, this->model().curSol()));

        GlobalPosition Knormal;
        K.mv(scvf.unitOuterNormal(), Knormal);
        Knormal*=-1;
        PrimaryVariables values(Knormal*gradu);
        return values;
    }

    // \}

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
        PrimaryVariables priVars(0);
        if(testCase_ == 1)
            priVars[pressureIdx] = 1.0;//exact(element.geometry().center());
        else if(testCase_ == 3)
            priVars[pressureIdx] = 0.0;//exact(element.geometry().center());
        else
            priVars[pressureIdx] = 1.0;//exact(element.geometry().center());

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

            int bla = 0;

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

        Scalar press = 0.0;
        if(testCase_ != 3)
        {
            press = 1.0 + sin(pi_*x)*sin(pi_*y);
        }

        return press;

    }

    Dune::FieldVector<Scalar,dim> exactGrad (const GlobalPosition& evalPos, const GlobalPosition& regionPos) const
    {

        Scalar x = evalPos[0];
        Scalar y = evalPos[1];

        GlobalPosition grad(0);

        grad[0] =  pi_*cos(pi_*x)*sin(pi_*y);
        grad[1] =  pi_*sin(pi_*x)*cos(pi_*y);

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

    std::string name_;
    Scalar pi_;
    unsigned int testCase_;
    Dumux::ResultEvaluation<TypeTag> result;
    static constexpr Scalar eps_ = 3e-6;
};
} //end namespace

#endif
