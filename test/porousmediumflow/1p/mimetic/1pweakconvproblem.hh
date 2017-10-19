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
#ifndef DUMUX_ONEPWEAKCONV_PROBLEM_HH
#define DUMUX_ONEPWEAKCONV_PROBLEM_HH


#include <dumux/porousmediumflow/implicit/problem.hh>
#if PROBLEM==1
#include <dumux/porousmediumflow/1p/mimetic/model.hh>
#include "resultevaluationmimetic.hh"
#else
#include <dumux/implicit/cellcentered/tpfa/properties.hh>
#include <dumux/implicit/cellcentered/mpfa/properties.hh>
#include <dumux/porousmediumflow/1p/implicit/model.hh>
#include "resultevaluationcc.hh"
#endif

#include <dumux/material/components/unit.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>

#include "1pweakconvspatialparams.hh"
#include <dune/geometry/quadraturerules.hh>

namespace Dumux
{
template <class TypeTag>
class OnePWeakConvProblem;

namespace Capabilities
{
    template<class TypeTag>
    struct isStationary<OnePWeakConvProblem<TypeTag>>
    { static const bool value = true; };
}

namespace Properties
{
#if PROBLEM==1
NEW_TYPE_TAG(OnePWeakConvProblem, INHERITS_FROM(OnePMimetic));
#else
NEW_TYPE_TAG(OnePWeakConvProblem, INHERITS_FROM(CCMpfaModel, OneP));

SET_PROP(OnePWeakConvProblem, MpfaMethod)
{
    static const MpfaMethods value = MpfaMethods::lMethod;
};
#endif

SET_PROP(OnePWeakConvProblem, Fluid)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef FluidSystems::LiquidPhase<Scalar, Dumux::Unit<Scalar> > type;
};

// Set the grid type
SET_TYPE_PROP(OnePWeakConvProblem, Grid, Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming>);

// Set the problem property
SET_TYPE_PROP(OnePWeakConvProblem, Problem, OnePWeakConvProblem<TypeTag> );

// Set the spatial parameters
SET_TYPE_PROP(OnePWeakConvProblem, SpatialParams, OnePWeakConvSpatialParams<TypeTag> );

// Enable gravity
SET_BOOL_PROP(OnePWeakConvProblem, ProblemEnableGravity, false);

SET_BOOL_PROP(OnePWeakConvProblem, EnableGlobalFVGeometryCache, true);

SET_BOOL_PROP(OnePWeakConvProblem, EnableGlobalFluxVariablesCache, true);
SET_BOOL_PROP(OnePWeakConvProblem, EnableGlobalVolumeVariablesCache, true);

SET_TYPE_PROP(OnePWeakConvProblem, LinearSolver, Dumux::SuperLUBackend<TypeTag> );
}


template <class TypeTag>
class OnePWeakConvProblem : public ImplicitPorousMediaProblem<TypeTag>
{
    typedef ImplicitPorousMediaProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };
    enum {
        // indices of the primary variables
        conti0EqIdx = Indices::conti0EqIdx,
        pressureIdx = Indices::pressureIdx,
    };

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::Intersection Intersection;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, SubControlVolume) SubControlVolume;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    typedef Dune::QuadratureRule<Scalar, dim> Quad;
    typedef typename Quad::iterator QuadIterator;
    typedef typename Dune::ReferenceElements<Scalar, dim> ReferenceElements;
    typedef typename Dune::ReferenceElement<Scalar, dim> ReferenceElement;

public:
    OnePWeakConvProblem(TimeManager &timeManager, const GridView &gridView)
    : ParentType(timeManager, gridView)
    {
        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag,
                                             std::string,
                                             Problem,
                                             Name);

        pi_ = 4.0*atan(1.0);

        testCase_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag,
                                                int,
                                                Problem,
                                                TestCase);

        if(testCase_ == 1)
        {
            alpha_ = 0.564519410563541;

            a_[0] = 1.0;
            a_[1] = -4.859104753080609;
            a_[2] = -4.859104753080598;
            a_[3] = -0.966437911685206;

            b_[0] = -12.041395892779573;
            b_[1] = -6.069941532218735;
            b_[2] = -6.069941532218734;
            b_[3] = -0.283720433459071;
        }
        else if(testCase_ == 2)
        {
            alpha_ = 0.614183857106858;
            a_[0] = 1.0;
            a_[1] = -0.427486203896442;
            a_[2] = -0.427486203896443;
            a_[3] = -0.760426111461169;
            b_[0] = -1.054628816800985;
            b_[1] = 0.214224702673554;
            b_[2] = 0.214224702673553;
            b_[3] = -0.649510087067298;
        }
        else if(testCase_ == 3)
        {
            alpha_ = 0.886633751528399;
            a_[0] = 1.0;
            a_[1] = -0.014387134650317;
            a_[2] = -0.014387134650317;
            a_[3] = 0.754437916582875;
            b_[0] = -0.370574113879094;
            b_[1] = 0.002197322794969;
            b_[2] = 0.002197322794969,
            b_[3] = -0.656381872494571;
        }

        beta_ = 7.0/16.0*pi_;
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

        auto K = this->spatialParams().permeability(element, scv,
                this->model().elementSolution(element, this->model().curSol()));

        int eIdx = this->elementMapper().index(element);

        //get the Gaussian quadrature rule for intervals
        const ReferenceElement& referenceElement = ReferenceElements::general(element.geometry().type());
        const Quad &quad = Dune::QuadratureRules<Scalar, dim>::rule(referenceElement.type(), 10);
        //std::cout << "Using quadrature rule of order: " << quad.order() << std::endl;
        const QuadIterator qend = quad.end();

        QuadIterator qp = quad.begin();
        for(; qp != qend; ++qp)
        {
            //std::cout << qp->position() << std::endl;
            GlobalPosition globalPos = element.geometry().global(qp->position());

            GlobalPosition gradRad = exactGradRad(globalPos);
            GlobalPosition dr = dr_(globalPos);
            GlobalPosition dtheta = dtheta_(globalPos);
            GlobalPosition ddr = ddr_(globalPos);
            GlobalPosition ddtheta = ddtheta_(globalPos);
            auto Hessian = HessianRad(globalPos);

            Scalar sxx = Hessian[0][0]*dr[0]*dr[0] + 2.0*Hessian[0][1]*dr[0]*dtheta[0] + Hessian[1][1]*dtheta[0]*dtheta[0]
                         + gradRad[0]*ddr[0] + gradRad[1]*ddtheta[0];
            Scalar syy = -sxx;

            Scalar integrationElement = element.geometry().integrationElement(qp->position());
            values[conti0EqIdx] += -(K[0][0]*sxx + K[1][1]*syy)*qp->weight()*integrationElement;
        }

        values[conti0EqIdx] /=element.geometry().volume();

        return values;
    }
    // \}
    /*!
     * \name Boundary conditions
     */
    // \{

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        values.setAllDirichlet();

        return values;
    }

    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(exact(globalPos));

        return values;
    }


    PrimaryVariables neumannAtPos(const GlobalPosition& globalPos) const
    {
        return PrimaryVariables(0);
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
        PrimaryVariables priVars(exact(globalPos));
        return priVars;
    }

    // \}

    void resultEvaluation()
    {
        result.evaluate(*this);

        std::ofstream file;
        std::string outname = name_;
        outname.append("_rates.txt");
        file.open(outname, std::ios::out | std::ios::app);
        if (file.fail())
            throw std::ios_base::failure(std::strerror(errno));

        file    << result.absL2Error  << "\t " << result.relativeL2Error  << "\t\t "<< result.absH1Error << "\t "
                << result.absH1ErrorApproxMin << "\t\t " << result.absH1ErrorDiffMin << "\t "
                << result.bilinearFormApprox <<"\t " << result.bilinearFormDiffApprox << "\t "
                << result.residualTermApprox <<"\t " << result.residualTermDiff << "\t "
                << result.residualTermError << "\t "
                << this->newtonController().newtonNumSteps() << "\t "
#if PROBLEM==1
                << -1 << "\t "
#else
                << this->model().jacobianAssembler().matrix().nonzeroes() << "\t "
#endif
                << result.hMax << "\t " << result.hMin << "\t "
                << result.uMax << "\t " << result.uMin << "\t "
                << this->gridView().size(0) << "\t " << std::endl;

        std::cout.setf(std::ios_base::scientific, std::ios_base::floatfield);
        std::cout.precision(2);

        std::cout
                << "\t absL2Error \t absH1ErrorDiffMin \t bilinearFormApprox \t bilinearFormDiffApprox \t hMax \t\t hMin \t\t numEle"
                << std::endl;
        std::cout << "\t " << result.absL2Error  << "\t "
                << result.absH1ErrorDiffMin << "\t\t "
                <<  result.bilinearFormApprox <<"\t\t " << result.bilinearFormDiffApprox << "\t\t\t "
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

    bool shouldWriteOutput() const
    {
        return
            this->timeManager().willBeFinished();
    }

    void postTimeStep()
    {
        ParentType::postTimeStep();

        if(this->shouldWriteOutput())
            this->resultEvaluation();
    }


    Scalar exact (const GlobalPosition& globalPos) const
    {
        Scalar theta = theta_(globalPos);
        Scalar r = r_(globalPos);

        int region = region_(globalPos);

        Scalar press = 10.0 + pow(r,alpha_)*(a_[region-1]*cos(alpha_*theta) + b_[region-1]*sin(alpha_*theta));

        return press;

    }

    Dune::FieldVector<Scalar,dim> exactGrad (const GlobalPosition& globalPos, const GlobalPosition& regionPos) const
    {
        Scalar theta = theta_(globalPos);
        Scalar r = r_(globalPos);

        GlobalPosition grad(0);

        GlobalPosition gradRad = exactGradRad(globalPos);
        GlobalPosition dr = dr_(globalPos);
        GlobalPosition dtheta = dtheta_(globalPos);

        grad[0] =  gradRad[0]*dr[0] + gradRad[1]*dtheta[0];
        grad[1] =  gradRad[0]*dr[1] + gradRad[1]*dtheta[1];

        return grad;
    }

    Dune::FieldVector<Scalar,dim> exactGradRad(const GlobalPosition& globalPos) const
    {
        Scalar theta = theta_(globalPos);
        Scalar r = r_(globalPos);

        int region = region_(globalPos);

        GlobalPosition grad(0);

        grad[0] =  alpha_*pow(r,(alpha_ - 1.0))*(a_[region-1]*cos(alpha_*theta) + b_[region-1]*sin(alpha_*theta));
        grad[1] =  pow(r,alpha_)*(alpha_*b_[region-1]*cos(alpha_*theta) - a_[region-1]*alpha_*sin(alpha_*theta));

        return grad;
    }

    Dune::FieldMatrix<Scalar,dim,dim> HessianRad(const GlobalPosition& globalPos) const
    {
        Scalar theta = theta_(globalPos);
        Scalar r = r_(globalPos);

        int region = region_(globalPos);

        Dune::FieldMatrix<Scalar,dim,dim> Hessian(0);

        Hessian[0][0] = (alpha_ - 1.0)*alpha_*pow(r,(alpha_ - 2.0))*(a_[region-1]*cos(alpha_*theta) + b_[region-1]*sin(alpha_*theta));
        Hessian[0][1] = alpha_*alpha_*pow(r,(alpha_ - 1))*(b_[region-1]*cos(alpha_*theta) - a_[region-1]*sin(alpha_*theta));
        Hessian[1][0] = Hessian[0][1];
        Hessian[1][1] = -alpha_*alpha_*pow(r,alpha_)*(a_[region-1]*cos(alpha_*theta) + b_[region-1]*sin(alpha_*theta));

        return Hessian;
    }

    Scalar r_(const GlobalPosition& globalPos) const
    {
        Scalar x = globalPos[0];
        Scalar y = globalPos[1];

        return std::sqrt(x*x + y*y);
    }

    GlobalPosition dr_(const GlobalPosition& globalPos) const
    {
        Scalar x = globalPos[0];
        Scalar y = globalPos[1];

        GlobalPosition dr(0);

        dr[0] = x/r_(globalPos);
        dr[1] = y/r_(globalPos);

        return dr;
    }

    GlobalPosition ddr_(const GlobalPosition& globalPos) const
    {
        Scalar x = globalPos[0];
        Scalar y = globalPos[1];

        GlobalPosition ddr(0);

        ddr[0] = 1.0/r_(globalPos) - (x*x)/std::pow(r_(globalPos),3);
        ddr[1] = 1.0/r_(globalPos) - (y*y)/std::pow(r_(globalPos),3);

        return ddr;
    }

    Scalar theta_(const GlobalPosition& globalPos) const
    {
        Scalar theta = 0.0;
        Scalar x = globalPos[0];
        Scalar y = globalPos[1];

        int region = region_(globalPos);

        if(region == 1 || region == 2) theta = std::atan2(y,x);
        else if(region == 3 || region == 4) theta = 2*pi_ + std::atan2(y,x);

        return theta;
    }

    GlobalPosition dtheta_(const GlobalPosition& globalPos) const
    {
        Scalar x = globalPos[0];
        Scalar y = globalPos[1];
        Scalar r = r_(globalPos);

        GlobalPosition dtheta(0);

        dtheta[0] = -y/(std::pow(r,2.0));
        dtheta[1] =  x/(std::pow(r,2.0));

        return dtheta;
    }

    GlobalPosition ddtheta_(const GlobalPosition& globalPos) const
    {
        Scalar x = globalPos[0];
        Scalar y = globalPos[1];
        Scalar r = r_(globalPos);

        GlobalPosition ddtheta(0);

        ddtheta[0] = 2.0*x*y/(std::pow(r,4.0));
        ddtheta[1] =  -2.0*x*y/(std::pow(r,4.0));

        return ddtheta;
    }


    int region_(const GlobalPosition& globalPos) const
    {
        int region = 1;

        Scalar x = globalPos[0];
        Scalar y = globalPos[1];

        if((y - tan(beta_)*x)<=0 && y>=0) region = 1;
        else if((y - tan(beta_)*x)>0 && y>=0) region = 2;
        else if((y - tan(beta_)*x)>0 && y<0) region = 3;
        else region = 4;

        return region;
    }

private:
    std::string name_;
    unsigned int testCase_;
    Dune::FieldVector<Scalar,4> a_;
    Dune::FieldVector<Scalar,4> b_;
    Scalar alpha_;
    Scalar beta_;
    Scalar pi_;
    Dumux::ResultEvaluation<TypeTag> result;
};
} //end namespace

#endif
