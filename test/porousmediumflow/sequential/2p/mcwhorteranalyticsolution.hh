// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *         See the file COPYING for full copying permissions.                *
 *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
#ifndef DUMUX_MCWHORTER_ANALYTIC_HH
#define DUMUX_MCWHORTER_ANALYTIC_HH

#include <dune/common/math.hh>
#include <dumux/porousmediumflow/2p/sequential/properties.hh>

/**
 * \file
 * \ingroup SequentialTwoPTests
 * \brief  Analytic solution of the McWhorter problem
 */

namespace Dumux
{
/**
 * \ingroup SequentialTwoPTests
 * \brief Analytic solution of
 * the McWhorter problem
 *
 * for naming of variables see "An Improved Semi-Analytical Solution for Verification
 * of Numerical Models of Two-Phase Flow in Porous Media"
 * (R. Fucik, J. Mikyska, M. Benes, T. H. Illangasekare; 2007)
 */

template<typename TypeTag>
class McWhorterAnalytic
{
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    using SpatialParams = GetPropType<TypeTag, Properties::SpatialParams>;

    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using FluidState = GetPropType<TypeTag, Properties::FluidState>;

    using SolutionTypes = GetProp<TypeTag, Properties::SolutionTypes>;
    using PrimaryVariables = typename SolutionTypes::PrimaryVariables;

    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    enum
    {
        dimworld = GridView::dimensionworld
    };
    enum
    {
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
        saturationIdx = Indices::saturationIdx
    };

    using BlockVector = Dune::BlockVector<Dune::FieldVector<Scalar, 1> >;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

private:

    // functions needed for analytical solution

    void initializeAnalytic()
    {
        int size = problem_.gridView().size(0);
        analyticSolution_.resize(size);
        analyticSolution_=0;
        errorGlobal_.resize(size);
        errorGlobal_ = 0;
        errorLocal_.resize(size);
        errorLocal_ = 0;

        return;
    }

    void calcSatError()
    {
        Scalar globalVolume = 0;
        Scalar errorNorm = 0.0;

        for (const auto& element : elements(problem_.gridView()))
        {
            // get entity
            int index = problem_.variables().index(element);

            Scalar sat = problem_.variables().cellData(index).saturation(wPhaseIdx);

            Scalar volume = element.geometry().volume();

            Scalar error = analyticSolution_[index] - sat;

            errorLocal_[index] = error;

            if (sat > swr_ + 1e-6)
            {
            globalVolume += volume;
            errorNorm += (volume*volume*error*error);
            }
        }

        if (globalVolume > 0.0 && errorNorm > 0.0)
        {
            using std::sqrt;
            errorNorm = sqrt(errorNorm)/globalVolume;
            errorGlobal_ = errorNorm;
        }
        else
        {
            errorGlobal_ = 0;
        }

        return;
    }

    void prepareAnalytic()
    {
        const auto& dummyElement = *problem_.gridView().template begin<0>();
        const auto& fluidMatrixInteraction = problem_->spatialParams().fluidMatrixInteractionAtPos(dummyElement.geometry().center());

        swr_ = fluidMatrixInteraction.effToAbsParams().swr();
        snr_ = fluidMatrixInteraction.effToAbsParams().snr();
        porosity_ = problem_.spatialParams().porosity(dummyElement);
        permeability_ = problem_.spatialParams().intrinsicPermeability(dummyElement)[0][0];
        PrimaryVariables initVec;
        problem_.initial(initVec, dummyElement);
        sInit_ = initVec[saturationIdx];
        Scalar s0 =(1 - snr_ - swr_);
//        std::cerr<<"\n swr, snr: "<< swr_ << snr_ << "\n";
//        std::cerr<<"\n sInit "<< sInit_<< "\n";

 //       h_= (s0 - sInit_)/intervalNum_;
        h_= (s0)/intervalNum_;
//         std::cout<<"h_= "<<h_<<std::endl;

        // define saturation range for analytic solution
        satVec_= swr_;
        for (int i=1; i<pointNum_; i++)
        {
            satVec_[i]=satVec_[i-1]+h_;
        }
        FluidState fluidState;
        fluidState.setTemperature(problem_.temperature(dummyElement));
        fluidState.setPressure(wPhaseIdx, problem_.referencePressure(dummyElement));
        fluidState.setPressure(nPhaseIdx, problem_.referencePressure(dummyElement));
        Scalar viscosityW = FluidSystem::viscosity(fluidState, wPhaseIdx);
        Scalar viscosityNW = FluidSystem::viscosity(fluidState, nPhaseIdx);

        // get fractional flow function vector
        for (int i=0; i<pointNum_; i++)
        {
            fractionalW_[i] = fluidMatrixInteraction.krw(satVec_[i])/viscosityW;
            fractionalW_[i] /= (fractionalW_[i] + fluidMatrixInteraction.krn(satVec_[i])/viscosityNW);
        }

        // get capillary pressure derivatives
        dpcdsw_=0;

        for (int i=0; i<pointNum_; i++)
            dpcdsw_[i] = fluidMatrixInteraction.dpc_dsw(satVec_[i]);
//         std::cout<<"dpcdsw = "<<dpcdsw_<<std::endl;

        // set initial fW
        if (sInit_ == 0)
        fInit_=0;
        else
        fInit_=fractionalW_[0];

        fractionalW_[0]=0;

        // normalize fW
        // with r_ = qt/q0
        // qt: total volume flux, q0: displacing phase flux at boundary
        // --> r_ = 1 for unidirectional displacement; r_ = 0 for impermeable boundary
        for (int i=0; i<pointNum_; i++)
        {
            fn_[i]= r_ * (fractionalW_[i] - fInit_)/ (1 - r_ * fInit_);
        }

        // std::cout<<"r_ = "<<r_<<std::endl;
        // std::cout<<"fn_ = "<<fn_<<std::endl;

        // diffusivity function
        for (int i=0; i<pointNum_; i++)
            d_[i] = fractionalW_[i]*(-dpcdsw_[i])*(fluidMatrixInteraction.krn(satVec_[i])/viscosityNW)*permeability_;

        // std::cout<<"fractionalW_ = "<<fractionalW_<<std::endl;
        // std::cout<<"permeability_ = "<<permeability_<<std::endl;
        // std::cout<<"d_ = "<<d_<<std::endl;


        // gk_: part of fractional flow function
        // initial guess for gk_
        for (int i=0; i<pointNum_; i++)
        {
            gk_[i] = d_[i]/(1-fn_[i]);
        }

        gk_[0] = 0;

//         std::cout<<"gk_ = "<<gk_<<std::endl;

        return;
    }

    void updateExSol()
    {
        Scalar time = problem_.timeManager().time() + problem_.timeManager().timeStepSize();

        // with displacing phase flux at boundary q0 = A * 1/sqrt(t)
        // Akm1, Ak: successive approximations of A
        double Ak = 0;
        double Akm1 = 0;
        double diff = 1e100;

        // partial numerical integrals a_, b_
        a_=0, b_=0;
        fp_=0;

        // approximation of integral I
        double I0 = 0;
        double Ii = 0;

        int k = 0;

        using std::pow;
        using Dune::power;
        using std::abs;

        while (diff> tolAnalytic_)
        {
            k++;
//             std::cout<<" k = "<<k<<std::endl;
            if (k> 50000)
            {
                std::cout<<"Analytic solution: Too many iterations!"<<std::endl;
                break;
            }

            Akm1=Ak;
            I0=0;
            for (int i=0; i<intervalNum_; i++)
            {
                a_[i] = 0.5 * h_ * sInit_ *(gk_[i] + gk_[i+1])+ power(h_, 2) / 6* ((3* i + 1) * gk_[i]
                        + (3 * i + 2) * gk_[i+1]);
                b_[i] = 0.5 * h_ * (gk_[i] + gk_[i+1]);
                I0 += (a_[i] - sInit_ * b_[i]);
            }
            // std::cout<<" I0 = "<<I0<<std::endl;

            gk_[0]=0;
            for (int i=1; i<pointNum_; i++)
            {
                Ii=0;
                for (int j = i; j<intervalNum_; j++)
                Ii += (a_[j] - satVec_[i] * b_[j]);
                //gk_[i] = d_[i] + gk_[i]*(fn_[i] + Ii/I0); // method A
                gk_[i] = (d_[i] + gk_[i]*fn_[i])/(1 - Ii/I0); // method B
            }

            // with f(sInit) = 0: relationship between A and sInit
            Ak = pow((0.5*porosity_/power((1 - fInit_), 2)*I0), 0.5);
            diff=abs(Ak - Akm1);
            // std::cout<<"diff = "<<diff<<std::endl;
        }

        // std::cout<<" b_ = "<<b_<<std::endl;
        // std::cout<<" Ak = "<<Ak<<std::endl;


        // fp_: first derivative of f
        for (int i = 0; i<pointNum_; i++)
        {
            for (int j = i; j<intervalNum_; j++)
            fp_[i] += b_[j]/I0;
        }

        // std::cout<<" fp_ = "<<fp_<<std::endl;
        // std::cout<<fInit_<<std::endl;

        for (int i = 0; i<pointNum_; i++)
        {
            xf_[i]= 2 * Ak * (1 - fInit_ * r_)/ porosity_ * fp_[i]* pow(time, 0.5);
        }

//         std::cout<<" xf_ = "<<xf_<<std::endl;

        // iterate over vertices and get analytical saturation solution
        for (const auto& element : elements(problem_.gridView()))
        {
            // get global coordinate of cell center
            const GlobalPosition& globalPos = element.geometry().center();

//            std::cout<<"globalPos = "<<globalPos[0]<<", x0 = "<<xf_[0]<<"\n";

            // find index of current vertex
            int index = problem_.variables().index(element);

            // find x_f next to global coordinate of the vertex
            int xnext = 0;
            for (int i=intervalNum_; i>0; i--)
            {
                if (globalPos[0] <= xf_[i])
                {
                    xnext = i;
                    break;
                }
            }

            // account for the area not yet reached by the front
            if (globalPos[0]> xf_[0])
            {
                analyticSolution_[index] = sInit_;
                continue;
            }

            if (globalPos[0] <= xf_[0])
            {
                analyticSolution_[index] = satVec_[xnext];
                continue;
            }

//             std::cout<<"Analytical = "<<satVec_[xnext]<<std::endl;
        }

        // call function to calculate the saturation error
        calcSatError();
        return;
    }

public:
    void calculateAnalyticSolution()
    {
        initializeAnalytic();

        updateExSol();
    }

    BlockVector AnalyticSolution() const
    {
        return analyticSolution_;
    }

    //Write saturation and pressure into file
     template<class MultiWriter>
     void addOutputVtkFields(MultiWriter &writer)
     {
         int size = problem_.gridView().size(0);
         BlockVector *analyticSolution = writer.allocateManagedBuffer (size);
         BlockVector *errorGlobal = writer.allocateManagedBuffer (size);
         BlockVector *errorLocal = writer.allocateManagedBuffer (size);

         *analyticSolution = analyticSolution_;
         *errorGlobal = errorGlobal_;
         *errorLocal = errorLocal_;

         writer.attachCellData(*analyticSolution, "saturation (exact solution)");
         writer.attachCellData(*errorGlobal, "global error");
         writer.attachCellData(*errorLocal, "local error");

         return;
     }


     void initialize()
     {
         initializeAnalytic();
         prepareAnalytic();
     }

    McWhorterAnalytic(Problem& problem, Scalar tolAnalytic = 1e-14) :
        problem_(problem), analyticSolution_(0), errorGlobal_(0), errorLocal_(0), tolAnalytic_(tolAnalytic)
    {}

private:
    Problem& problem_;

    BlockVector analyticSolution_;
    BlockVector errorGlobal_;
    BlockVector errorLocal_;

    Scalar tolAnalytic_;
    Scalar r_;

    Scalar swr_;
    Scalar snr_;
    Scalar porosity_;
    Scalar sInit_;
    Scalar permeability_;
    enum
    {   intervalNum_ = 1000, pointNum_ = intervalNum_+1};
    Dune::FieldVector<Scalar, pointNum_> satVec_;
    Dune::FieldVector<Scalar,pointNum_> fractionalW_;
    Dune::FieldVector<Scalar, pointNum_> dpcdsw_;
    Dune::FieldVector<Scalar, pointNum_> fn_;
    Dune::FieldVector<Scalar, pointNum_> d_;
    Dune::FieldVector<Scalar, pointNum_> gk_;
    Dune::FieldVector<Scalar, pointNum_> xf_;
    Scalar fInit_;
    Scalar h_;
    Dune::FieldVector<Scalar,intervalNum_> a_;
    Dune::FieldVector<Scalar,intervalNum_> b_;
    Dune::FieldVector<Scalar,pointNum_> fp_;

};
}
#endif
