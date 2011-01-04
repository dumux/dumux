// $Id$
/*****************************************************************************
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
#ifndef DUMUX_MCWHORTER_ANALYTIC_HH
#define DUMUX_MCWHORTER_ANALYTIC_HH

/**
 * @file
 * @brief  Analytic solution of
 * the McWhorter problem
 * @author Markus Wolff, Anneli Schöniger
 */

namespace Dumux
{
/**
 * @brief Analytic solution of
 * the McWhorter problem
 *
 * for naming of variables see "An Improved Semi-Analytical Solution for Verification
 * of Numerical Models of Two-Phase Flow in Porous Media"
 * (R. Fucik, J. Mikyska, M. Benes, T. H. Illangasekare; 2007)
 */

template<class TypeTag>
class McWhorterAnalytic
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SpatialParameters)) SpatialParameters;
    typedef typename SpatialParameters::MaterialLaw MaterialLaw;
    typedef typename MaterialLaw::Params MaterialLawParams;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidState)) FluidState;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices)) Indices;

    enum
    {
        dim = GridView::dimension, dimworld = GridView::dimensionworld
    };
    enum
    {
        wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx
    };

    typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > BlockVector;
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef Dune::FieldVector<Scalar, dimworld> GlobalPosition;

public:

    // functions needed for analytical solution

    void initializeAnalytic()
    {
        analyticSolution_.resize(size_);
        analyticSolution_=0;
        error_.resize(size_);
        error_=0;
        elementVolume_.resize(size_);
        elementVolume_=0;

        return;
    }

    void calcSatError(BlockVector &Approx)
    {
        error_=0;
        elementVolume_=0;

        ElementIterator eItEnd = problem_.gridView().template end<0> ();

        for (ElementIterator eIt = problem_.gridView().template end<0> (); eIt!= eItEnd; ++eIt)
        {
            // get entity
            const Element& element = *eIt;
            int index = problem_.variables().index(*eIt);
            elementVolume_[index]= element.geometry().volume();
            // std::cout<<"elementVolume_ = "<<elementVolume_[index]<<std::endl;
        }

        double globalVolume = elementVolume_.one_norm();
         //std::cout<<"globalVolume = "<<globalVolume<<std::endl; TODO: globalVolume always zero!

        for (int i=0; i<size_; i++)
        {
            error_[i]=analyticSolution_[i]-Approx[i];
             //std::cout<<"error_ = "<<error_[i]<<std::endl;
             //std::cout<<"analyticSolution_ = "<<analyticSolution_[i]<<std::endl;
             //std::cout<<"Approx = "<<Approx[i]<<std::endl;
        }

        // std::cout<<"error_ = "<<error_<<std::endl;

        double diffNorm = error_.two_norm();
        // std::cout<<"diffNorm = "<<diffNorm<<std::endl;

        for (int i=0; i<size_; i++)
        {
            error_[i] = (globalVolume) ? diffNorm * pow((elementVolume_[i]/globalVolume), 0.5) : 0;
        }

        return;
    }

    void prepareAnalytic()
    {
        const MaterialLawParams& materialLawParams(problem_.spatialParameters().materialLawParams(dummyGlobal_, dummyElement_));

        swr_ = materialLawParams.Swr();
        snr_ = materialLawParams.Snr();
        porosity_ = problem_.spatialParameters().porosity(dummyGlobal_, dummyElement_);
        permeability_ = problem_.spatialParameters().intrinsicPermeability(dummyGlobal_, dummyElement_)[0][0];
        sInit_ = problem_.initSat(dummyGlobal_, dummyElement_);
        Scalar s0 =(1 - snr_ - swr_);
        time_=0;

        h_= (s0 - sInit_)/intervalNum_;
//         std::cout<<"h_= "<<h_<<std::endl;

        // define saturation range for analytic solution
        satVec_= swr_;
        for (int i=1; i<pointNum_; i++)
        {
            satVec_[i]=satVec_[i-1]+h_;
        }
        FluidState fluidState;
        Scalar temp = problem_.temperature(dummyGlobal_, dummyElement_);
        Scalar press = problem_.referencePressure(dummyGlobal_, dummyElement_);
        Scalar viscosityW = FluidSystem::phaseViscosity(wPhaseIdx, temp, press, fluidState);
        Scalar viscosityNW = FluidSystem::phaseViscosity(nPhaseIdx, temp, press, fluidState);


        // get fractional flow function vector
        for (int i=0; i<pointNum_; i++)
        {
            fractionalW_[i] = MaterialLaw::krw(materialLawParams, satVec_[i])/viscosityW;
            fractionalW_[i] /= (fractionalW_[i] + MaterialLaw::krn(materialLawParams, satVec_[i])/viscosityNW);
        }

        // get capillary pressure derivatives
        dpcdsw_=0;

        for (int i=0; i<pointNum_; i++)
        {
            dpcdsw_[i] = MaterialLaw::dpC_dSw(materialLawParams, satVec_[i]);
                    }
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
        {d_[i] = fractionalW_[i]*(-dpcdsw_[i])*(MaterialLaw::krn(materialLawParams, satVec_[i])/viscosityNW)*permeability_;
         }

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
                a_[i] = 0.5 * h_ * sInit_ *(gk_[i] + gk_[i+1])+ pow(h_, 2) / 6* ((3* i + 1) * gk_[i]
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
            Ak = pow((0.5*porosity_/pow((1 - fInit_), 2)*I0), 0.5);
            diff=fabs(Ak - Akm1);
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
            xf_[i]= 2 * Ak * (1 - fInit_ * r_)/ porosity_ * fp_[i]* pow(time_, 0.5);
        }

//         std::cout<<" xf_ = "<<xf_<<std::endl;

        // iterate over vertices and get analytical saturation solution
        ElementIterator eItEnd = problem_.gridView().template end<0> ();
        for (ElementIterator eIt = problem_.gridView().template begin<0> (); eIt!= eItEnd; ++eIt)
        {
            // get global coordinate of cell center
            const GlobalPosition& globalPos = eIt->geometry().center();

//            std::cout<<"globalPos = "<<globalPos[0]<<", x0 = "<<xf_[0]<<"\n";

            // find index of current vertex
            int index = problem_.variables().index(*eIt);

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
        calcSatError(problem_.variables().saturation());
        return;
    }

    void calculateAnalyticSolution()
    {
        time_ += problem_.timeManager().timeStepSize();
        updateExSol();
    }

    //Write saturation and pressure into file
     template<class MultiWriter>
     void addOutputVtkFields(MultiWriter &writer)
     {
         BlockVector *analyticSolution = writer.template createField<Scalar, 1> (size_);
         BlockVector *error = writer.template createField<Scalar, 1> (size_);

         *analyticSolution = analyticSolution_;
         *error = error_;

         writer.addCellData(analyticSolution, "saturation (exact solution)");
         writer.addCellData(error, "error_");

         return;
     }

    McWhorterAnalytic(Problem& problem, Scalar tolAnalytic = 1e-14) :
        problem_(problem), analyticSolution_(0), error_(0), elementVolume_(0), size_(problem.gridView().size(0)), dummyElement_(
                *(problem_.gridView().template begin<0> ())), dummyGlobal_(GlobalPosition(1)), tolAnalytic_(tolAnalytic)
    {
        initializeAnalytic();
        prepareAnalytic();
    }

protected:


private:
    Problem& problem_;

    BlockVector analyticSolution_;
    BlockVector error_;
    BlockVector elementVolume_;

    int size_;

    const Element& dummyElement_;
    const GlobalPosition& dummyGlobal_;

    Scalar tolAnalytic_;
    Scalar r_;

    Scalar swr_;
    Scalar snr_;
    Scalar porosity_;
    Scalar sInit_;
    Scalar time_;
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
