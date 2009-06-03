// $Id: impes.hh 972 2009-01-12 10:15:57Z lauser $

#ifndef DUNE_IMPES_BUCKLEYLEVERETT_ANALYTICAL_HH
#define DUNE_IMPES_BUCKLEYLEVERETT_ANALYTICAL_HH

#include "dumux/fractionalflow/impes/impes.hh"

/**
 * @file
 * @brief  IMPES scheme
 * @author Markus Wolff, Anneli Schöniger
 */

namespace Dune
{
/**
 * \ingroup fracflow
 * @brief IMplicit Pressure Explicit Saturation (IMPES) scheme for the solution of
 * the McWhorter problem
 *
 * for naming of variables see "An Improved Semi-Analytical Solution for Verification
 * of Numerical Models of Two-Phase Flow in Porous Media"
 * (R. Fucik, J. Mikyska, M. Benes, T. H. Illangasekare; 2007)
 */

template<class GridView, class Diffusion, class Transport, class VC> class IMPESMcWAnalytic: public IMPES<
        GridView, Diffusion, Transport, VC>
{
    typedef Dune::FractionalFlow<GridView, Diffusion, Transport, VC>
            FractionalFlow;
typedef    typename FractionalFlow::RepresentationType PressType;
    typedef typename FractionalFlow::Scalar Scalar;
    enum
    {   dim=GridView::dimension,dimworld = GridView::dimensionworld};
    typedef Dune::BlockVector<FieldVector<Scalar, 1> > BlockVector;
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::Grid Grid;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::template Codim<0>::EntityPointer ElementPointer;
    typedef Dune::FieldVector<Scalar,dim> LocalPosition;
    typedef Dune::FieldVector<Scalar,dimworld> GlobalPosition;

public:
    typedef typename Transport::RepresentationType RepresentationType;

    // functions needed for analytical solution

    void initializeAnalytic()
    {
        analyticSolution.resize(size_);
        analyticSolution=0;
        error.resize(size_);
        error=0;
        elementVolume.resize(size_);
        elementVolume=0;

        return;
    }

    void calcSatError(BlockVector &Approx)
    {
        error=0;
        elementVolume=0;

        ElementIterator eendit = Transport::gridView.template end<0>();

        for (ElementIterator it = Transport::gridView.template begin<0>(); it!= eendit; ++it)
        {
            // get entity
            const Element& element = *it;
            int index = this->transProblem.variables().indexTransport(*it);
            elementVolume[index]= element.geometry().volume();
            // std::cout<<"elementVolume = "<<elementVolume[index]<<std::endl;
        }

        double globalVolume = elementVolume.one_norm();
        // std::cout<<"globalVolume = "<<globalVolume<<std::endl;

        for (int i=0; i<size_; i++)
        {
            error[i]=analyticSolution[i]-Approx[i];
            // std::cout<<"error = "<<error[i]<<std::endl;
            // std::cout<<"analyticSolution = "<<analyticSolution[i]<<std::endl;
            // std::cout<<"Approx = "<<Approx[i]<<std::endl;
        }

        // std::cout<<"error = "<<error<<std::endl;

        double diffNorm = error.two_norm();
        // std::cout<<"diffNorm = "<<diffNorm<<std::endl;

        for (int i=0; i<size_; i++)
        {
            error[i] = diffNorm * pow((elementVolume[i]/globalVolume), 0.5);
        }

        return;
    }

    void prepareAnalytic()
    {
        swr_ = this->transProblem.soil().Sr_w(dummyGlobal_, dummyElement_, dummyLocal_);
        snr_ = this->transProblem.soil().Sr_n(dummyGlobal_, dummyElement_, dummyLocal_);
        porosity_ = this->transProblem.soil().porosity(dummyGlobal_, dummyElement_, dummyLocal_);
        permeability_ = this->transProblem.soil().K(dummyGlobal_, dummyElement_, dummyLocal_)[0][0];
        sInit_=this->transProblem.initSat(dummyGlobal_, dummyElement_, dummyLocal_);

        time_=0;

        h_= (1-sInit_)/intervalNum_;
        // std::cout<<"h_= "<<h_<<std::endl;

        // define saturation range for analytic solution
        satVec_= 0;
        for (int i=1; i<pointNum_; i++)
        {
            satVec_[i]=satVec_[i-1]+h_;
        }

        // get fractional flow function vector
        for (int i=0; i<pointNum_; i++)
        {
            fractionalW_[i] = this->transProblem.materialLaw().fractionalW(satVec_[i], dummyGlobal_, dummyElement_, dummyLocal_);
        }

        // get capillary pressure derivatives
        dpcdsw_=0;

        for (int i=0; i<pointNum_; i++)
        {
            dpcdsw_[i] = this->transProblem.materialLaw().dPdS(satVec_[i], dummyGlobal_, dummyElement_, dummyLocal_);
        }
        // std::cout<<"dpcdsw = "<<dpcdsw_<<std::endl;

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
        {
            d_[i] = fractionalW_[i]*this->transProblem.materialLaw().mobN(1-satVec_[i], dummyGlobal_, dummyElement_, dummyLocal_)*(-dpcdsw_[i])*permeability_;
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

        // std::cout<<"gk_ = "<<gk_<<std::endl;

        return;
    }

    void setTime(double &dt)
    {
        time_+=dt;

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
            // std::cout<<" k = "<<k<<std::endl;
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

        for (int i = 0; i<pointNum_; i++)
        {
            xf_[i]= 2 * Ak * (1 - fInit_ * r_)/ porosity_ * fp_[i]* pow(time_, 0.5);
        }

        // std::cout<<" xf_ = "<<xf_[pointNum_]<<std::endl;

        // iterate over vertices and get analytical saturation solution
        ElementIterator eendit = Transport::gridView.template end<0>();
        for (ElementIterator it = Transport::gridView.template begin<0>(); it!= eendit; ++it)
        {
            // find index of current vertex
            int index = this->transProblem.variables().indexTransport(*it);

            // get global coordinate of vertex
            Dune::FieldVector<Scalar,dimworld> global = it->geometry().corner(0);

            // find x_f next to global coordinate of the vertex
            int xnext = 0;
            for (int i=intervalNum_; i>0; i--)
            {
                if (global[0]<xf_[i])
                {
                    xnext = i;
                    break;
                }
            }

            // account for the area not yet reached by the front
            if (global[0] >= xf_[0])
            {
                analyticSolution[index] = sInit_;
                continue;
            }

            if (global[0] < xf_[0])
            {
                analyticSolution[index] = satVec_[xnext];
                continue;
            }

            // std::cout<<"Analytical = "<<satVec_[xnext]<<std::endl;
        }

        // call function to calculate the saturation error
        calcSatError(this->transProblem.variables().saturation());

        return;
    }

    void postProcessUpdate(double t, double& dt)
    {
        setTime(dt);
        updateExSol();
    }

    virtual void vtkout(const char* name, int k) const
    {
        // vtkout for numerical and analytical solution
        VTKWriter<GridView> vtkwriter(Transport::gridView);
        char fname[128];
        sprintf(fname, "%s-%05d", name, k);
        vtkwriter.addCellData(this->transProblem.variables().saturation(), "saturation");
        vtkwriter.addCellData(this->transProblem.variables().pressure(), "total pressure p~");
        vtkwriter.addCellData(analyticSolution, "saturation (exact solution)");
        vtkwriter.addCellData(error, "error");
        vtkwriter.write(fname, VTKOptions::ascii);

        return;
    }

    //! Construct an IMPES object.
    IMPESMcWAnalytic(Diffusion& diffusion, Transport& transport, int flag = 0, int nIt = 2,
            Scalar maxDef = 1e-5, Scalar om = 1):
    IMPES<GridView, Diffusion, Transport, VC> (diffusion, transport, flag, nIt, maxDef, om),
    analyticSolution(0), error(0), elementVolume(0), size_(this->transProblem.variables().gridSizeTransport()),
    dummyLocal_(LocalPosition(1)), dummyGlobal_(GlobalPosition(1)), dummyElement_(*(Transport::gridView.template begin<0>())),
    r_(0)
    {
        initializeAnalytic();
        prepareAnalytic();
    }

protected:
    BlockVector analyticSolution;
    BlockVector error;
    BlockVector elementVolume;

private:
    int size_;
    const LocalPosition& dummyLocal_;
    const GlobalPosition& dummyGlobal_;
    const Element& dummyElement_;
    Scalar swr_;
    Scalar snr_;
    Scalar porosity_;
    Scalar sInit_;
    Scalar time_;
    Scalar permeability_;
    Scalar tolAnalytic_;
    enum
    {   intervalNum_ = 1000, pointNum_ = intervalNum_+1};
    FieldVector<Scalar, pointNum_> satVec_;
    FieldVector<Scalar,pointNum_> fractionalW_;
    FieldVector<Scalar, pointNum_> dpcdsw_;
    FieldVector<Scalar, pointNum_> fn_;
    FieldVector<Scalar, pointNum_> d_;
    FieldVector<Scalar, pointNum_> gk_;
    FieldVector<Scalar, pointNum_> xf_;
    Scalar fInit_;
    Scalar r_;
    Scalar h_;
    FieldVector<Scalar,intervalNum_> a_;
    FieldVector<Scalar,intervalNum_> b_;
    FieldVector<Scalar,pointNum_> fp_;
};
}
#endif
