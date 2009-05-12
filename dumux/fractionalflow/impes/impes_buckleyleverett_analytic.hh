// $Id: impes.hh 972 2009-01-12 10:15:57Z lauser $

#ifndef DUNE_IMPES_BUCKLEYLEVERETT_ANALYTICAL_HH
#define DUNE_IMPES_BUCKLEYLEVERETT_ANALYTICAL_HH

#include "dumux/fractionalflow/impes/impes.hh"

/**
 * @file
 * @brief  IMPES scheme
 * @author Bernd Flemisch, last changed by Markus Wolff
 */

namespace Dune
{
/**
 * \ingroup fracflow
 * @brief IMplicit Pressure Explicit Saturation (IMPES) scheme for the solution of
 * the Buckley-Leverett problem
 */

template<class Grid, class Diffusion, class Transport, class VC> class IMPESBLAnalytic: public IMPES<
        Grid, Diffusion, Transport, VC>
{
    typedef Dune::FractionalFlow<Grid, Diffusion, Transport, VC> FractionalFlow;
typedef    typename FractionalFlow::RepresentationType PressType;
    typedef typename FractionalFlow::Scalar Scalar;
    enum
    {   dim=Grid::dimension,dimworld = Grid::dimensionworld};
    typedef Dune::BlockVector<FieldVector<Scalar, 1> > BlockVector;
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename Grid::LevelGridView GridView;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename Grid::template Codim<0>::EntityPointer ElementPointer;

    template<int dim> struct ElementLayout
    {
        bool contains(Dune::GeometryType gt)
        {
            return gt.dim() == dim;
        }
    };
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridView,ElementLayout> ElementMapper;
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
        const GridView& gridView(grid.levelView(0));
        ElementIterator eendit = gridView.template end<0>();

        for (ElementIterator it = gridView.template begin<0>(); it!= eendit; ++it)
        {
            // get entity
            const Element& element = *it;
            int index = mapper.map(*it);
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

        //        std::cout<<"error = "<<error<<std::endl;

        double diffNorm = error.two_norm();
        std::cout<<"diffNorm = "<<diffNorm<<std::endl;

        for (int i=0; i<size_; i++)
        {
            error[i] = diffNorm * pow((elementVolume[i]/globalVolume), 0.5);
        }

        return;
    }

    void prepareAnalytic()
    {
        Swr_ = this->transProblem.soil().Sr_w(dummyGlobal_, dummyElement_, dummyLocal_);
        Snr_ = this->transProblem.soil().Sr_n(dummyGlobal_, dummyElement_, dummyLocal_);

        time_=0;

        SatVec_=Swr_;
        for (int i=1; i<pointNum_; i++)
        {
            SatVec_[i]=SatVec_[i-1]+(1-Snr_-Swr_)/intervalNum_;
        }

        //        std::cout<<"SatVec_ = "<<SatVec_<<std::endl;
        for (int i=0; i<pointNum_; i++)
        {
            fractionalW_[i] = this->transProblem.materialLaw().fractionalW(SatVec_[i],dummyGlobal_, dummyElement_, dummyLocal_);
        }

        //        std::cout<<"fractionalW_ = "<<fractionalW_<<std::endl;

        dfwdsw_=0;

        for (int i=1; i<intervalNum_; i++)
        {
            dfwdsw_[i]=(fractionalW_[i+1]-fractionalW_[i-1])/(SatVec_[i+1]-SatVec_[i-1]);
        }

        //        std::cout<<"dfwdsw_ = "<<dfwdsw_<<std::endl;
        for (int i=0; i<pointNum_; i++)
        {
            if (dfwdsw_[i]>dfwdsw_[i+1])
            {
                dfwdswmax_ = i;
                break;
            }
        }

        return;
    }

    void setTime(double &dt)
    {
        time_+=dt*CFL_;
        return;
    }

    void updateExSol()
    {
        //position of the fluid front
        xf_=0;
        for (int i=0; i<pointNum_; i++)
        {
            xf_[i]=vtot_*time_/this->transProblem.soil().porosity(dummyGlobal_, dummyElement_, dummyLocal_)*dfwdsw_[i];
        }

        //std::cout<<"xf_ = "<<xf_<<std::endl;
        int xhelp=pointNum_/3;
        int xhelpold = 0;
        int xhelpoldold = 0;
        int xfmax = 0;

        //position of maximum xf_
        for (int i=0; i<pointNum_; i++)
        {
            if (xf_[i]>xf_[i+1])
            {
                xfmax=i;
                break;
            }
        }

        //balancing of the areas ahead of the front and below the curve
        bool a = true;
        double A1;
        double A2;
        double b;
        int xhelp2=0;

        while (a)
        {
            if (this->diffProblem.soil().relPermFlag(dummyGlobal_, dummyElement_, dummyLocal_) == Matrix2p<Grid, Scalar>::linear)
            break;

            A1=0;

            for (int i=0; i<=xhelp-1; i++)
            {
                A1+=(SatVec_[i]-Swr_+SatVec_[i+1]-Swr_)*0.5*(xf_[i+1]-xf_[i]);
            }

            A2=0;

            for (int i=xhelp; i<=xfmax-1; i++)
            {
                A2+=(SatVec_[xfmax]-SatVec_[i]+SatVec_[xfmax]-SatVec_[i+1])*0.5*(xf_[i+1]-xf_[i]);
            }

            b=xf_[xfmax];
            xhelp2=xfmax;

            while (b>xf_[xhelp])
            {
                xhelp2+=1;
                b=xf_[xhelp2];
            }

            for (int i=xfmax; i<=xhelp2; i++)
            {
                A2+=(SatVec_[i]-SatVec_[xfmax]+SatVec_[i+1]-SatVec_[xfmax])*0.5*(xf_[i]-xf_[i+1]);
            }

            xhelpoldold=xhelpold;
            xhelpold=xhelp;

            if (fabs(A1)>fabs(A2))
            {
                xhelp = xhelp - 1;
            }

            if (fabs(A1)<fabs(A2))
            {
                xhelp = xhelp +1;
            }

            if (xhelp == xhelpoldold)
            {
                a=false;
            }
        }

        // std::cout<<"xf_ = "<<xf_<<s                    Scalar oneminusSnr = 1-Snr_;td::endl;
        //iterate over vertices and get analytical saturation solution

        const GridView& gridView(this->grid.levelView(0));
        ElementIterator eendit = gridView.template end<0>();

        for (ElementIterator it = gridView.template begin<0>(); it!= eendit; ++it)
        {
            //find index of current vertex
            int index = this->mapper.map(*it);

            //get global coordinate of vertex
            Dune::FieldVector<Scalar,dimworld> global = it->geometry().corner(0);

            //account for linear material law
            if (this->transProblem.soil().relPermFlag(dummyGlobal_, dummyElement_, dummyLocal_) == Matrix2p<Grid, Scalar>::linear)
            {
                if (global[0]<=xf_[1])
                {
                    analyticSolution[index] = 1-Snr_;
                }
                if (global[0]>xf_[1])
                {
                    analyticSolution[index] = Swr_;
                }
            }

            //non-linear material law

            else
            {
                //find x_f next to global coordinate of the vertex
                int xnext = 0;
                for (int i=intervalNum_; i>=0; i--)
                {
                    if (global[0]<xf_[i])
                    {
                        xnext = i;
                        break;
                    }
                }

                //account for the area not yet reached by the front
                if (global[0]> xf_[xhelp2])
                {
                    analyticSolution[index] = Swr_;
                    continue;
                }

                if (global[0] <= xf_[xhelp2])
                {
                    analyticSolution[index] = SatVec_[xnext];
                    continue;
                }
            }
        }

        //call function to calculate the saturation error
        calcSatError(this->variables.saturation);

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
        VTKWriter<typename Grid::LeafGridView> vtkwriter(grid.leafView());
        char fname[128];
        sprintf(fname, "%s-%05d", name, k);
        vtkwriter.addCellData(this->transProblem.variables.saturation, "saturation");
        vtkwriter.addCellData(this->transProblem.variables.pressure, "total pressure p~");
        vtkwriter.addCellData(analyticSolution, "saturation (exact solution)");
        vtkwriter.addCellData(error, "error");
        vtkwriter.write(fname, VTKOptions::ascii);

        return;
    }

    //! Construct an IMPES object.
    IMPESBLAnalytic(Diffusion& diffusion, Transport& transport, int flag = 0, int nIt = 2,
            Scalar maxDef = 1e-5, Scalar om = 1, Scalar cfl = 1, Scalar totalvelocity = 3e-7):
    IMPES<Grid, Diffusion, Transport, VC> (diffusion, transport, flag, nIt, maxDef, om),
    grid(this->variables.grid), mapper(grid.levelView(grid.maxLevel())), analyticSolution(0), error(0),
    elementVolume(0), size_(mapper.size()),CFL_(cfl), vtot_(totalvelocity), dummyElement_(*(grid.levelView(grid.maxLevel()).template begin<0>())), dummyLocal_(LocalPosition(1)), dummyGlobal_(GlobalPosition(1))
    {
        initializeAnalytic();
        prepareAnalytic();
    }

protected:
    const Grid& grid;
    ElementMapper mapper;
    BlockVector analyticSolution;
    BlockVector error;
    BlockVector elementVolume;

private:
    Scalar time_;
    int size_;
    Scalar CFL_;
    Scalar Swr_;
    Scalar Snr_;
    Scalar vtot_;
    enum
    {   intervalNum_ = 1000, pointNum_ = intervalNum_+1};
    FieldVector<Scalar, pointNum_> SatVec_;
    FieldVector<Scalar,pointNum_> fractionalW_;
    FieldVector<Scalar,pointNum_> dfwdsw_;
    FieldVector<Scalar,pointNum_> xf_;
    int dfwdswmax_;
    const Element& dummyElement_;
    const LocalPosition& dummyLocal_;
    const GlobalPosition& dummyGlobal_;

};
}
#endif
