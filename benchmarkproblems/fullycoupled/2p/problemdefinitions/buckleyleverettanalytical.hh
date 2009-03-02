#ifndef DUNE_BUCKLEYLEVERETTANALYTICAL_HH
#define DUNE_BUCKLEYLEVERETTANALYTICAL_HH

#include"dumux/twophase/exsolution.hh"

/**
 * @file
 * @brief Class adding the analytical solution of the Buckley Leverett problem
 * @author Markus Wolff
 */

namespace Dune {

template<class G, class RT> class BLWithAnalytical : public ExSolution<G, RT>,
                                                     public
    BuckleyLeverettProblem<G, RT> {

    enum {BrooksCorey=0};
    enum {n=G::dimension,dimworld = G::dimensionworld,m=2};
    typedef typename G::ctype DT;
    typedef typename G::template Codim <n>:: LeafIterator Iterator;
    typedef BlockVector<FieldVector<RT, m> > BVu;
    typedef BlockVector<FieldVector<RT, 1> > BV;

public:

    BLWithAnalytical(const G &g, DeprecatedTwoPhaseRelations &law = *(new DeprecatedLinearLaw), const FieldVector<DT,n> Left = 0,
                     const FieldVector<DT,n> Right = 0, int chooselaw = BrooksCorey,
                     RT totalvelocity = 3e-7) :
        ExSolution<G, RT>(g), BuckleyLeverettProblem<G, RT>(law, Left, Right,
                                                            chooselaw, true), vtot(totalvelocity) {
        prepareanalytical();
    }

    void prepareanalytical() {
        Swr=this->materialLaw_.wettingPhase.Sr();
        Snr=this->materialLaw_.nonwettingPhase.Sr();
        Porosity=this->Porosity_;

        time=0;

        SatVec=Swr;
        for (int i=1; i<pointnum; i++) {
            SatVec[i]=SatVec[i-1]+(1-Snr-Swr)/intervalnum;
        }
        //std::cout<<"SatVec = "<<SatVec_<<std::endl;
        for (int i=0; i<pointnum; i++) {
            fractionalW[i] = this->materialLaw_.fractionalW(SatVec[i]);
        }
        //std::cout<<"fractionalW = "<<fractionalW_<<std::endl;
        dfwdsw=0;
        for (int i=1; i<intervalnum; i++) {
            dfwdsw[i]=(fractionalW[i+1]-fractionalW[i-1])/(SatVec[i+1]
                                                           -SatVec[i-1]);
        }
        //std::cout<<"dfwdsw = "<<dfwdsw_<<std::endl;
        for (int i=0; i<pointnum; i++) {
            if (dfwdsw[i]>dfwdsw[i+1]) {
                dfwdswmax = i;
                break;
            }
        }
        return;
    }

    void updateExSol(double &dt, BVu &approxSol) {

        time+=dt;

        //possition of the fluid front
        xf=0;
        for (int i=0; i<pointnum; i++) {
            xf[i]=vtot*time/Porosity*dfwdsw[i];
        }
        //std::cout<<"xf = "<<xf_<<std::endl;
        int xhelp=pointnum/3;
        int xhelpold = 0;
        int xhelpoldold = 0;
        int xfmax = 0;
        ;

        //position of maximum xf
        for (int i=0; i<pointnum; i++) {
            if (xf[i]>xf[i+1]) {
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

        while (a) {
            if (this->materialLaw_.isLinear())
                break;

            A1=0;
            for (int i=0; i<=xhelp-1; i++) {
                A1+=(SatVec[i]-Swr+SatVec[i+1]-Swr)*0.5*(xf[i+1]-xf[i]);
            }
            A2=0;
            for (int i=xhelp; i<=xfmax-1; i++) {
                A2+=(SatVec[xfmax]-SatVec[i]+SatVec[xfmax]-SatVec[i+1])*0.5
                    *(xf[i+1]-xf[i]);
            }
            b=xf[xfmax];
            xhelp2=xfmax;
            while (b>xf[xhelp]) {
                xhelp2+=1;
                b=xf[xhelp2];
            }
            for (int i=xfmax; i<=xhelp2; i++) {
                A2+=(SatVec[i]-SatVec[xfmax]+SatVec[i+1]-SatVec[xfmax])*0.5
                    *(xf[i]-xf[i+1]);
            }
            xhelpoldold=xhelpold;
            xhelpold=xhelp;
            if (fabs(A1)>fabs(A2)) {
                xhelp = xhelp - 1;
            }
            if (fabs(A1)<fabs(A2)) {
                xhelp = xhelp +1;
            }
            if (xhelp == xhelpoldold) {
                a=false;
            }
        }
        //iterate over vertices and get analytical saturation solution
        for (Iterator it = this->grid.template leafbegin<n>(); it
                 != this->grid.template leafend<n>(); ++it) {

            //find index of current vertex
            int index = this->mapper.map(*it);

            //get global coordinate of vertex
            Dune::FieldVector<DT,dimworld> global = it->geometry().corner(0);

            //account for linear material law
            if (this->materialLaw_.isLinear()) {
                if (global[0]<=xf[1]) {
                    RT oneminusSnr = 1-Snr;
                    this->uExInVertex(oneminusSnr, index, 1);
                }
                if (global[0]>xf[1])
                    this->uExInVertex(Swr, index, 1);
            }

            //non-linear material law
            else {
                //find x_f next to global coordinate of the vertex
                int xnext = 0;
                for (int i=intervalnum; i>0; i--) {
                    if (global[0]<xf[i]) {
                        xnext = i;
                        break;
                    }
                }

                //account for the area not yet reached by the front
                if (global[0] > xf[xhelp2]) {
                    this->uExInVertex(Swr, index, 1);
                    continue;
                }

                if (global[0] <= xf[xhelp2]) {
                    this->uExInVertex(SatVec[xnext], index, 1);
                    continue;
                }
            }
        }
        //get Sw from Sn -> pwSn-formulation is used <- for error calculation
        BVu Approx(this->mapper.size());
        for (int i = 0; i < this->mapper.size(); i++) {
            Approx[i][1] = 1-approxSol[i][1];
        }
        //call function to calculate the saturation error
        this->calcSatError(Approx);

        return;
    }

    RT uExOutVertex(int &ElementIndex, int VariableIndex) const {
        return this->uEx[ElementIndex][VariableIndex];
    }

private:
    RT Swr;
    RT Snr;
    RT Porosity;

    RT time;
    RT vtot;
    FieldVector<RT,4> parameters;
    enum {intervalnum = 1000, pointnum = intervalnum+1};
    FieldVector<RT, pointnum> SatVec;
    FieldVector<RT,pointnum> fractionalW;
    FieldVector<RT,pointnum > dfwdsw;
    FieldVector<RT,pointnum> xf;
    int dfwdswmax;
};

}
#endif
