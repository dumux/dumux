#ifndef DUNE_CHTRANSPORT_HH
#define DUNE_CHTRANSPORT_HH

#include <dune/common/helpertemplates.hh>
#include <dune/common/typetraits.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/bvector.hh>
#include "dumux/transport/transport.hh"
#include "dumux/transport/fv/numericalflux.hh"
#include "dumux/transport/fv/diffusivepart.hh"
#include "dumux/transport/transportproblem.hh"

namespace Dune
{

//! \ingroup transport
//! The finite volume model for the solution of the transport equation
template<class Grid, class Scalar, class VC, class Problem = TransportProblem<
        Grid, Scalar, VC> > class ChTransport: public Transport<Grid, Scalar,
        VC, Problem>
{

    template<int dim> struct ElementLayout
    {
        bool contains(Dune::GeometryType gt)
        {
            return gt.dim() == dim;
        }
    };

    //structure of an element in the dynamic grid
    //Sat:    saturation
    //x[i,j]: position of the i-front in time t_j. i,j \in {0,1}
    //m[i]:   propagation velocity of the i-front. i \in {0,1}
    struct ChNode
    {
        Scalar Sat;
        Dune::FieldVector<Dune::FieldVector<Scalar,2>,2> x;
        Dune::FieldVector<Scalar,2> m;
    };

    //node of the streamline-list
    //i:     index of the element the streamline crossed
    //xi[i]: streamline parameter of the entrance and exit points
    struct slNode
    {
        int index;
        Dune::FieldVector<Scalar,2> xi;
    };

    enum
    {
        dim = Grid::dimension
    };
    enum
    {
        dimWorld = Grid::dimensionworld
    };
    typedef BlockVector<Dune::FieldVector<Scalar,1> > PressType;
    typedef BlockVector<FieldVector<FieldVector<Scalar, Grid::dimension> , 2
            * Grid::dimension> > VelType;
typedef    typename Grid::Traits::template Codim<0>::Entity Entity;
    typedef typename Grid::LevelGridView GV;
    typedef typename GV::IndexSet IS;
    typedef typename GV::template Codim<0>::Iterator ElementIterator;
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<GV,ElementLayout> EM;
    typedef typename Grid::template Codim<0>::EntityPointer EntityPointer;
    typedef typename GV::IntersectionIterator IntersectionIterator;
    typedef typename std::template list<ChNode>::iterator ListIterator;
    typedef typename Grid::ctype ct;
    typedef Dune::FieldVector<Scalar,dim> VType;
    typedef BlockVector< Dune::FieldVector<Scalar,dim> > SlopeType;

    typedef Dune::FieldVector<Scalar,dim> LocalPosition;
    typedef Dune::FieldVector<Scalar,dimWorld> GlobalPosition;

public:
    typedef BlockVector< Dune::FieldVector<Scalar,1> > RepresentationType;
    /*!
     *  \param t time
     *  \param dt time step size to estimate
     *  \param update vector to be filled with the update
     *
     *  This method calculates the update vector of \f$\text{div}\, (f_\text{w} \boldsymbol{v}_t)\f$
     *  by using a characteristic method. The total velocity \f$\boldsymbol{v}_t)\f$ has to be given
     *  by the block vector \a velocity, containing values at each cell face. The fractional flow
     *  function \f$f_\text{w}\f$ is realized by the numerical flux function \a numericalFlux such
     *  that the \f$i\f$-th entry of the vector update is obtained by calculating
     *  \f[ \sum_{j \in \mathcal{N}(i)} v_{ij}g(S_i, S_j) - v_{ji}g(S_j, S_i), \f]
     *  where \f$\mathcal{N}(i)\f$ denotes the set of neighbors of the cell \f$i\f$ and
     *  \f$g\f$ stands for \a numericalFlux. The normal velocities \f$v_{ij}\f$ and \f$v_{ji}\f$
     *  are given by
     *  \f{align*} v_{ij} = \int_{\Gamma_{ij}} \max(\boldsymbol{v}_\text{t}{\cdot}\boldsymbol{n}_{ij}, \, 0), \qquad
     *  v_{ji} = \int_{\Gamma_{ij}} \min(\boldsymbol{v}_\text{t}{\cdot}\boldsymbol{n}_{ij}, \, 0), \f}
     *  where \f$\boldsymbol{n}_{ij}\f$ denotes the unit normal vector from cell \f$i\f$ towards cell \f$j\f$.
     *
     */
    int update(const Scalar t, Scalar& dt, RepresentationType& updateVec, Scalar& cFLFac);

    void initialTransport();

    /*! @brief constructor
     *
     * @param g a DUNE grid object
     * @param prob an object of class TransportProblem or derived
     * @param lev the grid level on which the Transport equation is to be solved.
     * @param diffPart an object of class DiffusivePart or derived. This determines the diffusive flux incorporated in the transport.
     * @param rec flag to switch on linear reconstruction (second order TVD)
     * @param amax alphamax parameter for slope limiter in TVD
     * @param numFl an object of class Numerical Flux or derived
     */
    ChTransport(Grid& g, Problem& prob,
            DiffusivePart<Grid,Scalar>& diffPart = *(new DiffusivePart<Grid, Scalar>), bool rec = false,
            double amax = 0.8, const NumericalFlux<Scalar>& numFl = *(new Upwind<Scalar>), int K=1000) :
    Transport<Grid, Scalar, VC, Problem>(g, prob),
    elementmapper(g.levelView(this->level())), reconstruct(rec),
    numFlux(numFl), diffusivePart(diffPart), alphamax(amax),K(K)
    {}

private:

    int calcm(GlobalPosition globalPos, ElementIterator eIt, LocalPosition localPos);
    int init(const Scalar dt,const Scalar lambda, EntityPointer& eIt);
    int calcStreamline(const Scalar dt,const Scalar lambda, EntityPointer& eIt);
    int solveRP(Scalar& mch, std::list<ChNode>& chrar,Scalar sat0,Scalar sat1, GlobalPosition globalPos, ElementIterator eIt, LocalPosition localPos);
    Scalar linearflux(Scalar sat, GlobalPosition globalPos, ElementIterator eIt, LocalPosition localPos);
    int fronttracking(const Scalar& dt, Scalar v, GlobalPosition globalPos, ElementIterator eIt, LocalPosition localPos);
    int gelb(RepresentationType& updateVec,int indexi);
    int blau(RepresentationType& updateVec,int indexi);
    int rot(RepresentationType& updateVec,RepresentationType& dxtt);
private:
    EM elementmapper;
    bool reconstruct;
    const NumericalFlux<Scalar>& numFlux;
    const DiffusivePart<Grid, Scalar>& diffusivePart;
    double alphamax;


    int K;
    std::list<ChNode> ch;
    std::list<slNode> sl;
};

//approximation of the fractional flow in a piecewise linear function
//evaluation at the parameter Scalar sat
template<class Grid,class Scalar,class VC, class Problem>
Scalar ChTransport<Grid,Scalar,VC, Problem>::linearflux(Scalar sat, GlobalPosition globalPos, ElementIterator eIt, LocalPosition localPos)
{
    Scalar satA=-1;
    //search of satA with sat \in [satA, satA+i/K]
    for(int i=0;i<K;++i)
    if(K*sat-i>=0)
    satA=i*1.0/K;

    //calculation of f(satA), f(satA+i/K)
    double fa=this->transProblem.materialLaw().fractionalW(satA, globalPos, *eIt, localPos);
    double fb=this->transProblem.materialLaw().fractionalW(satA+1.0/K, globalPos, *eIt, localPos);

    //evaluation at sat
    double fsat=fa+(sat-satA)*K*(fb-fa);
    return(fsat);
}

//Riemann problem solver
//mch:   propagation velocity of the front, if the solution is a shock
//chrar: list of dynamic grid nodes, if the solution is a rarefaction wave
//sat0:  value on the left side of the discontinuity
//sat1:  value on the right side of the discontinuity
template<class Grid,class Scalar,class VC, class Problem>
int ChTransport<Grid,Scalar,VC, Problem>::solveRP(Scalar& mch, std::list<ChNode>& chrar,Scalar sat0,Scalar sat1, GlobalPosition globalPos, ElementIterator eIt, LocalPosition localPos)
{
    //clear temporary dynamic grid list
    chrar.clear();

    std::vector<Scalar> sat(2,0);
    std::vector<Scalar> m(2,0);
    sat[0]=sat0;
    sat[1]=sat1;

    //left and right value are identical
    //function returns the propagation speed of the characteristic
    if(fabs(sat1-sat0)<1e-5)
    {
        mch=0;
        return (linearflux(sat0, globalPos, eIt, localPos));
    }

    //Rankine Hugoniot condition
    mch=(linearflux(sat1, globalPos, eIt, localPos)-linearflux(sat0, globalPos, eIt, localPos))/(sat1-sat0);
    m[0]=mch;
    m[1]=mch;

    //according to the low value on the left or right side
    //a convex or concave envelope will be constructed
    int low;
    if(sat[0]<sat[1])
    low=0;
    else
    low=1;
    int high=1-low;

    Scalar sl=sat[low];
    Scalar sh = sat[high];
    double fsl=linearflux(sl, globalPos, eIt, localPos);
    double fsh = linearflux(sh, globalPos, eIt, localPos);

    //construction of a list with all nodes between the low and high value
    int il=0,ih=0;

    for(int i=1;i<K;++i)
    {
        if(i-K*sat[low]>0 && il==0)
        il=i;
        if(K*sat[high]-i>0)
        ih=i;
    }

    //if there are no nodes the solution is a shock front
    if(il==0 || il> ih)
    return 0;

    //one node between sat0 and sat1
    if(il==ih)
    {
        Scalar satil = this->transProblem.materialLaw().fractionalW(il*1.0/K, globalPos, *eIt, localPos);
        //Rankine Hugoniot condition: shockfront
        if(sat0> sat1 && (fsl+mch*(il*1.0/K-sl) >= satil))
        return 0;
        else if(sat0 < sat1 && (fsl+mch*(il*1.0/K-sl) <= satil))
        return 0;
        //rarefaction wave approximation with two shock fronts

        else
        {
            m[low]=(satil-fsl)*K/(il-sl*K);
            m[high]=(satil-fsh)*K/(il-sh*K);

            ChNode cn;
            cn.Sat=il*1.0/K;
            cn.m[high]=m[high];
            cn.m[low]=m[low];
            chrar.push_front(cn);
        }
    }

    //more than one nodes between sat0 and sat1

    else
    {
        //len : number of the nodes
        int len=ih-il+1;

        //(b,fb) are the nodes and interpolation values between sat0 and sat1
        std::vector<double> b(len,0);
        std::vector<double> fb(len,0);
        for(int i=0;i<len;++i)
        {
            b[i]=(il+i)*1.0/K;
            fb[i]=this->transProblem.materialLaw().fractionalW(b[i], globalPos, *eIt, localPos);
        }

        int j1 = 0;
        int j2 = len-1;
        //the lower value is on the right side:
        //search of b[j] with [sat[low] b[j] b[j+1] .. b[len-1] sat[high]
        //build the convex envelope
        if(low == 1)
        {
            double n=(fb[1]-fsl)/(b[1]-sl);
            while(fsl+n*(b[j1]-sl) >= fb[j1])
            {
                j1++;
                if(j1==j2)
                n=0;
                else
                n=(fb[j1+1]-fsl)/(b[j1+1]-sl);
            }

            //there is no b[j]: shock front with propagation speed mch
            //(Rankine Hugoniot)
            if(j1==j2 && fsl+mch*(b[j1]-sl) >= fb[j1])
            return 0;
        }
        else
        {
            double n = (fb[j2-1]-fsh)/(b[j2-1]-sh);
            while(fsh+n*(b[j2]-sh) <= fb[j2])
            {
                j2--;
                if(j2==j1)
                n = 0;
                else
                n = (fb[j2-1]-fsh)/(b[j2-1]-sh);
            }
            if(j2 == j1 && fsh+mch*(b[j2]-sh) <= fb[j2])
            return 0;
        }

        m[low]=(fb[j1]-fsl)/(b[j1]-sl);
        m[high]=(fb[j2]-fsh)/(b[j2]-sh);

        ChNode c0;
        c0.Sat=b[j1];
        c0.m[high]=m[high];
        c0.m[low]=m[low];
        chrar.push_front(c0);
        //konstruiere neue Zellen für die zusätzlichen Schocks welche die
        //Verdünnungswelle approximieren
        //  std::cout<<b[j1]<<" ; "<<b[j2]<<std::endl;
        if(j1!=j2)
        {
            chrar.front().m[high]=(fb[j1+1]-fb[j1])*K;
            j1++;
            if(low == 1)
            while(j1<j2)
            {
                ChNode cn;
                cn.Sat=b[j1];
                cn.m[0]=(fb[j1+1]-fb[j1])*K;
                cn.m[1]=chrar.front().m[0];
                if(fabs(cn.m[0]-cn.m[1])>1e-5)
                chrar.push_front(cn);
                j1++;
            }
            else
            while(j1<j2)
            {
                ChNode cn;
                cn.Sat=b[j1];
                cn.m[0]=chrar.back().m[1];
                cn.m[1]=(fb[j1+1]-fb[j1])*K;
                if(fabs(cn.m[0]-cn.m[1])>1e-5)
                chrar.push_back(cn);
                j1++;
            }

            ChNode clen;
            clen.Sat=b[j2];
            clen.m[high]=m[high];

            if(low==1)
            {
                if(fabs(clen.m[0]-clen.m[1])>1e-5)
                {
                    clen.m[1]=chrar.front().m[0];
                    chrar.push_front(clen);
                }
            }
            else
            if(fabs(clen.m[0]-clen.m[1])>1e-5)
            {
                clen.m[0]=chrar.back().m[1];
                chrar.push_back(clen);
            }
        }
    }

    if(chrar.size()==1 && fabs(chrar.front().m[0]-chrar.back().m[1])<1e-5)
    chrar.clear();
    return 0;
}

template<class Grid,class Scalar,class VC, class Problem>
int ChTransport<Grid,Scalar,VC, Problem>::calcm(GlobalPosition globalPos, ElementIterator eIt, LocalPosition localPos)
{

    //berechne Steigungen der Charakteristiken
    ListIterator lit=ch.begin();
    ListIterator lit2(lit);
    lit2++;
    while(lit2!=ch.end())
    {
        std::list<ChNode> chrar;
        Scalar m;

        solveRP(m,chrar,lit->Sat,lit2->Sat,globalPos, eIt, localPos);

        if(chrar.empty())
        {
            if(m!=0)
            //Schock

            {
                lit->m[1]=m;
                lit2->m[0]=m;
            }
            else
            //Sättigungen sind identisch, Elemente können vereinigt werden
            //(kann in Abhängigkeit vom Ausgangselement vorkommen)

            {
                lit2->m[0]=lit->m[0];
                lit2->x[0][0]=lit->x[0][0];
                lit=ch.erase(lit);

                //  ch.front().x[0][0]=ch.front().x[0][1];
                //ch.back().x[0][1]=ch.back().x[0][0];
            }
        }
        else
        //Verdünnungswelle, berechne zusätzliche Schocks

        {
            lit->m[1]=chrar.front().m[0];
            lit2->m[0]=chrar.back().m[1];
            for(ListIterator kit=chrar.begin();kit!=chrar.end();++kit)
            for(int q=0;q<2;++q)
            kit->x[0][q]=lit2->x[0][0];

            ch.insert(lit2,chrar.begin(),chrar.end());
        }
        ch.front().m[0] = solveRP(m,chrar,ch.front().Sat,ch.front().Sat,globalPos, eIt, localPos);
        ch.back().m[1] = solveRP(m,chrar,ch.back().Sat,ch.back().Sat,globalPos, eIt, localPos);
        //      ch.front().m[0] = ch.front().m[1];
        //ch.back().m[1]=ch.back().m[0];

        lit=lit2;
        lit2++;
    }
    return 0;
}

template<class Grid,class Scalar,class VC, class Problem>
int ChTransport<Grid,Scalar,VC, Problem>::init(const Scalar dt,const Scalar lambda, EntityPointer& eIt)
{
    ch.clear();
    sl.clear();

    int indexi = elementmapper.map(*eIt);

    const GV& gridView = this->grid_.levelView(this->level());

    //Bilde Ausgangselement in der Charakteristikenliste
    ChNode i;
    i.Sat=this->transProblem.variables.saturation[indexi];

    //setze Geschwindigkeit und x-Koordinaten für das Ausgangselement
    IntersectionIterator endis = gridView.template iend(*eIt);
    for(IntersectionIterator is=gridView.template ibegin(*eIt);is!=endis;++is)
    {
        GeometryType gtf = is->geometryInInside().type();
        const Dune::FieldVector<ct,dim-1>& facelocal = Dune::ReferenceElements<Scalar,dim-1>::general(gtf).position(0,0);
        const VType& faceglobal = is->geometry().global(facelocal);
        //      const VType& faceInElem = eIt->geometry().local(faceglobal);
        int numberInSelf = is->indexInInside();

        //linke Kante
        if(numberInSelf == 0)
        i.x[0][0] = faceglobal[0];
        else
        i.x[0][1] = faceglobal[0];
    }
    ch.push_front(i);

    slNode s;
    s.index=indexi;
    s.xi=ch.front().x[0];
    sl.push_front(s);

    //std::vector<double> gdat(2,0);
    //beginne beim Ausgangselement und laufe solange nach links bzw rechts, bis eine Randkante kommt
    //trifft man auf einen Sättigungssprung, so nehme dies in die Charakteristikenliste mit auf
    for(int p=0;p<2;++p)
    {
        int abb=1;
        EntityPointer ep=eIt;
        Scalar satI=i.Sat;

        while(abb>0)
        {
            ChNode j;
            slNode s;
            s.index=indexi;

            IntersectionIterator endis = gridView.template iend(*ep);
            for(IntersectionIterator is=gridView.template ibegin(*ep);is!=endis;++is)
            {
                GeometryType gtf = is->geometryInInside().type();
                const Dune::FieldVector<ct,dim-1>& facelocal = Dune::ReferenceElements<Scalar,dim-1>::general(gtf).position(0,0);
                const VType& faceglobal = is->geometry().global(facelocal);
                //      const VType& faceInElem = eIt->geometry().local(faceglobal);
                int numberInSelf = is->indexInInside();

                s.xi[numberInSelf]=faceglobal[0];

                if (numberInSelf == p) // left face

                {
                    j.x[0][0] = faceglobal[0];
                    j.x[0][1] = j.x[0][0];
                    //Scalar satJ=0;
                    if(is->neighbor())
                    {
                        ep=is->outside();
                        indexi=elementmapper.map(*ep);
                        j.Sat=this->transProblem.variables.saturation[indexi];
                    }
                    else
                    {
                        //gdat[p] = faceglobal[0];

                        const Dune::FieldVector<ct,dim>& facelocalDim = Dune::ReferenceElements<ct,dim>::general(gtf).position(numberInSelf,1);
                        j.Sat =this->transProblem.dirichletSat(faceglobal,*ep,facelocalDim);
                        abb=-1;
                    }

                    if(fabs(satI-j.Sat)>1e-5)
                    {
                        if(p==0)
                        {
                            ch.front().x[0][0]=j.x[0][1];
                            ch.push_front(j);
                        }
                        else
                        {
                            ch.back().x[0][1]=j.x[0][0];
                            ch.push_back(j);
                        }
                        satI=j.Sat;
                    }
                }
            }
            if(p==0 && s.index != sl.front().index)
            sl.push_front(s);
            else
            if(p == 1 && s.index != sl.back().index)
            sl.push_back(s);
        }//end while
    }//end for

    ch.front().x[0][0] -= lambda*dt;
    ch.back().x[0][1] += lambda*dt;

    if(ch.size()<2)
    {
        ch.clear();
        sl.clear();
    }

    return 0;

}

template<class Grid,class Scalar,class VC, class Problem>
int ChTransport<Grid,Scalar,VC, Problem>::calcStreamline(const Scalar dt,const Scalar lambda, EntityPointer& eIt)
{
    const GV& gridView = this->grid_.levelView(this->level());

    ch.clear();
    sl.clear();

    EntityPointer ep(eIt);
    Scalar xi = 0;

    // cell geometry type
    Dune::GeometryType gt = eIt->geometry().type();
    // cell center in reference element
    const LocalPosition
        &localPos = Dune::ReferenceElements<Scalar,dim>::general(gt).position(0, 0);

    //
    const GlobalPosition& globalPos = eIt->geometry().global(localPos);

    int indexi = elementmapper.map(*eIt);
    const Dune::FieldVector<ct,dim>& local = Dune::ReferenceElements<ct,dim>::general(eIt->geometry().type()).position(0,0);
    Dune::FieldVector<ct,dim> global = eIt->geometry().global(local);
    VType xq(global);

    //Bilde Ausgangselement, ausgehend von dem Element dessen Sättigung berechnet werden soll
    ChNode i;
    i.Sat=this->transProblem.variables.saturation[indexi];
    i.x[0][0] = xi;
    i.x[0][1] = xi;
    ch.push_front(i);
    slNode j;
    j.index=indexi;
    j.xi[0]=xi;
    j.xi[1] = xi;
    sl.push_front(j);

    //Es wird so lange die Stromlinie zurückverfolgt, wie ein Partikel in der Zeit t zurücklegen konnte.
    //lambda stellt die höchste "Geschwindigkeit" eines Partikels dar ( max(f'(u)) )
    while(xi> -lambda*dt)
    {
        //std::cout<<xi<<std::endl;
        Dune::FieldVector<VType,2*dim> xx;
        Dune::FieldVector<VType,2*dim> vx;

        //Setzen der Positionen und Geschwindigkeiten an den Elementkanten
        IntersectionIterator endis = gridView.template iend(*ep);
        for(IntersectionIterator is=gridView.template ibegin(*ep);is!=endis;++is)
        {
            GeometryType gtf = is->geometryInInside().type();
            const Dune::FieldVector<ct,dim-1>& facelocal = Dune::ReferenceElements<Scalar,dim-1>::general(gtf).position(0,0);
            const VType& faceglobal = is->geometry().global(facelocal);
            //      const VType& faceInElem = eIt->geometry().local(faceglobal);
            int numberInSelf = is->indexInInside();

            xx[numberInSelf] = faceglobal;

            vx[numberInSelf] = this->transProblem.variables.vTotal(*ep,numberInSelf);
            vx[numberInSelf] /= this->transProblem.soil().porosity(globalPos, *eIt, localPos);
        }
        //Berechnung des Geschwindigkeitsgradienten und der Partikelgeschwindigkeit
        VType Ax(0);
        VType vxq(0);

        for(int i=0;i<dim;++i)
        {
            Ax[i] = (vx[2*i+1][i]-vx[2*i][i])/(xx[2*i+1][i]-xx[2*i][i]);
            vxq[i] = Ax[i]*(xq[i]-xx[2*i][i])+vx[2*i][i];

        }

        //Berechnung aus welcher Richtung das Partikel ins Element gelangt ist. Für negativen Gradienten/negative
        //Geschwindigkeit kam das Partikel von rechts/oben, andernfalls von links/unten.
        //kant speichert den entsprechenden Kantenindex pro Richtung.
        //dtx berechnet pro Richtung die Zeit, welche das Partikel von der entsprechenden Kante zur momentanen
        //Partikelposition gebraucht hat.
        VType dtx(-2*lambda*dt);
        Dune::FieldVector<int,dim> kant(-1);

        for(int i = 0;i<dim;++i)
        {
            if(fabs(Ax[i])<1e-5)
            //konstante Geschwindigkeit in i-Richtung innerhalb des Elements

            {
                if(vxq[i]<0)
                //v < 0.. Partikel muss rechts/oben eingetreten sein
                dtx[i] = (xx[2*i+1][i]-xq[i])/vxq[i];
                else if(vxq[i]>0)
                //Partikel muss links/unten eingetreten sein
                dtx[i] = (xx[2*i][i]-xq[i])/vxq[i];
            }
            else
            dtx[i] = 1./Ax[i]*log(vx[i*2+1][i]/vxq[i]);

            if(vxq[i]<0)
            kant[i] = 2*i+1;
            else if(vxq[i]>0)
            kant[i] = 2*i;
        }

        //Berechnung der Zeit von der Eintrittskante zum gewuenschten Position
        int m = 0; //Kantenindex
        Scalar dte = dtx[m]; //Time-of-Flight
        for(int i=1;i<dim;++i)
        if(dtx[i]>dte)
        {
            m = i;
            dte=dtx[i];
        }

        //Berechnung der Position der Eintrittsstelle
        int n = -1;
        if(dim>1)
        {
            for(int i = 0;i<dim;++i)

            if(fabs(dtx[i]-dte)<1e-7)
            {
                xq[i]=xx[kant[i]][i];
                if(i!=m)
                n = 17;
            }
            else
            if(fabs(Ax[i])>1e-5)
            xq[i] = xq[i]+1./Ax[i]*(vxq[i]*exp(Ax[i]*dte)-vx[i*dim][i]);
            else
            xq[i] = xq[i]+vxq[i]*dte;
        }
        else
        //im 1D ist Eintrittstelle = Kante
        xq = xx[kant[m]];

        //Berechnung der Eintrittszeit, Bildung der neuen Knotenpunkte für die Listen
        xi = xi+dte;
        ch.front().x[0][0]=xi;
        ChNode j;
        j.x[0][0]=xi;
        j.x[0][1] = xi;
        sl.front().xi[0]=xi;
        slNode s;
        s.xi=j.x[0];
        s.index = -1;

        for(IntersectionIterator is = gridView.template ibegin(*ep);is!=endis;++is)
        {
            int indexInInside = is->indexInInside();
            if(indexInInside == kant[m])
            if(is->neighbor())
            {
                ep = is->outside();

                if(n> 0)
                {
                    for(int i=0;i<dim;++i)
                    if(i!=m && fabs(dtx[i]-dte)<1e-7)
                    {
                        IntersectionIterator endir = gridView.template iend(*ep);
                        for(IntersectionIterator ir = gridView.template ibegin(*ep);ir!=endir;++ir)
                        {
                            int nIS = ir->indexInInside();

                            if(nIS == kant[i])
                            if(ir->neighbor())
                            {
                                ep = ir->outside();
                                indexi = elementmapper.map(*ep);
                                j.Sat = this->transProblem.variables.saturation[indexi];

                                s.index=indexi;
                            }
                            else
                            {
                                GeometryType gtf = ir->geometryInInside().type();
                                const Dune::FieldVector<ct,dim-1>& facelocal = Dune::ReferenceElements<Scalar,dim-1>::general(gtf).position(0,0);
                                const VType& faceglobal = ir->geometry().global(facelocal);
                                const Dune::FieldVector<ct,dim>& facelocalDim = Dune::ReferenceElements<ct,dim>::general(gtf).position(indexInInside,1);

                                j.Sat = this->transProblem.dirichletSat(faceglobal,*ep,facelocalDim);
                                if(xi>-lambda*dt)
                                {
                                    xi = -2*lambda*dt;
                                    j.x[0][0]=xi;
                                }
                            }//end if ir->neighbor()
                        }//end ir
                    }//end i!=m
                }//if n>0

                else
                {
                    indexi = elementmapper.map(*ep);
                    j.Sat = this->transProblem.variables.saturation[indexi];

                    s.index=indexi;
                }
            }//if is->neighbor()

            else
            {
                GeometryType gtf = is->geometryInInside().type();
                const Dune::FieldVector<ct,dim-1>& facelocal = Dune::ReferenceElements<Scalar,dim-1>::general(gtf).position(0,0);
                const VType& faceglobal = is->geometry().global(facelocal);
                const Dune::FieldVector<ct,dim>& facelocalDim = Dune::ReferenceElements<ct,dim>::general(gtf).position(indexInInside,1);
                j.Sat = this->transProblem.dirichletSat(faceglobal,*ep,facelocalDim);

                if(xi>-lambda*dt)
                {
                    xi = -2*lambda*dt;
                    j.x[0][0]=xi;
                }
            }
        }//end is
        if(s.index != -1)
        sl.push_front(s);

        if(fabs(j.Sat-ch.front().Sat)>1e-5)
        ch.push_front(j);

    }//end while
    ch.front().x[0][0] = xi;
    sl.front().xi[0] = xi;

    if(ch.size()<2)
    {
        ch.clear();
        sl.clear();
    }

    return 0;
}

template<class Grid,class Scalar,class VC, class Problem>
int ChTransport<Grid,Scalar,VC, Problem>::fronttracking(const Scalar& dt, Scalar v, GlobalPosition globalPos, ElementIterator eIt, LocalPosition localPos)
{
    calcm(globalPos, eIt, localPos);

    //berechne Postion der Schocks zum  Zeitpunkt dt
    ListIterator lit;

    for(lit=ch.begin();lit!=ch.end();++lit)
    for(int q=0;q<2;++q)
    lit->x[1][q]=lit->x[0][q]+lit->m[q]*dt*v;
    //da eventuell Makrozeitschritte berechnet werden:
    //tau=Makrozeitpunkt, t1=vorherige Makrozeitpunkt,dtau=tau-t1
    //x[0] gibt die Position der Shocks zum Zeitpunkt t1 an
    //x[1] gibt die Position der Shocks zum Zeipunkt tau an


    Scalar t1=0;
    while(t1<dt)
    {
        Scalar dtau=dt;
        ListIterator merk;
        //berechne kleinsten Zeitschritt, bei dem sich die Shocks schneiden.
        //merk gibt an, zu welcher Zelle die entsprechenden Shocks gehören
        for(lit=ch.begin();lit!=ch.end();++lit)
        if(lit->x[0][0]<=lit->x[0][1] && lit->x[1][0]>lit->x[1][1] )
        {
            Scalar dch = lit->m[0]-lit->m[1];
            Scalar tt=dt;
            if(dch!=0)
            tt=(lit->x[0][1]-lit->x[0][0])/dch;
            if(tt<dtau)
            {
                dtau=tt;
                merk=lit;
            }
        }

        if(t1+dtau>=dt)
        dtau=dt-t1;

        //berechne Position der Shocks zum Zeitpunkt tau
        //xt=x+ch*n*dtau
        //setze x für nächsten Zeitschritt
        for(lit=ch.begin();lit!=ch.end();++lit)
        for(int p=0;p<2;++p)
        {
            lit->x[1][p]=lit->x[0][p]+lit->m[p]*dtau*v;
            lit->x[0][p]=lit->x[1][p];
        }
        //Verkleinerung des Zeitfensters in dem Schnitte auftreten können.
        t1=t1+dtau;

        if(t1>=dt)
        return 0;
        ListIterator chend = ch.end();
        //berechne neue Charakteristik an der Schnittstelle
        lit=merk;

        if(merk!=ch.begin() && merk!=--chend)
        {
            ListIterator merk2(lit);
            merk--; //merk ist nun lit-1
            merk2++; //merk ist nun lit+1

            std::list<ChNode> chrar;
            Scalar m;

            //   std::cout<<merk->Sat<<" ; "<<merk2->Sat<<std::endl;

            solveRP(m,chrar,merk->Sat,merk2->Sat,globalPos, eIt, localPos);

            //setze neue Steigungen
            if(chrar.empty())
            {
                if(m!=-0)
                //Schock

                {
                    merk->m[1]=m;
                    merk2->m[0]=m;
                    ch.erase(lit);

                    merk2->x[0][0]=0.5*(merk->x[0][1]+merk2->x[0][0]);
                    merk->x[0][1] = merk2->x[0][0];
                }
                else
                //Sättigungen identisch (kommt eigtl nicht vor)

                {
                    merk2->m[0]=merk->m[0];
                    merk2->x[0][0]=merk->x[0][0];

                    ch.erase(lit);
                }
            }
            else
            //Verdünnungswelle

            {
                merk->m[1]=chrar.front().m[0];
                merk2->m[0]=chrar.back().m[1];
                for(ListIterator kit=chrar.begin();kit!=chrar.end();++kit)
                for(int p=0;p<2;++p)
                kit->x[0][p]=merk2->x[0][0];

                ch.erase(lit);
                ch.insert(merk2,chrar.begin(),chrar.end());
            }
            chrar.clear();
        }
        else
        ch.erase(lit);
        //  solveRP(ch.front().m[0],chrar,ch.front().Sat,ch.front().Sat);
        //solveRP(ch.back().m[1],chrar,ch.back().Sat,ch.back().Sat);
        //ch.front().m[0]=ch.front().m[1];
        //ch.back().m[1]=ch.back().m[0];


        //berechne Position der Schocks zum Zeitpunkt dt
        for(lit=ch.begin();lit!=ch.end();++lit)
        for(int q=0;q<2;++q)
        lit->x[1][q]=lit->x[0][q]+lit->m[q]*(dt-t1)*v;
    }

    return 0;
}

template<class Grid, class Scalar, class VC, class Problem>
int ChTransport<Grid,Scalar,VC, Problem>::blau(RepresentationType& updateVec, int indexi)
{

    Scalar sat = 0;

    if(ch.front().x[1][0] <= 0 && 0 < ch.front().x[1][1])
    {
        //die Zelle liegt vor dem ersten Schock -> linke Randsättigung
        sat=ch.front().Sat;
    }
    else if(ch.back().x[1][0] < 0 && 0 <= ch.back().x[1][1])
    //die Zelle liegt nach dem letzten Schock -> rechte Randsättigung
    sat=ch.back().Sat;

    else
    {
        ListIterator lit=ch.begin();
        ++lit;
        ListIterator lit2(lit);
        ++lit2;
        while(lit2!=ch.end() && sat == 0)
        {
            Scalar xt0=lit->x[1][0];
            Scalar xt1=lit->x[1][1];

            if(xt0 <= 0 && 0 <= xt1)
            {
                Scalar xt05=0.5*(xt0+xt1);
                Scalar st05=lit->Sat;

                ListIterator merk(lit),merk2(lit);
                merk--;
                merk2++;
                Scalar sx0,sx1,xx0,xx1;

                //berechnet die Lösung an der Stelle x über lineare Interpolation der
                //Stützstellen des dualen Gitters

                if(0<xt05)
                {
                    if(merk->x[1][0]==ch.front().x[1][0])
                    {
                        xx0 = merk->x[1][1];
                        //  sx0 = 0.5*(lit->Sat+merk->Sat);
                        sx0=merk->Sat;
                    }
                    else
                    {
                        xx0 = 0.5*(merk->x[1][1] + merk->x[1][0]);
                        sx0 = merk->Sat;
                    }
                    sat=(st05*(0-xx0)-sx0*(0-xt05))/(xt05-xx0);
                }
                else
                {
                    if(merk2->x[1][1] == ch.back().x[1][1])
                    {
                        xx1=merk2->x[1][0];
                        //	      sx1=0.5*(lit->Sat+merk2->Sat);
                        sx1 = merk2->Sat;
                    }
                    else
                    {
                        xx1 = 0.5*(merk2->x[1][0] + merk2->x[1][1]);
                        sx1 = merk2->Sat;
                    }
                    sat=(sx1*(0-xt05)-st05*(0-xx1))/(xx1-xt05);
                }
            }//end if(xt0<= x < xt1)
            ++lit;
            ++lit2;
        }//end while
    }
    updateVec[indexi]=sat - this->transProblem.variables.saturation[indexi];

    return 0;
}

template<class Grid, class Scalar, class VC, class Problem>
int ChTransport<Grid,Scalar,VC, Problem>::rot(RepresentationType& updateVec,RepresentationType& dxtt)
{

    typedef typename std::template list<slNode>::iterator SlistIterator;
    for(SlistIterator slit=sl.begin();slit!=sl.end();++slit)
    {
        Scalar x0=slit->xi[0];
        Scalar x1=slit->xi[1];

        Scalar dx=fabs(x1-x0);
        Scalar sat=0;
        Scalar dxt=0;

        ListIterator lit = ch.begin();
        while(dxt<dx && lit!=ch.end())
        {
            Scalar xt0=lit->x[1][0];
            Scalar xt1=lit->x[1][1];

            if((xt0 <= x0 && x0 <= xt1) || (xt0 <= x1 && x1 <= xt1) || (x0 <= xt0 && xt0 <= x1)/*|| (x0 <= xt1 && xt1 <= x1)*/)
            {
                sat+=fabs(std::min(xt1,x1)-std::max(xt0,x0))*lit->Sat;
                dxt+=fabs(std::min(xt1,x1)-std::max(xt0,x0));
            }

            ++lit;
        }

        if(dxt>0)
        {
            updateVec[slit->index]+=sat;
            dxtt[slit->index]+=dxt;
        }
    }
    return 0;
}

template<class Grid, class Scalar, class VC, class Problem>
int ChTransport<Grid,Scalar,VC, Problem>::update(const Scalar t, Scalar& dt, RepresentationType& updateVec, Scalar& cFLFactor)
{
    const GV& gridView = this->grid_.levelView(this->level());

    if(t==0)
    // dt=8.64e6;
    // dt = 2.16e7;
    //dt = 1e8;
    dt = 4.32e7;

    // set update vector to zero
    updateVec = 0;
    RepresentationType dxtt(updateVec);
    dxtt=0;

    ElementIterator eendit = gridView.template end<0>();
    if(Grid::dimension == 1)
    {

        ElementIterator eIt = gridView.template begin<0>();

        // cell geometry type
        Dune::GeometryType gt = eIt->geometry().type();
        // cell center in reference element
        const LocalPosition
            &localPos = Dune::ReferenceElements<Scalar,dim>::general(gt).position(0, 0);

        //
        const GlobalPosition& globalPos = eIt->geometry().global(localPos);

        Scalar maxa = 0;
        for(int l=1;l<K;++l)
        {
            Scalar a;
            a = (this->transProblem.materialLaw().fractionalW(l*1.0/K, globalPos, *eIt, localPos)-this->transProblem.materialLaw().fractionalW((l-1)*1.0/K, globalPos, *eIt, localPos))*K;
            if(a>maxa)
            maxa = a;
        }

        init(dt,maxa,eIt);

        Scalar fv = this->transProblem.variables.vTotal(*eIt,0)[0];
        fv /= this->transProblem.soil().porosity(globalPos, *eIt, localPos);

        fronttracking(dt, fv, globalPos,eIt, localPos);
        rot(updateVec,dxtt);
        //blau nur für streamlines möglich
    }
    else

    for (ElementIterator eIt = gridView.template begin<0>(); eIt != eendit; ++eIt)
    {
        // cell geometry type
        Dune::GeometryType gt = eIt->geometry().type();
        // cell center in reference element
        const LocalPosition
            &localPos = Dune::ReferenceElements<Scalar,dim>::general(gt).position(0, 0);

        //
        const GlobalPosition& globalPos = eIt->geometry().global(localPos);

        Scalar maxa = 0;
        for(int l=1;l<K;++l)
        {
            Scalar a;
            a = (this->transProblem.materialLaw().fractionalW(l*1.0/K, globalPos, *eIt, localPos)-this->transProblem.materialLaw().fractionalW((l-1)*1.0/K, globalPos, *eIt, localPos))*K;
            if(a>maxa)
            maxa = a;
        }

        calcStreamline(dt,1.4*maxa,eIt);

        if(ch.size()>0)
        {
            //  	  typedef typename std::template list<slNode>::iterator SlistIterator;
            //    	  for(SlistIterator slit=sl.begin();slit!=sl.end();++slit)
            //           {
            //     	std::cout<<"index:"<<slit->index<<std::endl;
            //    	std::cout<<"xi:"<<slit->xi<<std::endl;
            //     	std::cout<<"Sat:"<<this->transProblem.variables.saturation[slit->index]<<std::endl;
            //           }
            fronttracking(dt,1.0, globalPos, eIt, localPos);
            //  	  for(ListIterator klit = ch.begin(); klit != ch.end();++klit)
            //            {
            //      	std::cout<<"Sättigung: "<<klit->Sat<<std::endl;
            //      	std::cout<<"Position: "<<klit->x[1]<<std::endl;
            //      	std::cout<<"Steigung: "<<klit->m<<std::endl;
            //      	std::cout<<std::endl;
            //            }
            rot(updateVec,dxtt);
            // blau(updateVec,indexi);
        }
    }

    for (ElementIterator eIt = gridView.template begin<0>(); eIt != eendit; ++eIt)
    {
        int indexi = elementmapper.map(*eIt);
        if(dxtt[indexi]>0)
        {
            updateVec[indexi] /= dxtt[indexi];
            updateVec[indexi] -= this->transProblem.variables.saturation[indexi];
        }
    }
    updateVec /= dt;

    dt/=cFLFactor;

    return 0;

}

template<class Grid, class Scalar, class VC, class Problem>
void ChTransport<Grid, Scalar, VC, Problem>::initialTransport()
{
    //	std::cout<<"initsat = "<<&this->transProblem.variables.saturation<<std::endl;
    // iterate through leaf grid an evaluate c0 at cell center

#define IMP
#ifdef IMP

    ElementIterator eendit = this->grid_.template lend<0>(this->level());
    for (ElementIterator eIt = this->grid_.template lbegin<0>(this->level()); eIt != eendit; ++eIt)
    {
        // get geometry type
        Dune::GeometryType gt = eIt->geometry().type();

        // get cell center in reference element
        const Dune::FieldVector<ct,dim>
        &local = Dune::ReferenceElements<ct,dim>::general(gt).position(0, 0);

        // get global coordinate of cell center
        Dune::FieldVector<ct,dimWorld> global = eIt->geometry().global(local);

        // initialize cell concentration
        this->transProblem.variables.saturation[elementmapper.map(*eIt)] = this->transProblem.initSat(global, *eIt, local);
    }
#else
    char impFileName[120];
    strcpy(impFileName, "data2.dgf");

    FILE* dataFile;
    dataFile=fopen(impFileName, "r");
    //      dataFile.open(impFileName.c_str());
    for(int i=0;i<this->transProblem.variables.saturation.size();++i)
    fscanf(dataFile,"%lf",&this->transProblem.variables.saturation[i]);
    fclose(dataFile);
#endif
    return;
}

}
#endif
