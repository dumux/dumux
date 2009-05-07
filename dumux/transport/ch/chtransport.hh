// $Id$

#ifndef DUNE_CHTRANSPORT_HH
#define DUNE_CHTRANSPORT_HH

#include <dune/common/helpertemplates.hh>
#include <dune/common/typetraits.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/bvector.hh>
#include "dumux/transport/transport.hh"
#include "dumux/transport/transportproblem.hh"
#include "dumux/transport/ch/fluxfunction.hh"

/**
 * @file
 * @brief  Base class for defining an instance of a numerical transport model
 * @brief  Forward characteristics method
 * @author Annika Fuchs; last changed by Yufei Cao
 * \defgroup transport Transport
 */

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
     *  \f[ \sum_{j \in \mathcal{N}(i)} v_{ij}grid(S_i, S_j) - v_{ji}grid(S_j, S_i), \f]
     *  where \f$\mathcal{N}(i)\f$ denotes the set of neighbors of the cell \f$i\f$ and
     *  \f$grid\f$ stands for \a numericalFlux. The normal velocities \f$v_{ij}\f$ and \f$v_{ji}\f$
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
     * @param grid a DUNE grid object
     * @param problem an object of class TransportProblem or derived
     * @param fluxFunc an object of flux function
     */
    ChTransport(Grid& grid, Problem& problem, FluxFunction<Grid,Scalar>& fluxFunc = *(new FluxFunction<Grid,Scalar>), int K=1000) :
    Transport<Grid, Scalar, VC, Problem>(grid, problem),
    elementmapper(grid.levelView(this->level())), fluxFunc_(fluxFunc), K(K)
    {}

private:

    int calcm(GlobalPosition globalPos, ElementIterator eIt, LocalPosition localPos);
    int init(const Scalar dt,const Scalar lambda, EntityPointer& eIt);
    int calcStreamline(const Scalar dt,const Scalar lambda, EntityPointer& eIt);
    int solveRP(Scalar& mch, std::list<ChNode>& chrar,Scalar sat0,Scalar sat1, GlobalPosition globalPos, ElementIterator eIt, LocalPosition localPos);
    Scalar linearflux(Scalar sat, GlobalPosition globalPos, ElementIterator eIt, LocalPosition localPos);
    int fronttracking(const Scalar& dt, Scalar v, GlobalPosition globalPos, ElementIterator eIt, LocalPosition localPos);
    int approxSol(RepresentationType& updateVec,RepresentationType& dxtt);
private:
    EM elementmapper;
    const FluxFunction<Grid, Scalar>& fluxFunc_;
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
    double fa = fluxFunc_(satA, globalPos, *eIt, localPos);
    double fb = fluxFunc_(satA+1.0/K, globalPos, *eIt, localPos);

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
        Scalar satil = fluxFunc_(il*1.0/K, globalPos, *eIt, localPos);
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
            fb[i]= fluxFunc_(b[i], globalPos, *eIt, localPos);
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
	//lower value on the left side:
	//build concave envelope
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

	//it's a rarefaction wave
	// m contains the maximal and minimal propagation speed
        m[low]=(fb[j1]-fsl)/(b[j1]-sl);
        m[high]=(fb[j2]-fsh)/(b[j2]-sh);

        ChNode c0;
        c0.Sat=b[j1];
        c0.m[high]=m[high];
        c0.m[low]=m[low];
        chrar.push_front(c0);

	//build new nodes with addidional shocks
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
    //there is one new node but the propagation speed left and right are the same.
    //it's a numerical error. it's a shock.
    if(chrar.size()==1 && fabs(chrar.front().m[0]-chrar.back().m[1])<1e-5)
    chrar.clear();
    return 0;
}


//calcm() calculates the propagation speed of discontinuities in the dynamic grid
template<class Grid,class Scalar,class VC, class Problem>
int ChTransport<Grid,Scalar,VC, Problem>::calcm(GlobalPosition globalPos, ElementIterator eIt, LocalPosition localPos)
{

  
    ListIterator lit=ch.begin();
    ListIterator lit2(lit);
    lit2++;
    //slope over the elements of the dynamic grid
    while(lit2!=ch.end())
    {
        std::list<ChNode> chrar;
        Scalar m;

	//solving the Riemannproblem
        solveRP(m,chrar,lit->Sat,lit2->Sat,globalPos, eIt, localPos);

        if(chrar.empty())
        {
	  //shock wave
	  if(m!=0)
          {
	    lit->m[1]=m;
	    lit2->m[0]=m;
	  }
	  //saturations left and right are idential.
	  //unifying the two elements 
	  else
          {
                lit2->m[0]=lit->m[0];
                lit2->x[0][0]=lit->x[0][0];
                lit=ch.erase(lit);	      
	  }
        }
	//rarefaction wave 
	//new nodes for addidional shocks
        else
        {
            lit->m[1]=chrar.front().m[0];
            lit2->m[0]=chrar.back().m[1];
            for(ListIterator kit=chrar.begin();kit!=chrar.end();++kit)
            for(int q=0;q<2;++q)
            kit->x[0][q]=lit2->x[0][0];

            ch.insert(lit2,chrar.begin(),chrar.end());
        }

        lit=lit2;
        lit2++;
    }
    //in boundary elements set the propagation speed of boundary fronts so
    //that discontinuity curves run in a parallel way
    ch.front().m[0] = ch.front().m[1];
    ch.back().m[1]=ch.back().m[0];
    
return 0;
}


//for one dimensional problems you don't need streamlines. 
//init() build the dynamic grid, which contains the exact solution of
//the discretized problem.
template<class Grid,class Scalar,class VC, class Problem>
int ChTransport<Grid,Scalar,VC, Problem>::init(const Scalar dt,const Scalar lambda, EntityPointer& eIt)
{
    ch.clear();
    sl.clear();

    int indexi = elementmapper.map(*eIt);

    const GV& gridView = this->grid_.levelView(this->level());

    //build the basic element of the dynamic grid
    ChNode i;
    i.Sat=this->transProblem.variables.saturation[indexi];

    //set position and propagation speed of the basic element
    IntersectionIterator endis = gridView.template iend(*eIt);
    for(IntersectionIterator is=gridView.template ibegin(*eIt);is!=endis;++is)
    {
        GeometryType gtf = is->geometryInInside().type();
        const Dune::FieldVector<ct,dim-1>& facelocal = Dune::ReferenceElements<Scalar,dim-1>::general(gtf).position(0,0);
        const VType& faceglobal = is->geometry().global(facelocal);
        int numberInSelf = is->indexInInside();

        if(numberInSelf == 0)
        i.x[0][0] = faceglobal[0];
        else
        i.x[0][1] = faceglobal[0];
    }
    ch.push_front(i);

    //use the streamline list for this problem to use the same approximation function 
    slNode s;
    s.index=indexi;
    s.xi=ch.front().x[0];
    sl.push_front(s);

    //running through the elements of the basic grid in one direction by hitting a boundary point
    //then the same procedure in the other direction
    //build new elements for the dynamic grid if there's a discontinuity jump
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
                int numberInSelf = is->indexInInside();

                s.xi[numberInSelf]=faceglobal[0];

		//p==0 -> running through in left direction
		//p==1 -> running through in right direction
                if (numberInSelf == p) 
                {
                    j.x[0][0] = faceglobal[0];
                    j.x[0][1] = j.x[0][0];
                 
		    //there's a neighbor element
                    if(is->neighbor())
                    {
                        ep=is->outside();
                        indexi=elementmapper.map(*ep);
                        j.Sat=this->transProblem.variables.saturation[indexi];
                    }
		    //it's the boundary point. Break off the running in this direction
                    else
                    {
                        
                        const Dune::FieldVector<ct,dim>& facelocalDim = Dune::ReferenceElements<ct,dim>::general(gtf).position(numberInSelf,1);
                        j.Sat =this->transProblem.dirichletSat(faceglobal,*ep,facelocalDim);
                        abb=-1;
                    }

		    //if there's a saturation jump build a new node for the dynamic grid
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

    //set the boundary points of the dynamic grid far enough to cover the static grid all the time
    ch.front().x[0][0] -= lambda*dt;
    ch.back().x[0][1] += lambda*dt;


    if(ch.size()<2)
    {
        ch.clear();
        sl.clear();
    }

    return 0;

}


//in more dimensions you need a streamline method to translate the more dimensional problem in a one dimensional.
//calcStreamline() approximates the streamlines in the velocity field by backtracking particles.
//in each node one particle is positioned in element midpoint at point of time dt.
template<class Grid,class Scalar,class VC, class Problem>
int ChTransport<Grid,Scalar,VC, Problem>::calcStreamline(const Scalar dt,const Scalar lambda, EntityPointer& eIt)
{
    const GV& gridView = this->grid_.levelView(this->level());

    ch.clear();
    sl.clear();

    EntityPointer ep(eIt);
    Scalar xi = 0;

    Dune::GeometryType gt = eIt->geometry().type();
    const LocalPosition
        &localPos = Dune::ReferenceElements<Scalar,dim>::general(gt).position(0, 0);
    const GlobalPosition& globalPos = eIt->geometry().global(localPos);

    int indexi = elementmapper.map(*eIt);
    const Dune::FieldVector<ct,dim>& local = Dune::ReferenceElements<ct,dim>::general(eIt->geometry().type()).position(0,0);
    Dune::FieldVector<ct,dim> global = eIt->geometry().global(local);
    VType xq(global);

    //build basic nodes of the dynamic grid and the streamline information list
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

    //the particle is backtracked till -lambda*dt. so it's sure, that each element has a solution at time dt.
    //lambda is the minimum time multiplicated with a safety factor.
    while(xi> -lambda*dt)
    {
        Dune::FieldVector<VType,2*dim> xx;
        Dune::FieldVector<VType,2*dim> vx;

        //set positions and velocities of the edges of the current element
        IntersectionIterator endis = gridView.template iend(*ep);
        for(IntersectionIterator is=gridView.template ibegin(*ep);is!=endis;++is)
        {
            GeometryType gtf = is->geometryInInside().type();
            const Dune::FieldVector<ct,dim-1>& facelocal = Dune::ReferenceElements<Scalar,dim-1>::general(gtf).position(0,0);
            const VType& faceglobal = is->geometry().global(facelocal);
            int numberInSelf = is->indexInInside();

            xx[numberInSelf] = faceglobal;

            vx[numberInSelf] = this->transProblem.variables.vTotal(*ep,numberInSelf);
            vx[numberInSelf] /= this->transProblem.soil().porosity(globalPos, *eIt, localPos);
        }
	//calculation of the gradient of velocity and the velocity of the particle
        VType Ax(0);
        VType vxq(0);

        for(int i=0;i<dim;++i)
        {
            Ax[i] = (vx[2*i+1][i]-vx[2*i][i])/(xx[2*i+1][i]-xx[2*i][i]);
            vxq[i] = Ax[i]*(xq[i]-xx[2*i][i])+vx[2*i][i];

        }

	//calculation of the direction the particle comes from. If vxq<0, the particle comes from the right or upper neighbor element,
	//otherwise from the left or lower neighbor element
        //kant contains the index of the accordant edges
        //dtx contains the time of flight, from the accordant edge to the current position.
        VType dtx(-2*lambda*dt);
        Dune::FieldVector<int,dim> kant(-1);

        for(int i = 0;i<dim;++i)
        {
	  //constant velocity in i-th direction
            if(fabs(Ax[i])<1e-5)
            {
                if(vxq[i]<0)
		  dtx[i] = (xx[2*i+1][i]-xq[i])/vxq[i];
                else if(vxq[i]>0)
		  dtx[i] = (xx[2*i][i]-xq[i])/vxq[i];
            }
	    //gradient unequal zero in i-th direction
            else
	      dtx[i] = 1./Ax[i]*log(vx[i*2+1][i]/vxq[i]);

	    //set the accordant edge index
            if(vxq[i]<0)
	      kant[i] = 2*i+1;
            else if(vxq[i]>0)
	      kant[i] = 2*i;
        }

        //calculation of time of flight
	//m contains the edge index of the actual entry edge
	//dte contains the actual time of flight
        int m = 0; 
        Scalar dte = dtx[m]; 
        for(int i=1;i<dim;++i)
        if(dtx[i]>dte)
        {
            m = i;
            dte=dtx[i];
        }

        //calculating the position of the entry point
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
	//in one dimension position of entry point and entry edge are identical
        else
	  xq = xx[kant[m]];

        //calcuting the entry time 
        xi = xi+dte;
	//building new nodes for the dynamic grid and the streamline information list
        ch.front().x[0][0]=xi;
        ChNode j;
        j.x[0][0]=xi;
        j.x[0][1] = xi;
        sl.front().xi[0]=xi;
        slNode s;
        s.xi=j.x[0];
        s.index = -1;

	//searching saturation jumps
        for(IntersectionIterator is = gridView.template ibegin(*ep);is!=endis;++is)
        {
            int indexInInside = is->indexInInside();

            if(indexInInside == kant[m])
	      if(is->neighbor())
	      {
                ep = is->outside();
		//the entry is not a edge of the element but a element point
		//the neighbor element is the diagonal element
                if(n > 0)
                {
                    for(int i=0;i<dim;++i)
		      if(i!=m && fabs(dtx[i]-dte)<1e-7)
		      {
                        IntersectionIterator endir = gridView.template iend(*ep);
                        for(IntersectionIterator ir = gridView.template ibegin(*ep);ir!=endir;++ir)
                        {
			  //setting informations of the diagonal element
                            int nIS = ir->indexInInside();

                            if(nIS == kant[i])
			      if(ir->neighbor())
			      {
                                ep = ir->outside();
                                indexi = elementmapper.map(*ep);
                                j.Sat = this->transProblem.variables.saturation[indexi];

                                s.index=indexi;
			      }
			    //boundary conditions
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

		//neighbor element is the entry element
                else
                {
                    indexi = elementmapper.map(*ep);
                    j.Sat = this->transProblem.variables.saturation[indexi];

                    s.index=indexi;
                }
            }//if is->neighbor()
	    
	    //particle comes from boundary
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

    //boundary front of the streamline section
    ch.front().x[0][0] = xi;
    sl.front().xi[0] = xi;

    if(ch.size()<2)
    {
        ch.clear();
        sl.clear();
    }

    return 0;
}


//fronttracking() tracks the discontinuities and calculates new shock fronts if there are crossings.
//at the end at point of time dt there is an exact solution 

template<class Grid,class Scalar,class VC, class Problem>
int ChTransport<Grid,Scalar,VC, Problem>::fronttracking(const Scalar& dt, Scalar v, GlobalPosition globalPos, ElementIterator eIt, LocalPosition localPos)
{
    //calculating the propagation speeds of the discontinuities
    calcm(globalPos, eIt, localPos);

    //position of each shock front at point of time dt
    ListIterator lit;
    for(lit=ch.begin();lit!=ch.end();++lit)
    for(int q=0;q<2;++q)
    lit->x[1][q]=lit->x[0][q]+lit->m[q]*dt*v;

    //t1 is the current point of time where an exact solution is given.
    Scalar t1=0;
    while(t1<dt)
    {
        Scalar dtau=dt;
        ListIterator merk;
        //if the left front of an element is at the right side of the right
	//front at dtau, there's a crossing in [t1,t1+dtau]
	//calculation of the minimal time of crossings 
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
	
	//no further crossings in [t1,dtau]
        if(t1+dtau>=dt)
        dtau=dt-t1;

        //calculation of the position at point of time t1+dtau
        //xt=x+m*dtau
        //set the next point of time an exact solution will be given: t1 = t1+dtau
        for(lit=ch.begin();lit!=ch.end();++lit)
        for(int p=0;p<2;++p)
        {
            lit->x[1][p]=lit->x[0][p]+lit->m[p]*dtau*v;
            lit->x[0][p]=lit->x[1][p];
        }
        t1=t1+dtau;

	//no further crossings. at point of time dt is an exact solution given
        if(t1>=dt)
	  return 0;

	//calculating new front / fronts at the crossing point
        ListIterator chend = ch.end();
        lit=merk;
	//it's not a boundary element
        if(merk!=ch.begin() && merk!=--chend)
        {
            ListIterator merk2(lit);
            merk--; 
            merk2++; 

            std::list<ChNode> chrar;
            Scalar m;

	    //solving the Riemannproblem with initial data merk.sat, merk2.sat
            solveRP(m,chrar,merk->Sat,merk2->Sat,globalPos, eIt, localPos);

            //set propagation speed
            if(chrar.empty())
            {
	      //shock, deleting the element with crossing fronts 
                if(m!=-0)
                {
                    merk->m[1]=m;
                    merk2->m[0]=m;
                    ch.erase(lit);

                    merk2->x[0][0]=0.5*(merk->x[0][1]+merk2->x[0][0]);
                    merk->x[0][1] = merk2->x[0][0];
                }
		//saturation is identical, unifying the elements
                else
		{
		  merk2->m[0]=merk->m[0];
		  merk2->x[0][0]=merk->x[0][0];
		  
		  ch.erase(lit);
                }
            }
	    //rarefaction wave. deleting the element with crossing fronts
	    //inserting the additional elements for approximation
            else
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
	//deleting the boundary element with crossing fronts
        else
	  ch.erase(lit);

        //calculation the shock front positions at point of time dt
        for(lit=ch.begin();lit!=ch.end();++lit)
	  for(int q=0;q<2;++q)
	    lit->x[1][q]=lit->x[0][q]+lit->m[q]*(dt-t1)*v;
    }

    return 0;
}


//approxSol() approximates the solution of the dynamic grid as piecewise constant at the static grid.
//The function works with integral conservation.
template<class Grid, class Scalar, class VC, class Problem>
int ChTransport<Grid,Scalar,VC, Problem>::approxSol(RepresentationType& updateVec,RepresentationType& dxtt)
{

    typedef typename std::template list<slNode>::iterator SlistIterator;
    for(SlistIterator slit=sl.begin();slit!=sl.end();++slit)
    {
        Scalar x0=slit->xi[0];
        Scalar x1=slit->xi[1];

        Scalar dx=fabs(x1-x0);
        Scalar sat=0;
        Scalar dxt=0;

	//integral conservation over the different grids. 
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

	//sat contains the saturation
	//dxtt contains the width of the streamline section in the element
        if(dxt>0)
        {
            updateVec[slit->index]+=sat;
            dxtt[slit->index]+=dxt;
        }
    }
    return 0;
}

template<class Grid, class Scalar, class VC, class Problem>
int ChTransport<Grid,Scalar,VC, Problem>::update(const Scalar t, Scalar& dt, RepresentationType& updateVec, Scalar& cFLFac = 1)
{
    const GV& gridView = this->grid_.levelView(this->level());

    // the timestep is given by parameter 'dt'
    //if the problem allow it, use only one timestep!!!

    // set update vector to zero
    updateVec = 0;
    RepresentationType dxtt(updateVec);
    dxtt=0;

    //one dimensional case
    if(Grid::dimension == 1)
    {

        ElementIterator eIt = gridView.template begin<0>();
        Dune::GeometryType gt = eIt->geometry().type();
        const LocalPosition
            &localPos = Dune::ReferenceElements<Scalar,dim>::general(gt).position(0, 0);
	const GlobalPosition& globalPos = eIt->geometry().global(localPos);

	//calculation of the maximal propagation velocity of a shock front
        Scalar maxa = 0;
        for(int l=1;l<K;++l)
        {
            Scalar a;
            a = (fluxFunc_(l*1.0/K, globalPos, *eIt, localPos)-fluxFunc_((l-1)*1.0/K, globalPos, *eIt, localPos))*K;
            if(a>maxa)
            maxa = a;
        }

	//build dynamic grid
        init(dt,maxa,eIt);

	//constant velocity in the velocity field
        Scalar fv = this->transProblem.variables.vTotal(*eIt,0)[0];
        fv /= this->transProblem.soil().porosity(globalPos, *eIt, localPos);

	//tracking the discontinuities to get an exact solution on the dynamic grid
        fronttracking(dt, fv, globalPos,eIt, localPos);
	//approximating the solution on the dynamic grid as piecewise constant 
	//on the static grid
        approxSol(updateVec,dxtt);
        
    }
    //more dimensional case
    else
    {
      //calculation streamlines for each element
      ElementIterator eendit = gridView.template end<0>();
      for (ElementIterator eIt = gridView.template begin<0>(); eIt != eendit; ++eIt)
      {
        Dune::GeometryType gt = eIt->geometry().type();
        const LocalPosition
            &localPos = Dune::ReferenceElements<Scalar,dim>::general(gt).position(0, 0);
        const GlobalPosition& globalPos = eIt->geometry().global(localPos);

	//calculating the maximal propagation speed of a shock front
        Scalar maxa = 0;
        for(int l=1;l<K;++l)
        {
            Scalar a;
            a = (fluxFunc_(l*1.0/K, globalPos, *eIt, localPos)-fluxFunc_((l-1)*1.0/K, globalPos, *eIt, localPos))*K;
            if(a>maxa)
            maxa = a;
        }
        calcStreamline(dt,1.4*maxa,eIt);

	//if there are saturation jumps, solving the exact solutions
	// and approximating it on the static grid
        if(ch.size()>0)
        {
            fronttracking(dt,1.0, globalPos, eIt, localPos);
            approxSol(updateVec,dxtt);
       
        }
      }
    }
    //the different solutions of the streamline are weightend over the length of the
    //accordant streamline sections
    //updateVec = (Sat_new-Sat_old)/dt
    ElementIterator eendit = gridView.template end<0>();
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

    dt/=cFLFac;

    return 0;

}

template<class Grid, class Scalar, class VC, class Problem>
void ChTransport<Grid, Scalar, VC, Problem>::initialTransport()
{

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
    return;
}

}
#endif
