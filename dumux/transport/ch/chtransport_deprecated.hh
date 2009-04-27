// $Id$

#ifndef DUNE_CHTRANSPORT_HH
#define DUNE_CHTRANSPORT_HH

#include <dune/common/helpertemplates.hh>
#include <dune/common/typetraits.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/istl/bvector.hh>
#include "dumux/transport/transport_deprecated.hh"
#include "dumux/transport/fv/numericalflux.hh"
#include "dumux/transport/fv/diffusivepart.hh"

/**
 * @file
 * @brief  Base class for defining an instance of a numerical transport model
 * @brief  Forward characteristics method
 * @author Annika Fuchs
 * \defgroup transport Transport
 */

namespace Dune {

  //! \ingroup transport
  //! The finite volume model for the solution of the transport equation
  template<class G, class RT, class VC> class ChTransport : 
    public Transport< G, RT, VC> {


    template<int dim> struct ElementLayout {
      bool contains(Dune::GeometryType gt) {
	return gt.dim() == dim;
      }
    };


    //structure of an element in the dynamic grid
    //Sat:    saturation
    //x[i,j]: position of the i-front in time t_j. i,j \in {0,1}
    //m[i]:   propagation velocity of the i-front. i \in {0,1}
    struct ChNode
    {
      RT Sat;
      Dune::FieldVector<Dune::FieldVector<RT,2>,2 > x;
      Dune::FieldVector<RT,2> m;
    };

    //node of the streamline-list
    //i:     index of the element the streamline crossed
    //xi[i]: streamline parameter of the entrance and exit points 
    struct slNode
    {
      int index;
      Dune::FieldVector<RT,2> xi;
    };


  enum {dim = G::dimension};
  enum {dimworld = G::dimensionworld};
  typedef BlockVector< Dune::FieldVector<RT,1> > PressType;
  typedef BlockVector< FieldVector<FieldVector<RT, G::dimension>, 2*G::dimension> >
			VelType;
  typedef typename G::Traits::template Codim<0>::Entity Entity;
  typedef typename G::LevelGridView GV;
  typedef typename GV::IndexSet IS;
  typedef typename GV::template Codim<0>::Iterator Iterator;
  typedef Dune::MultipleCodimMultipleGeomTypeMapper<GV,ElementLayout> EM;
  typedef typename G::template Codim<0>::EntityPointer EntityPointer;
    typedef typename G::LevelGridView::IntersectionIterator IntersectionIterator;
    typedef typename std::template list<ChNode>::iterator ListIterator;
  typedef typename G::ctype ct;
 typedef Dune::FieldVector<RT,dim> VType;
  typedef BlockVector< Dune::FieldVector<RT,dim> > SlopeType;

public:
    typedef BlockVector< Dune::FieldVector<RT,1> > RepresentationType;
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
    int update(const RT t, RT& dt, RepresentationType& updateVec, RT& cFLFac);

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
    ChTransport(G& g, TransportProblem<G, RT, VC>& prob, int lev = 0,
			DiffusivePart<G,RT>& diffPart = *(new DiffusivePart<G, RT>), bool rec = false,
		    double amax = 0.8, const NumericalFlux<RT>& numFl = *(new Upwind<RT>), int K=1000) :
	  Transport<G, RT, VC>(g, prob, lev),
				elementmapper( g.levelView(lev)), gridview(g.levelView(lev)),
				indexset(gridview.indexSet()), reconstruct(rec),
	  numFlux(numFl), diffusivePart(diffPart), alphamax(amax),K(K),grid(g)
				{}

private:

    int calcm();
    int init(const RT dt,const RT lambda, EntityPointer& it);
    int calcStreamline(const RT dt,const RT lambda, EntityPointer& it);
    int solveRP(RT& mch, std::list<ChNode>& chrar,RT sat0,RT sat1);
    RT linearflux(RT sat);
    int fronttracking(const RT dt, const RT v);
    int approxSol(RepresentationType& updateVec,RepresentationType& dxtt);
private:
	EM elementmapper;
	const GV& gridview;
	const IS& indexset;
	bool reconstruct;
	const NumericalFlux<RT>& numFlux;
	const DiffusivePart<G, RT>& diffusivePart;
	double alphamax;	
	G& grid;

	int K;
    	std::list<ChNode> ch;
	std::list<slNode> sl;
};




 //approximation of the fractional flow in a piecewise linear function
 //evaluation at the parameter RT sat
 template<class G,class RT,class VC> 
 RT ChTransport<G,RT,VC>::linearflux(RT sat)
  {
    RT satA=-1;
    //search of satA with sat \in [satA, satA+i/K]
    for(int i=0;i<K;++i)
      if(K*sat-i>=0)
	satA=i*1.0/K;

    //calculation of f(satA), f(satA+i/K)
    double fa=this->transproblem.materialLaw.fractionalW(satA);
    double fb=this->transproblem.materialLaw.fractionalW(satA+1.0/K);
    
    //evaluation at sat
    double fsat=fa+(sat-satA)*K*(fb-fa);
    return(fsat);
  }

  //Riemann problem solver
  //mch:   propagation velocity of the front, if the solution is a shock
  //chrar: list of dynamic grid nodes, if the solution is a rarefaction wave
  //sat0:  value on the left side of the discontinuity
  //sat1:  value on the right side of the discontinuity
  template<class G,class RT,class VC>
  int ChTransport<G,RT,VC>::solveRP(RT& mch, std::list<ChNode>& chrar,RT sat0,RT sat1)
  {
    //clear temporary dynamic grid list
    chrar.clear();
  
    std::vector<RT> sat(2,0);
    std::vector<RT> m(2,0);
    sat[0]=sat0;
    sat[1]=sat1;

    //left and right value are identical
    //function returns the propagation speed of the characteristic
    if(fabs(sat1-sat0)<1e-5)   
    {
      mch=0;
      return (linearflux(sat0));
    }
    
    //Rankine Hugoniot condition
    mch=(linearflux(sat1)-linearflux(sat0))/(sat1-sat0);
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

    RT sl=sat[low];
    RT sh = sat[high];
    double fsl=linearflux(sl);
    double fsh = linearflux(sh);
    
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
    if(il==0 || il > ih)
      return 0;
    
    //one node between sat0 and sat1
    if(il==ih)
    {
      RT satil = this->transproblem.materialLaw.fractionalW(il*1.0/K);
      //Rankine Hugoniot condition: shockfront	
      if(sat0 > sat1 && (fsl+mch*(il*1.0/K-sl) >= satil))
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
	fb[i]=this->transproblem.materialLaw.fractionalW(b[i]);
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
      // m contains the maximal and minimal propagation speedy
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
  }
      
  


  //calcm() calculates the propagation speed of discontinuities in the dynamic grid
  template<class G,class RT,class VC>
  int ChTransport<G,RT,VC>::calcm()
  {

    ListIterator lit=ch.begin();
    ListIterator lit2(lit);
    lit2++;
    //slope over the elements of the dynamic grid
    while(lit2!=ch.end())
    {
      std::list<ChNode> chrar;
      RT m;
      
      //solving the Riemannproblem
      solveRP(m,chrar,lit->Sat,lit2->Sat);
     
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
  template<class G,class RT,class VC>
  int ChTransport<G,RT,VC>::init(const RT dt,const RT lambda, EntityPointer& it)
  {
    ch.clear();
    sl.clear();
   
    int indexi = elementmapper.map(*it);
    
    //build the basic element of the dynamic grid
    ChNode i;
    i.Sat=this->transproblem.variables.saturation[indexi];
  
    //set position and propagation speed of the basic element
    IntersectionIterator endis=it->ilevelend();
    for(IntersectionIterator is=it->ilevelbegin();is!=endis;++is)
    {
      GeometryType gtf = is->geometryInInside().type();
      const Dune::FieldVector<ct,dim-1>& facelocal = Dune::ReferenceElements<RT,dim-1>::general(gtf).position(0,0);
      const VType& faceglobal = is->geometry().global(facelocal);
      int numberInSelf = is->indexInInside();

      //linke Kante
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
      EntityPointer ep=it;
      RT satI=i.Sat;

  
      while(abb>0)
      {
	ChNode j;
	slNode s;
	s.index=indexi;
	
	IntersectionIterator endis = ep->ilevelend();
	for(IntersectionIterator is = ep->ilevelbegin();is!=endis;++is)
        {

	  GeometryType gtf = is->geometryInInside().type();
	  const Dune::FieldVector<ct,dim-1>& facelocal = Dune::ReferenceElements<RT,dim-1>::general(gtf).position(0,0);
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
	      j.Sat=this->transproblem.variables.saturation[indexi]; 
	    }
	    //it's the boundary point. Break off the running in this direction
	    else
	    {
	      const Dune::FieldVector<ct,dim>& facelocalDim = Dune::ReferenceElements<ct,dim>::general(gtf).position(numberInSelf,1);
	      j.Sat =this->transproblem.dirichlet(faceglobal,*ep,facelocalDim);
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
    else
      calcm();

    return 0;

  }



  //in more dimensions you need a streamline method to translate the more dimensional problem in a one dimensional.
  //calcStreamline() approximates the streamlines in the velocity field by backtracking particles.
  //in each node one particle is positioned in element midpoint at point of time dt.
  template<class G,class RT,class VC>
  int ChTransport<G,RT,VC>::calcStreamline(const RT dt,const RT lambda, EntityPointer& it)
  {
    ch.clear();
    sl.clear();

    EntityPointer ep(it);
    RT xi = 0;
  
    int indexi = elementmapper.map(*it);
    const Dune::FieldVector<ct,dim>& local = Dune::ReferenceElements<ct,dim>::general(it->geometry().type()).position(0,0);
    Dune::FieldVector<ct,dim> global = it->geometry().global(local);
    VType xq(global);

    //build basic nodes of the dynamic grid and the streamline information list
    ChNode i;   
    i.Sat=this->transproblem.variables.saturation[indexi];
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
    while(xi > -lambda*dt)
    {    
      Dune::FieldVector<VType,2*dim> xx;
      Dune::FieldVector<VType,2*dim> vx;
     
      //set positions and velocities of the edges of the current element
      IntersectionIterator endis=ep->ilevelend();
      for(IntersectionIterator is=ep->ilevelbegin();is!=endis;++is)
      {
	GeometryType gtf = is->geometryInInside().type();
	const Dune::FieldVector<ct,dim-1>& facelocal = Dune::ReferenceElements<RT,dim-1>::general(gtf).position(0,0);
	const VType& faceglobal = is->geometry().global(facelocal);

	int numberInSelf = is->indexInInside();
	xx[numberInSelf] = faceglobal;
      
	vx[numberInSelf] = this->transproblem.variables.vTotal(*ep,numberInSelf);
       	vx[numberInSelf] /= this->transproblem.porosity();
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
	if(fabs(Ax[i])<1e-5)
	  //constant velocity in i-th direction
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
      RT dte = dtx[m];            
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
      for(IntersectionIterator is = ep->ilevelbegin();is!=endis;++is)
      {
	int numberInSelf = is->indexInInside();
	if(numberInSelf == kant[m])
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
		  IntersectionIterator endir = ep->ilevelend();
		  for(IntersectionIterator ir = ep->ilevelbegin();ir!=endir;++ir)
		  {
		    //setting informations of the diagonal element
		    int nIS = ir->indexInInside();
		    
		    if(nIS == kant[i])
		      if(ir->neighbor())
		      {
			ep = ir->outside();
			indexi = elementmapper.map(*ep);
			j.Sat = this->transproblem.variables.saturation[indexi];		      
		   
			s.index=indexi;
		      }
		    //boundary conditions
		      else
		      {
			GeometryType gtf = ir->geometryInInside().type();
			const Dune::FieldVector<ct,dim-1>& facelocal = Dune::ReferenceElements<RT,dim-1>::general(gtf).position(0,0);
			const VType& faceglobal = ir->geometry().global(facelocal);
			const Dune::FieldVector<ct,dim>& facelocalDim = Dune::ReferenceElements<ct,dim>::general(gtf).position(numberInSelf,1);
			  
			j.Sat = this->transproblem.dirichlet(faceglobal,*ep,facelocalDim);
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
	      j.Sat = this->transproblem.variables.saturation[indexi];	      
	      
	      s.index=indexi;
	    }
	  }//if is->neighbor()

	//particle comes from boundary
	  else
	  {
	    GeometryType gtf = is->geometryInInside().type();
	    const Dune::FieldVector<ct,dim-1>& facelocal = Dune::ReferenceElements<RT,dim-1>::general(gtf).position(0,0);
	    const VType& faceglobal = is->geometry().global(facelocal);
	    const Dune::FieldVector<ct,dim>& facelocalDim = Dune::ReferenceElements<ct,dim>::general(gtf).position(numberInSelf,1);
	    j.Sat = this->transproblem.dirichlet(faceglobal,*ep,facelocalDim);

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
    else
      calcm();

    return 0;
  }

  //fronttracking() tracks the discontinuities and calculates new shock fronts if there are crossings.
  //at the end at point of time dt there is an exact solution 
  template<class G,class RT,class VC>
  int ChTransport<G,RT,VC>::fronttracking(const RT dt,const RT v = 1)
  {       
    //position of each shock front at point of time dt
    ListIterator lit;
    for(lit=ch.begin();lit!=ch.end();++lit)
      for(int q=0;q<2;++q)
	lit->x[1][q]=lit->x[0][q]+lit->m[q]*dt*v;


    //t1 is the current point of time where an exact solution is given.
    RT t1=0;
    while(t1<dt)       
    {   
      RT dtau=dt;
      ListIterator merk;
      //if the left front of an element is at the right side of the right
      //front at dtau, there's a crossing in [t1,t1+dtau]
      //calculation of the minimal time of crossings 
      for(lit=ch.begin();lit!=ch.end();++lit)
	if(lit->x[0][0]<=lit->x[0][1] && lit->x[1][0]>lit->x[1][1] )
	{
	  RT dch = lit->m[0]-lit->m[1];
	  RT tt=dt;
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
	merk--; //merk ist nun lit-1
	merk2++; //merk ist nun lit+1
     
	std::list<ChNode> chrar;
	RT m;
	
	//solving the Riemannproblem with initial data merk.sat, merk2.sat
	solveRP(m,chrar,merk->Sat,merk2->Sat);
	
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
  template<class G, class RT, class VC>
  int ChTransport<G,RT,VC>::approxSol(RepresentationType& updateVec,RepresentationType& dxtt)
  {

    typedef typename std::template list<slNode>::iterator  SlistIterator;
    for(SlistIterator slit=sl.begin();slit!=sl.end();++slit)
    {
      RT x0=slit->xi[0];
      RT x1=slit->xi[1];

      RT dx=fabs(x1-x0);
      RT sat=0;
      RT dxt=0;
  
      //integral conservation over the different grids. 
      ListIterator lit = ch.begin();
      while(dxt<dx && lit!=ch.end())
      { 
	RT xt0=lit->x[1][0];
	RT xt1=lit->x[1][1];

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


  template<class G, class RT, class VC>
  int ChTransport<G,RT,VC>::update(const RT t, RT& dt, RepresentationType& updateVec, RT& cFLFactor) 
  {
    //setting of the timestep 
    //if the problem allow it, use only one timestep!!!
    if(t==0)
      dt = 4.32e7;

    //calculation of the maximal propagation velocity of a shock front
    RT maxa = 0;
    for(int l=1;l<K;++l) 
    {
      RT a;
      a = (this->transproblem.materialLaw.fractionalW(l*1.0/K)-this->transproblem.materialLaw.fractionalW((l-1)*1.0/K))*K;
      if(a>maxa)
	maxa = a;
    }
       
    // set update vector to zero   
    updateVec = 0;    
    RepresentationType dxtt(updateVec);
    dxtt=0;

    //one dimensional case
    if(G::dimension == 1)
    {
      Iterator it = this->grid.template lbegin<0>(this->level());

      //build dynamic grid
      init(dt,maxa,it);

      //constant velocity in the velocity field
      RT fv = this->transproblem.variables.vTotal(*it,0)[0];
      fv /= this->transproblem.porosity();     

      //tracking the discontinuities to get an exact solution on the dynamic grid
      fronttracking(dt,fv);
      //approximating the solution on the dynamic grid as piecewise constant 
      //on the static grid
      approxSol(updateVec,dxtt);	             
    }
    //more dimensional case
    else 
    {
      //calculation streamlines for each element
      Iterator eendit = this->grid.template lend<0>(this->level());
      for (Iterator it = this->grid.template lbegin<0>(this->level()); it != eendit; ++it)
      {
	int indexi = elementmapper.map(*it);   
	calcStreamline(dt,1.4*maxa,it);

	//if there are saturation jumps, solving the exact solutions
	// and approximating it on the static grid
	if(ch.size()>0)
	{
	  fronttracking(dt);	 
	  approxSol(updateVec,dxtt);

	}
      }
    }
    //the different solutions of the streamline are weightend over the length of the
    //accordant streamline sections
    //updateVec = (Sat_new-Sat_old)/dt
    Iterator eendit = this->grid.template lend<0>(this->level());
    for (Iterator it = this->grid.template lbegin<0>(this->level()); it != eendit; ++it)
    {
      int indexi = elementmapper.map(*it);
      if(dxtt[indexi]>0)      
      {
	updateVec[indexi] /= dxtt[indexi];
	updateVec[indexi] -= this->transproblem.variables.saturation[indexi];
      }
    }
    updateVec /= dt;
      
    dt/=cFLFactor;
   
    return 0;
    
  } 

  template<class G, class RT, class VC> 
  void ChTransport<G, RT, VC>::initialTransport() {
    //	std::cout<<"initsat = "<<&this->transproblem.variables.saturation<<std::endl;
    // iterate through leaf grid an evaluate c0 at cell center

    Iterator eendit = this->grid.template lend<0>(this->level());
    for (Iterator it = this->grid.template lbegin<0>(this->level()); it != eendit; ++it) {
      // get geometry type
      Dune::GeometryType gt = it->geometry().type();
      
      // get cell center in reference element
      const Dune::FieldVector<ct,dim>
	&local = Dune::ReferenceElements<ct,dim>::general(gt).position(0, 0);
      
      // get global coordinate of cell center
      Dune::FieldVector<ct,dimworld> global = it->geometry().global(local);

      // initialize cell concentration
      this->transproblem.variables.saturation[elementmapper.map(*it)] = this->transproblem.initSat(global, *it, local);
    }

    return;
  }
  


}
#endif
