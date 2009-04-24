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
    int gelb(RepresentationType& updateVec,int indexi);   
    int blau(RepresentationType& updateVec,int indexi);   
    int rot(RepresentationType& updateVec,RepresentationType& dxtt);
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
  }
      
  


  template<class G,class RT,class VC>
  int ChTransport<G,RT,VC>::calcm()
  {

    //berechne Steigungen der Charakteristiken
    ListIterator lit=ch.begin();
    ListIterator lit2(lit);
    lit2++;
    while(lit2!=ch.end())
    {
      std::list<ChNode> chrar;
      RT m;
      
      solveRP(m,chrar,lit->Sat,lit2->Sat);
     
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
      ch.front().m[0] = solveRP(m,chrar,ch.front().Sat,ch.front().Sat);
      ch.back().m[1] = solveRP(m,chrar,ch.back().Sat,ch.back().Sat);
      //      ch.front().m[0] = ch.front().m[1];
      //ch.back().m[1]=ch.back().m[0];

      lit=lit2;
      lit2++;
    }
    return 0;
  }

  template<class G,class RT,class VC>
  int ChTransport<G,RT,VC>::init(const RT dt,const RT lambda, EntityPointer& it)
  {
    ch.clear();
    sl.clear();
   
    int indexi = elementmapper.map(*it);
    
    //Bilde Ausgangselement in der Charakteristikenliste
    ChNode i;
    i.Sat=this->transproblem.variables.saturation[indexi];
  
    //setze Geschwindigkeit und x-Koordinaten für das Ausgangselement
    IntersectionIterator endis=it->ilevelend();
    for(IntersectionIterator is=it->ilevelbegin();is!=endis;++is)
    {
      GeometryType gtf = is->geometryInInside().type();
      const Dune::FieldVector<ct,dim-1>& facelocal = Dune::ReferenceElements<RT,dim-1>::general(gtf).position(0,0);
      const VType& faceglobal = is->geometry().global(facelocal);
      //      const VType& faceInElem = it->geometry().local(faceglobal);
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
      EntityPointer ep=it;
      RT satI=i.Sat;

  
      while(abb>0)
      {
	ChNode j;
	slNode s;
	s.index=indexi;
	
	IntersectionIterator endis = IntersectionIteratorGetter<G,LevelTag>::end(*ep);
	for(IntersectionIterator is = IntersectionIteratorGetter<G,LevelTag>::begin(*ep);is!=endis;++is)
        {

	  GeometryType gtf = is->geometryInInside().type();
	  const Dune::FieldVector<ct,dim-1>& facelocal = Dune::ReferenceElements<RT,dim-1>::general(gtf).position(0,0);
	  const VType& faceglobal = is->geometry().global(facelocal);
	  //  const VType& faceInElem = ep->geometry().local(faceglobal);
	  int numberInSelf = is->indexInInside();	  
			
	  s.xi[numberInSelf]=faceglobal[0];

	  if (numberInSelf == p) // left face
	  { 		 
	    j.x[0][0] = faceglobal[0];
	    j.x[0][1] = j.x[0][0];
	    //RT satJ=0;
	    if(is->neighbor())
	    {
	      ep=is->outside();
	      indexi=elementmapper.map(*ep);
	      j.Sat=this->transproblem.variables.saturation[indexi]; 
	    }
	    else
	    {
	      //gdat[p] = faceglobal[0];

	      const Dune::FieldVector<ct,dim>& facelocalDim = Dune::ReferenceElements<ct,dim>::general(gtf).position(numberInSelf,1);
	      j.Sat =this->transproblem.dirichlet(faceglobal,*ep,facelocalDim);
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
    else
      calcm();

    return 0;

  }


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

    //Bilde Ausgangselement, ausgehend von dem Element dessen Sättigung berechnet werden soll
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

    //Es wird so lange die Stromlinie zurückverfolgt, wie ein Partikel in der Zeit t zurücklegen konnte. 
    //lambda stellt die höchste "Geschwindigkeit" eines Partikels dar ( max(f'(u)) )
    while(xi > -lambda*dt)
    {    
      //std::cout<<xi<<std::endl;
      Dune::FieldVector<VType,2*dim> xx;
      Dune::FieldVector<VType,2*dim> vx;
     
      //Setzen der Positionen und Geschwindigkeiten an den Elementkanten
      IntersectionIterator endis=IntersectionIteratorGetter<G,LevelTag>::end(*ep);
      for(IntersectionIterator is=IntersectionIteratorGetter<G,LevelTag>::begin(*ep);is!=endis;++is)
      {
	GeometryType gtf = is->geometryInInside().type();
	const Dune::FieldVector<ct,dim-1>& facelocal = Dune::ReferenceElements<RT,dim-1>::general(gtf).position(0,0);
	const VType& faceglobal = is->geometry().global(facelocal);

	int numberInSelf = is->indexInInside();
	xx[numberInSelf] = faceglobal;
      
	vx[numberInSelf] = this->transproblem.variables.vTotal(*ep,numberInSelf);
       	vx[numberInSelf] /= this->transproblem.porosity();
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
      int m = 0;                  //Kantenindex
      RT dte = dtx[m];            //Time-of-Flight
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
      

      for(IntersectionIterator is = IntersectionIteratorGetter<G,LevelTag>::begin(*ep);is!=endis;++is)
      {
	int numberInSelf = is->indexInInside();
	if(numberInSelf == kant[m])
	  if(is->neighbor())
	  {
	    ep = is->outside();
	
	    if(n > 0) 
	    {	
	      for(int i=0;i<dim;++i)
		if(i!=m && fabs(dtx[i]-dte)<1e-7)
		{
		  IntersectionIterator endir = IntersectionIteratorGetter<G,LevelTag>::end(*ep);
		  for(IntersectionIterator ir = IntersectionIteratorGetter<G,LevelTag>::begin(*ep);ir!=endir;++ir)
		  {
		    int nIS = ir->indexInInside();
		    
		    if(nIS == kant[i])
		      if(ir->neighbor())
		      {
			ep = ir->outside();
			indexi = elementmapper.map(*ep);
			j.Sat = this->transproblem.variables.saturation[indexi];		      
		   
			s.index=indexi;
		      }
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
	    else 
	    {
	      indexi = elementmapper.map(*ep);
	      j.Sat = this->transproblem.variables.saturation[indexi];	      
	      
	      s.index=indexi;
	    }
	  }//if is->neighbor()
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


  template<class G,class RT,class VC>
  int ChTransport<G,RT,VC>::fronttracking(const RT dt,const RT v = 1)
  {       

    //berechne Postion der Schocks zum  Zeitpunkt dt
    ListIterator lit;
  
    for(lit=ch.begin();lit!=ch.end();++lit)
      for(int q=0;q<2;++q)
	lit->x[1][q]=lit->x[0][q]+lit->m[q]*dt*v;
    //da eventuell Makrozeitschritte berechnet werden: 
    //tau=Makrozeitpunkt, t1=vorherige Makrozeitpunkt,dtau=tau-t1
    //x[0] gibt die Position der Shocks zum Zeitpunkt t1 an
    //x[1] gibt die Position der Shocks zum Zeipunkt tau an


    RT t1=0;
    while(t1<dt)       
    {   
      RT dtau=dt;
      ListIterator merk;
      //berechne kleinsten Zeitschritt, bei dem sich die Shocks schneiden. 
      //merk gibt an, zu welcher Zelle die entsprechenden Shocks gehören
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
	RT m;
	
	//   std::cout<<merk->Sat<<" ; "<<merk2->Sat<<std::endl;
	
	solveRP(m,chrar,merk->Sat,merk2->Sat);
	
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



 template<class G, class RT, class VC>
  int ChTransport<G,RT,VC>::blau(RepresentationType& updateVec, int indexi)
  {    
     
    RT sat = 0;
    
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
	RT xt0=lit->x[1][0];
	RT xt1=lit->x[1][1];
	  	  
	if(xt0 <= 0 && 0 <= xt1)      
	{
	  RT xt05=0.5*(xt0+xt1);
	  RT st05=lit->Sat;
	
	  ListIterator merk(lit),merk2(lit);
	  merk--;
	  merk2++;
	  RT sx0,sx1,xx0,xx1;
	  
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
    updateVec[indexi]=sat - this->transproblem.variables.saturation[indexi];
   

    return 0;
  }


  template<class G, class RT, class VC>
  int ChTransport<G,RT,VC>::rot(RepresentationType& updateVec,RepresentationType& dxtt)
  {

    typedef typename std::template list<slNode>::iterator  SlistIterator;
    for(SlistIterator slit=sl.begin();slit!=sl.end();++slit)
    {
      RT x0=slit->xi[0];
      RT x1=slit->xi[1];

      RT dx=fabs(x1-x0);
      RT sat=0;
      RT dxt=0;
  
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
    if(t==0)
      // dt=8.64e6;
      // dt = 2.16e7;
        dt = 1e8;
      //dt = 4.32e7;

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

     
    Iterator eendit = this->grid.template lend<0>(this->level());
    if(G::dimension == 1)
    {
      Iterator it = this->grid.template lbegin<0>(this->level());
      init(dt,maxa,it);

	       
      RT fv = this->transproblem.variables.vTotal(*it,0)[0];
      fv /= this->transproblem.porosity();     

      fronttracking(dt,fv);
      rot(updateVec,dxtt);	       
	//blau nur für streamlines möglich
	}
    else 
    
      for (Iterator it = this->grid.template lbegin<0>(this->level()); it != eendit; ++it)
      {
	int indexi = elementmapper.map(*it);
   
	calcStreamline(dt,1.4*maxa,it);

	if(ch.size()>0)
	{
//  	  typedef typename std::template list<slNode>::iterator SlistIterator;
//    	  for(SlistIterator slit=sl.begin();slit!=sl.end();++slit)
//           {
//     	std::cout<<"index:"<<slit->index<<std::endl;
//    	std::cout<<"xi:"<<slit->xi<<std::endl;
//     	std::cout<<"Sat:"<<this->transproblem.variables.saturation[slit->index]<<std::endl;
//           }
	  fronttracking(dt);	 
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
    
#define IMP
#ifdef IMP 

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
#else
		char impFileName[120];
		strcpy(impFileName, "data2.dgf");
		
		FILE* dataFile;
		dataFile=fopen(impFileName, "r");
		//      dataFile.open(impFileName.c_str());
		for(int i=0;i<this->transproblem.variables.saturation.size();++i)
		  fscanf(dataFile,"%lf",&this->transproblem.variables.saturation[i]);
		fclose(dataFile);  
#endif
    return;
  }
  


}
#endif
