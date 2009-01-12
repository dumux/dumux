// $Id$

#ifndef DUNE_MIXEDFUNCTION_HH
#define DUNE_MIXEDFUNCTION_HH

//C++ includes
#include<new>
#include<iostream>
#include<vector>
#include<list>
#include<map>
#include<set>

// Dune includes
#include<dune/common/fvector.hh>
#include<dune/common/exceptions.hh>
#include<dune/common/tuples.hh>
#include<dune/common/stdstreams.hh>
#include<dune/grid/common/grid.hh>
#include<dune/grid/common/mcmgmapper.hh>
#include<dune/grid/common/universalmapper.hh>
#include<dune/grid/io/file/vtk/vtkwriter.hh>
#include<dune/istl/bvector.hh>
#include<dune/istl/operators.hh>
#include<dune/istl/bcrsmatrix.hh>
#include<dune/istl/owneroverlapcopy.hh>
#include<dune/disc/shapefunctions/lagrangeshapefunctions.hh>
// same directory includes
#include<dune/disc/functions/functions.hh>
#include<dune/disc/functions/p0function.hh>
#include"dumux/functions/p1functionextended.hh"
#include"dumux/shapefunctions/CRshapefunctions.hh"

/**
* @file
* @brief  defines a class for piecewise linear finite element functions
* @author Peter Bastian
*/
namespace Dune
{
/** @addtogroup DISC_Functions
*
* @{
*/
/**
* @brief defines a class for piecewise linear finite element functions
*
*/

//! compute 1-overlap on non-overlapping grid
template<class G, class GV, class VM, class LC>
class MixedExtendOverlap {

    // types
    typedef typename G::ctype DT;
    enum {n=G::dimension};
    typedef typename GV::IndexSet IS;
    typedef typename G::template Codim<0>::Entity Entity;
    typedef typename GV::template Codim<0>::Iterator Iterator;
    typedef typename GV::template Codim<n>::Iterator VIterator;
    typedef typename G::template Codim<0>::EntityPointer EEntityPointer;
    typedef typename G::Traits::GlobalIdSet IDS;
    typedef typename IDS::IdType IdType;
    typedef std::set<IdType> GIDSet;
    typedef std::pair<IdType,int> Pair;
    typedef std::set<int> ProcSet;

    // A DataHandle class to exchange border rows
    class IdExchange
    : public CommDataHandleIF<IdExchange,Pair> {
    public:
        //! export type of data for message buffer
        typedef Pair DataType;

        //! returns true if data for this codim should be communicated
        bool contains (int dim, int codim) const
        {
            return (codim==1);
        }

        //! returns true if size per entity of given dim and codim is a constant
        bool fixedsize (int dim, int codim) const
        {
            return false;
        }

        /*! how many objects of type DataType have to be sent for a given entity

        Note: Only the sender side needs to know this size.
        */
        template<class EntityType>
        size_t size (EntityType& e) const
        {
            return myids[vertexmapper.map(e)].size();
        }

        //! pack data from user to message buffer
        template<class MessageBuffer, class EntityType>
        void gather (MessageBuffer& buff, const EntityType& e) const
        {
            int alpha=vertexmapper.map(e);
            GIDSet& thisset = myids[alpha];
            for (typename GIDSet::iterator i=thisset.begin(); i!=thisset.end(); ++i)
            {
                buff.write(Pair(*i,grid.comm().rank())); // I have these global ids
                owner[*i] = grid.comm().rank();
            }
            myprocs[alpha].insert(grid.comm().rank());
        }

        /*! unpack data from message buffer to user

        n is the number of objects sent by the sender
        */
        template<class MessageBuffer, class EntityType>
        void scatter (MessageBuffer& buff, const EntityType& e, size_t n)
        {
            int alpha=vertexmapper.map(e);
            GIDSet& thisset = myids[alpha];
            int source;
            for (size_t i=0; i<n; i++)
            {
                Pair x;
                buff.read(x);
                thisset.insert(x.first);
                source=x.second;
                if (owner.find(x.first)==owner.end())
                    owner[x.first] = source;
                else
                    owner[x.first] = std::min(owner[x.first],source);
            }
            myprocs[alpha].insert(source);
        }

        //! constructor
        IdExchange (const G& g, const VM& vm, std::map<int,GIDSet>& ids, std::map<int,ProcSet>& procs,
                std::map<IdType,int>& o)
        : grid(g), vertexmapper(vm), myids(ids), myprocs(procs), owner(o)
        {}

    private:
        const G& grid;
        const VM& vertexmapper;
        std::map<int,GIDSet>& myids;
        std::map<int,ProcSet>& myprocs;
        std::map<IdType,int>& owner;
    };

    // A DataHandle class to exchange border rows
    class BorderLinksExchange
    : public CommDataHandleIF<BorderLinksExchange,IdType>{
    public:
        //! export type of data for message buffer
        typedef IdType DataType;

        //! returns true if data for this codim should be communicated
        bool contains (int dim, int codim) const
        {
            return (codim==1);
        }

        //! returns true if size per entity of given dim and codim is a constant
        bool fixedsize (int dim, int codim) const
        {
            return false;
        }

        /*! how many objects of type DataType have to be sent for a given entity

        Note: Only the sender side needs to know this size.
        */
        template<class EntityType>
        size_t size (EntityType& e) const
        {
            return borderlinks[vertexmapper.map(e)].size();
        }

        //! pack data from user to message buffer
        template<class MessageBuffer, class EntityType>
        void gather (MessageBuffer& buff, const EntityType& e) const
        {
            GIDSet& myset = borderlinks[vertexmapper.map(e)];
            for (typename GIDSet::iterator i=myset.begin(); i!=myset.end(); ++i)
                buff.write(*i);
        }

        /*! unpack data from message buffer to user

        n is the number of objects sent by the sender
        */
        template<class MessageBuffer, class EntityType>
        void scatter (MessageBuffer& buff, const EntityType& e, size_t n)
        {
            GIDSet& myset = borderlinks[vertexmapper.map(e)];
            for (size_t i=0; i<n; i++)
            {
                DataType x;
                buff.read(x);
                myset.insert(x);
            }
        }

        //! constructor
        BorderLinksExchange (const G& g, std::map<int,GIDSet>& bl, const VM& vm)
        : grid(g), borderlinks(bl), vertexmapper(vm)
        {}

    private:
        const G& grid;
        std::map<int,GIDSet>& borderlinks;
        const VM& vertexmapper;
    };

    public:

        enum Attributes {slave=OwnerOverlapCopyAttributeSet::copy,
            master=OwnerOverlapCopyAttributeSet::owner,
            overlap=OwnerOverlapCopyAttributeSet::overlap};

        typedef IndexInfoFromGrid<IdType,int> MixedIndexInfoFromGrid;

        //! fill data structure with information needed by ISTL
        void fillIndexInfoFromGrid (const G& grid, const GV& gridview, const VM& vertexmapper, MixedIndexInfoFromGrid& info)
        {
            // build a map of sets where each local index is assigned
            // a set of global ids which are neighbors of this vertex
            // At the same time assign to each local index to a set of processors
            // and at the same time determine the owner of the gid
            std::map<int,GIDSet> myids;
            std::map<int,ProcSet> myprocs;
            std::map<IdType,int> owner;
            Iterator eendit = gridview.template end<0>();
            for (Iterator it = gridview.template begin<0>(); it!=eendit; ++it)
            {
                Dune::GeometryType gt = it->geometry().type();
                const typename Dune::ReferenceElementContainer<DT,n>::value_type&
                refelem = ReferenceElements<DT,n>::general(gt);

                if (it->partitionType()==InteriorEntity)
                    for (int i=0; i<refelem.size(n); i++)
                    {
                        if (it->template entity<n>(i)->partitionType()==BorderEntity)
                        {
                            int alpha = vertexmapper.template map<n>(*it,i);
                            GIDSet& thisset = myids[alpha];
                            for (int j=0; j<refelem.size(n); j++)
                            {
                                IdType beta = grid.globalIdSet().template subId<n>(*it,j);
                                thisset.insert(beta);
                            }
                        }
                    }
            }
            IdExchange datahandle(grid,vertexmapper,myids,myprocs,owner);
            lc.template communicate<IdExchange>(datahandle,InteriorBorder_InteriorBorder_Interface,ForwardCommunication);

            // build map from global id to local index
            std::map<IdType,int> gid2index;
            for (typename std::map<int,GIDSet>::iterator i=myids.begin(); i!=myids.end(); ++i)
                for (typename GIDSet::iterator j=(i->second).begin(); j!=(i->second).end(); ++j)
                    gid2index[*j] = -1; // indicates "not assigned yet"
            VIterator vendit = gridview.template end<n>();
            for (VIterator it = gridview.template begin<n>(); it!=gridview.template end<n>(); ++it)
            {
                IdType beta = grid.globalIdSet().id(*it);
                if (gid2index.find(beta)!=gid2index.end())
                {
                    int alpha = vertexmapper.map(*it);
                    gid2index[beta] = alpha; // assign existing local index
                }
            }
            int extraDOFs = 0;
            for (typename std::map<IdType,int>::iterator i=gid2index.begin(); i!=gid2index.end(); ++i)
                if (i->second==-1)
                {
                    i->second = vertexmapper.size()+extraDOFs; // assign new local index
                    extraDOFs++;
                }

            // build a set of all neighboring processors
            ProcSet neighbors;
            for (typename std::map<int,ProcSet>::iterator i=myprocs.begin(); i!=myprocs.end(); ++i)
                for (typename ProcSet::iterator j=(i->second).begin(); j!=(i->second).end(); ++j)
                    if (*j!=grid.comm().rank())
                        neighbors.insert(*j);

            // now all the necessary information is in place

            // application: for all neighbors build a list of global ids
            for (typename ProcSet::iterator p=neighbors.begin(); p!=neighbors.end(); ++p)
            {
                GIDSet remote;
                for (typename std::map<int,ProcSet>::iterator i=myprocs.begin(); i!=myprocs.end(); ++i)
                    if ((i->second).find(*p)!=(i->second).end())
                    {
                        GIDSet& thisset = myids[i->first];
                        for (typename GIDSet::iterator j=thisset.begin(); j!=thisset.end(); ++j)
                            remote.insert(*j);
                    }
            }


            // fill the info object
            std::set< Tuple<IdType,int,int> > ownindices;
            for (typename std::map<int,GIDSet>::iterator i=myids.begin(); i!=myids.end(); ++i)
                for (typename GIDSet::iterator j=(i->second).begin(); j!=(i->second).end(); ++j)
                {
                    int a=slave;
                    if (owner[*j]==grid.comm().rank()) a=master;
                    info.addLocalIndex(Tuple<IdType,int,int>(*j,gid2index[*j],a));
                }
            std::set< Tuple<int,IdType,int> > remoteindices;
            for (typename std::map<int,ProcSet>::iterator i=myprocs.begin(); i!=myprocs.end(); ++i)
            {
                GIDSet& thisset = myids[i->first];
                for (typename GIDSet::iterator j=thisset.begin(); j!=thisset.end(); ++j)
                    for (typename ProcSet::iterator p=(i->second).begin(); p!=(i->second).end(); ++p)
                    {
                        int a=slave;
                        if (owner[*j]==(*p)) a=master;
                        if (*p!=grid.comm().rank()) info.addRemoteIndex(Tuple<int,IdType,int>(*p,*j,a));
                    }
            }

            // clear what is not needed anymore to save memory
            myids.clear();
            gid2index.clear();
            myprocs.clear();
            owner.clear();
            neighbors.clear();

            return;
        }


        //! fill data structures needed for extension
        void extend (const G& grid, const GV& gridview, const VM& vertexmapper,
                std::map<int,GIDSet>& borderlinks, int& extraDOFs, std::map<IdType,int>& gid2index)
        {
            // initialize output parameters
            borderlinks.clear();
            extraDOFs = 0;
            gid2index.clear();

            // build local borderlinks from mesh
            Iterator eendit = gridview.template end<0>();
            for (Iterator it = gridview.template begin<0>(); it!=eendit; ++it)
            {
                Dune::GeometryType gt = it->geometry().type();
                const typename Dune::ReferenceElementContainer<DT,n>::value_type&
                refelem = ReferenceElements<DT,n>::general(gt);

                // generate set of neighbors in global ids for border vertices
                if (it->partitionType()==InteriorEntity)
                    for (int i=0; i<refelem.size(n); i++)
                        if (it->template entity<n>(i)->partitionType()==BorderEntity)
                        {
                            int alpha = vertexmapper.template map<n>(*it,i);
                            GIDSet& myset = borderlinks[alpha];
                            for (int j=0; j<refelem.size(n); j++)
                                if (i!=j)
                                {
                                    IdType beta = grid.globalIdSet().template subId<n>(*it,j);
                                    myset.insert(beta);
                                    //                           std::cout << g.comm().rank() << ": "
                                    //                                     << "borderlink " << alpha
                                    //                                     << " " << vertexmapper.template map<n>(*it,j)
                                    //                                     << " " << beta
                                    //                                     << std::endl;
                                }
                        }
            }

            // exchange neighbor info for border vertices
            BorderLinksExchange datahandle(grid,borderlinks,vertexmapper);
            lc.template communicate<BorderLinksExchange>(datahandle,
                    InteriorBorder_InteriorBorder_Interface,
                    ForwardCommunication);

            // initialize inverse map with ids we have
            for (typename std::map<int,GIDSet>::iterator i=borderlinks.begin(); i!=borderlinks.end(); ++i)
                for (typename GIDSet::iterator j=(i->second).begin(); j!=(i->second).end(); ++j)
                    gid2index[*j] = -1;

            // check with ids we already have in the grid to find out extra vertices
            VIterator vendit = gridview.template end<n>();
            for (VIterator it = gridview.template begin<n>(); it!=gridview.template end<n>(); ++it)
            {
                IdType beta = grid.globalIdSet().id(*it);
                if (gid2index.find(beta)!=gid2index.end())
                {
                    int alpha = vertexmapper.map(*it);
                    gid2index[beta] = alpha;
                }
            }

            // assign index to extra DOFs
            extraDOFs = 0;
            for (typename std::map<IdType,int>::iterator i=gid2index.begin(); i!=gid2index.end(); ++i)
                if (i->second==-1)
                {
                    i->second = vertexmapper.size()+extraDOFs;
                    extraDOFs++;
                }

            //           for (typename std::map<int,GIDSet>::iterator i=borderlinks.begin(); i!=borderlinks.end(); ++i)
            //             for (typename GIDSet::iterator j=(i->second).begin(); j!=(i->second).end(); ++j)
            //               std::cout << grid.comm().rank() << ": " << "comm borderlink " << i->first
            //                         << " " << gid2index[*j] << " " << *j << std::endl;
        }

        MixedExtendOverlap (LC lcomm)
        : lc(lcomm)
        {}

    private:
        LC lc;
};

//! class for Mixed finite element functions on a grid
/*! This class implements the interface of a DifferentiableGridFunction
with piecewise linear elements using a Mixed basis. It is implemented
using the general shape functions, thus it should work for all element types
and dimensions.

In addition to the DifferentiableGridFunction interface Mixed functions can be initialized
from a C0GridFunction via Mixed interpolation. Dereferencing delivers
the coefficient vector.
*/
template<class G, class RT, class GV, class LC, int m=1>
class MixedFunction : virtual public ElementwiseCInfinityFunction<G,RT,m>,
virtual public H1Function<typename G::ctype,RT,G::dimension,m>,
virtual public C0GridFunction<G,RT,m>
{
    //! get domain field type from the grid
    typedef typename G::ctype DT;

    //! get domain dimension from the grid
    enum {n=G::dimension};

    //! get entity from the grid
    typedef typename G::template Codim<0>::Entity Entity;
    typedef typename GV::IndexSet IS;

    //! Parameter for mapper class
    template<int dim>
    struct ElementAndFaceLayout
    {
        bool contains (Dune::GeometryType gt)
        {
            return (gt.dim() == n-1 || gt.dim() == n);
        }
    };

    //! make copy constructor private
    MixedFunction (const MixedFunction&);

    // types
    typedef typename G::Traits::GlobalIdSet IDS;
    typedef typename IDS::IdType IdType;
    typedef std::set<IdType> GIDSet;

public:
    typedef FieldVector<RT,m> BlockType;
    typedef BlockVector<BlockType> RepresentationType;
    typedef MultipleCodimMultipleGeomTypeMapper<G,IS,ElementAndFaceLayout> ElementAndFaceMapper;
    typedef typename MixedExtendOverlap<G,GV,ElementAndFaceMapper,LC>::MixedIndexInfoFromGrid MixedIndexInfoFromGrid;

    //! allocate data
    MixedFunction (const G& g, const GV& gv, LC lcomm, bool extendoverlap=false)
    : grid_(g), gridview(gv), is(gv.indexSet()), mapper_(g,is), lc(lcomm), oldcoeff(0)
    {
        // check if overlap extension is possible
        if (extendoverlap && g.overlapSize(0)>0)
            DUNE_THROW(GridError,"MixedFunction: extending overlap requires nonoverlapping grid");

        // no extra DOFs so far
        extraDOFs = 0;
        extendOverlap = extendoverlap;

        // overlap extension
        if (extendoverlap)
        {
            // set of neighbors in global ids for border vertices
            std::map<int,GIDSet> borderlinks;
            std::map<IdType,int> gid2index;

            // compute extension
            MixedExtendOverlap<G,GV,ElementAndFaceMapper,LC> extender(lc);
            extender.extend(g,gridview,mapper_,borderlinks,extraDOFs,gid2index);
        }

        // allocate the vector
        oldcoeff = 0;
        try {
            coeff = new RepresentationType(mapper_.size()+extraDOFs);
        }
        catch (std::bad_alloc) {
            std::cerr << "not enough memory in MixedFunction" << std::endl;
            throw; // rethrow exception
        }
        dverb << "making FE function with " << mapper_.size()+extraDOFs << " components"
        << "(" << extraDOFs << " extra degrees of freedom)" << std::endl;
    }

    //! deallocate the vector
    ~MixedFunction ()
    {
        delete coeff;
        if (oldcoeff!=0) delete oldcoeff;
    }

    //! evaluate single component comp at global point x
    /*! Evaluate a single component of the vector-valued
    function.
    @param[in] comp number of component to be evaluated
    @param[in] x    position to be evaluated
    \return         value of the component
    \todo Not implemented yet!
    */
    virtual RT eval (int comp, const Dune::FieldVector<DT,n>& x) const
    {
        DUNE_THROW(NotImplemented, "global eval not implemented yet");
        return 0;
    }

    //! evaluate all components at point x and store result in y
    /*! Evaluation function for all components at once.
    @param[in]  x    position to be evaluated
    @param[out] y    result vector to be filled
    \todo Not implemented yet!
    */
    virtual void evalall (const Dune::FieldVector<DT,n>& x, Dune::FieldVector<RT,m>& y) const
    {
        DUNE_THROW(NotImplemented, "global eval not implemented yet");
    }

    //! evaluate partial derivative
    /*! Evaluate partial derivative of a component of the vector-valued function.
    @param[in]  comp    number of component that should be differentiated
    @param[in]  d       vector giving order of derivative for each variable
    @param[in]  x       position where derivative is to be evaluated
    \return             value of the derivative
    \todo Not implemented yet!
    */
    virtual RT derivative (int comp, const Dune::FieldVector<int,n>& d, const Dune::FieldVector<DT,n>& x) const
    {
        DUNE_THROW(NotImplemented, "global derivative not implemented yet");
    }

    //! return number of partial derivatives that can be taken
    /*! A DifferentiableFunction can say how many derivatives exist
    and can be safely evaluated.
    */
    virtual int order () const
    {
        return 1; // up to now only one derivative is implemented
    }

    //! evaluate single component comp in the entity e at local coordinates xi
    /*! Evaluate the function in an entity at local coordinates.
    @param[in]  comp   number of component to be evaluated
    @param[in]  e      reference to grid entity of codimension 0
    @param[in]  xi     point in local coordinates of the reference element of e
    \return            value of the component
    */
    virtual RT evallocal (int comp, const Entity& e, const Dune::FieldVector<DT,n>& xi) const
    {
        RT value=0;
        Dune::GeometryType gt = e.geometry().type(); // extract type of element
        for (int i=0; i<Dune::CRShapeFunctions<DT,RT,n>::general(gt,1).size(); ++i)
            value += Dune::CRShapeFunctions<DT,RT,n>::general(gt,1)[i].evaluateFunction(0,xi)*(*coeff)[mapper_.template map<1>(e,i)][comp];
            return value;
    }

    //! evaluate all components  in the entity e at local coordinates xi
    /*! Evaluates all components of a function at once.
    @param[in]  e      reference to grid entity of codimension 0
    @param[in]  xi     point in local coordinates of the reference element of e
    @param[out] y      vector with values to be filled
    */
    virtual void evalalllocal (const Entity& e, const Dune::FieldVector<DT,G::dimension>& xi,
            Dune::FieldVector<RT,m>& y) const
            {
        Dune::GeometryType gt = e.geometry().type(); // extract type of element
        y = 0;
        for (int i=0; i<Dune::CRShapeFunctions<DT,RT,n>::general(gt,1).size(); ++i)
        {
            RT basefuncvalue=Dune::CRShapeFunctions<DT,RT,n>::general(gt,1)[i].evaluateFunction(0,xi);
            int index = mapper_.template map<1>(e,i);
            for (int c=0; c<m; c++)
                y[c] += basefuncvalue * (*coeff)[index][c];
        }
            }

    //! evaluate derivative in local coordinates
    /*! Evaluate the partial derivative a the given position
    in local coordinates in an entity.
    @param[in]  comp    number of component that should be differentiated
    @param[in]  d       vector giving order of derivative for each variable
    @param[in]  e       reference to grid entity of codimension 0
    @param[in]  xi      point in local coordinates of the reference element of e
    \return             value of the derivative
    */
    virtual RT derivativelocal (int comp, const Dune::FieldVector<int,n>& d,
            const Entity& e, const Dune::FieldVector<DT,n>& xi) const
            {
        int dir=-1;
        int order=0;
        for (int i=0; i<n; i++)
        {
            order += d[i];
            if (d[i]>0) dir=i;
        }
        assert(dir != -1);
        if (order!=1) DUNE_THROW(GridError,"can only evaluate one derivative");

        RT value=0;
        Dune::GeometryType gt = e.geometry().type(); // extract type of element
        const typename Dune::CRShapeFunctionSetContainer<DT,RT,n>::value_type&
        sfs=Dune::CRShapeFunctions<DT,RT,n>::general(gt,1);
        Dune::FieldMatrix<DT,n,n> jac = e.geometry().jacobianInverseTransposed(xi);
        for (int i=0; i<sfs.size(); ++i)
        {
            Dune::FieldVector<DT,n> grad(0),temp;
            for (int l=0; l<n; l++)
                temp[l] = sfs[i].evaluateDerivative(0,l,xi);
            jac.umv(temp,grad); // transform gradient to global ooordinates
            value += grad[dir] * (*coeff)[mapper_.template map<1>(e,i)][comp];
        }
        return value;
            }

    //! return const reference to coefficient vector
    /*! Dereferencing a finite element function returns the
    coefficient representation of the finite element function.
    This is the const version.
    */
    const RepresentationType& operator* () const
    {
        return (*coeff);
    }

    //! return reference to coefficient vector
    /*! Dereferencing a finite element function returns the
    coefficient representation of the finite element function.
    This is the non-const version.
    */
    RepresentationType& operator* ()
    {
        return (*coeff);
    }


    //! deliver communication object
    void fillIndexInfoFromGrid (MixedIndexInfoFromGrid& info)
    {
        MixedExtendOverlap<G,GV,ElementAndFaceMapper,LC> extender(lc);
        extender.fillIndexInfoFromGrid(grid_,gridview,mapper_,info);
    }


    /** empty method to maintain symmetry
    For vertex data nothing is required in preAdapt but for other
    finite element functions this method is necessary.
    */
    void preAdapt ()
    {
    }


private:
    // a reference to the grid
    const G& grid_;

    // reference to index set on the grid (might be level or leaf)
    const GV& gridview;
    const IS& is;

    // we need a mapper
    ElementAndFaceMapper mapper_;

    // level or leafwise communication object
    LC lc;

    // extra DOFs from extending nonoverlapping to overlapping grid
    int extraDOFs;
    bool extendOverlap;

    // and a dynamically allocated vector
    RepresentationType* coeff;

    // saved pointer in update phase
    RepresentationType* oldcoeff;
};


/** \brief Mixed finite element function on the leaf grid

\param G The grid
\param RT The type used for the component values of the function
\param m Vector-valued functions: number of components
*/
template<class G, class RT, int m=1>
class LeafMixedFunction : public MixedFunction<G,RT,typename G::LeafGridView,LeafCommunicate<G>,m>
{
public:
    /** \brief Constructor for a given grid
    \todo Please doc the second argument
    */
    LeafMixedFunction (const G& grid, bool extendoverlap=false)
    : MixedFunction<G,RT,typename G::LeafGridView,LeafCommunicate<G>,m>(grid,grid.leafView(),LeafCommunicate<G>(grid),extendoverlap)
    {}
};


/** \brief Mixed finite element function on a given level grid

\param G The grid
\param RT The type used for the component values of the function
\param m Vector-valued functions: number of components
*/
template<class G, class RT, int m=1>
class LevelMixedFunction : public MixedFunction<G,RT,typename G::LevelGridView,LevelCommunicate<G>,m>
{
public:
    /** \brief Constructor for a given grid
    \todo Please doc the third argument
    */
    LevelMixedFunction (const G& grid, int level, bool extendoverlap=false)
    : MixedFunction<G,RT,typename G::LevelGridView,LevelCommunicate<G>,m>(grid,grid.levelView(level),LevelCommunicate<G>(grid,level),extendoverlap)
    {}
};

/** @} */

}
#endif
