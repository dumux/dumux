// $Id$

#ifndef DUNE_P1FUNCTIONEXTENDED_HH
#define DUNE_P1FUNCTIONEXTENDED_HH

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
#include<dune/disc/functions/functions.hh>
#include<dune/disc/functions/p0function.hh>

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

template<class G>
class LeafCommunicate
{
public:
    LeafCommunicate (const G& g)
        : grid(g)
    {}

    template<class DataHandle>
    void communicate (DataHandle& data, InterfaceType iftype, CommunicationDirection dir) const
    {
        grid.template communicate<DataHandle>(data,iftype,dir);
    }

private:
    const G& grid;
};

/** \todo Please doc me! */

template<class G>
class LevelCommunicate
{
public:
    LevelCommunicate (const G& g, int l)
        : grid(g), level(l)
    {}

    template<class DataHandle>
    void communicate (DataHandle& data, InterfaceType iftype, CommunicationDirection dir) const
    {
        grid.template communicate<DataHandle>(data,iftype,dir,level);
    }

private:
    const G& grid;
    int level;
};

//! compute 1-overlap on non-overlapping grid
template<class G, class GV, class VM, class LC>
class P1ExtendOverlap {

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

    /** \todo Please doc me! */

    // A DataHandle class to exchange border rows
    class OwnerExchange
        : public CommDataHandleIF<OwnerExchange,Pair> {
    public:
        //! export type of data for message buffer
        typedef Pair DataType;

        //! returns true if data for this codim should be communicated
        bool contains (int dim, int codim) const
        {
            return (codim==dim);
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
                    typename std::map<IdType,int>::const_iterator oIt = owner.find(*i);
                    if (oIt != owner.end())
                        {
                            buff.write(Pair(oIt->first,oIt->second));
                        }
                }
        }

        /*! unpack data from message buffer to user

          n is the number of objects sent by the sender
        */
        template<class MessageBuffer, class EntityType>
        void scatter (MessageBuffer& buff, const EntityType& e, size_t n)
        {
            //int alpha=vertexmapper.map(e);
            for (size_t i=0; i<n; i++)
                {
                    Pair x;
                    buff.read(x);
                    typename std::map<IdType,int>::const_iterator oIt = owner.find(x.first);
                    if (oIt != owner.end())
                        {
                            owner[x.first] = std::max(owner[x.first], x.second);
                        }
                }
        }

        //! constructor
        OwnerExchange (const G& g, const VM& vm, std::map<int,GIDSet>& ids, std::map<IdType,int>& o, const std::map<int,IdType>& i2g)
            : grid(g), vertexmapper(vm), myids(ids), owner(o), index2gid(i2g)
        {}

    private:
        const G& grid;
        const VM& vertexmapper;
        std::map<int,GIDSet>& myids;
        std::map<IdType,int>& owner;
        const std::map<int,IdType>& index2gid;
    };

    /** \todo Please doc me! */

    // A DataHandle class to exchange border rows
    class IdExchange
        : public CommDataHandleIF<IdExchange,Pair> {
    public:
        //! export type of data for message buffer
        typedef Pair DataType;

        //! returns true if data for this codim should be communicated
        bool contains (int dim, int codim) const
        {
            return (codim==dim);
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
                    int rankIdx = grid.comm().rank();
                    typename std::map<IdType,int>::const_iterator g2iIt = gid2index.find(*i);
                    int localIdx = g2iIt->second;
                    typename std::map<int,int>::const_iterator d2OIt = doubled2Original.find(localIdx);
                    if (d2OIt != doubled2Original.end())
                        localIdx = d2OIt->second;
                    if (partitionType[localIdx] == GhostEntity) {
                        rankIdx = -1;//grid.comm().size();
                    }
                    buff.write(Pair(*i,rankIdx)); // I have these global ids
                    owner[*i] = rankIdx;
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
                        owner[x.first] = std::max(owner[x.first],source);
                }
            myprocs[alpha].insert(source);
        }

        //! constructor
        IdExchange (const G& g, const GV& gv, const VM& vm, std::map<int,GIDSet>& ids, std::map<int,ProcSet>& procs,
                    std::map<IdType,int>& o, const std::map<IdType,int>& g2i, const std::map<int,IdType>& i2g, int eDOFs,
                    const std::map<int,int>& d2O)
            : grid(g), gridview(gv), vertexmapper(vm), myids(ids), myprocs(procs), owner(o),
              gid2index(g2i), index2gid(i2g), extraDOFs(eDOFs), doubled2Original(d2O)
        {
            partitionType.resize(vertexmapper.size()+extraDOFs);
            VIterator vendit = gridview.template end<n>();
            for (VIterator it = gridview.template begin<n>(); it!=vendit; ++it)
                {
                    partitionType[vertexmapper.map(*it)] = it->partitionType();
                }

            //          for (typename std::map<IdType,int>::const_iterator it = gid2index.begin(); it != gid2index.end(); ++it)
            //              std::cout << grid.comm().rank() << ": global = " << it->first << ", local = " << it->second << std::endl;

            //          for (typename std::map<int,IdType>::const_iterator it = index2gid.begin(); it != index2gid.end(); ++it)
            //              std::cout << grid.comm().rank() << ": local = " << it->first << ", global = " << it->second << std::endl;
            //
            //          for (typename std::map<int,int>::const_iterator it = doubled2Original.begin(); it != doubled2Original.end(); ++it)
            //              std::cout << grid.comm().rank() << ": doubled = " << it->first << ", original = " << it->second << std::endl;
        }

    private:
        const G& grid;
        const GV& gridview;
        const VM& vertexmapper;
        std::map<int,GIDSet>& myids;
        std::map<int,ProcSet>& myprocs;
        std::map<IdType,int>& owner;
        const std::map<IdType,int>& gid2index;
        const std::map<int,IdType>& index2gid;
        int extraDOFs;
        const std::map<int,int>& doubled2Original;
        std::vector<PartitionType> partitionType;
    };

    /** \todo Please doc me! */

    // A DataHandle class to exchange border rows
    class BorderLinksExchange
        : public CommDataHandleIF<BorderLinksExchange,IdType>{
    public:
        //! export type of data for message buffer
        typedef IdType DataType;

        //! returns true if data for this codim should be communicated
        bool contains (int dim, int codim) const
        {
            return (codim==dim);
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

    typedef IndexInfoFromGrid<IdType,int> P1IndexInfoFromGrid;

    //! fill data structure with information needed by ISTL
    void fillIndexInfoFromGrid (const G& grid, const GV& gridview, const VM& vertexmapper, P1IndexInfoFromGrid& info)
    {

        std::map<IdType,int> gid2index;
        std::map<int,IdType> index2gid;
        std::map<int,GIDSet> borderlinks;
        int extraDOFs = 0;
        extend(grid, gridview, vertexmapper, borderlinks, extraDOFs, gid2index, index2gid);

        std::map<int,int> doubled2Original;
        VIterator ghostEndIt = gridview.template end<n>();
        for (VIterator ghostIt = gridview.template begin<n>(); ghostIt!=ghostEndIt; ++ghostIt)
            {
                if (ghostIt->partitionType() != GhostEntity)
                    continue;

                int ghostIndex = vertexmapper.map(*ghostIt);

                const IdType& ghostGlobalId = index2gid[ghostIndex];

                VIterator realEndIt = gridview.template end<n>();
                for (VIterator realIt = gridview.template begin<n>(); realIt!=realEndIt; ++realIt)
                    {
                        if (realIt->partitionType() == InteriorEntity)
                            continue;

                        int realIndex = vertexmapper.map(*realIt);

                        if (realIndex == ghostIndex)
                            continue;

                        const IdType& realGlobalId = index2gid[realIndex];

                        if (realGlobalId == ghostGlobalId) {
                            doubled2Original[std::max(realIndex, ghostIndex)] = std::min(realIndex, ghostIndex);
                            break;
                        }
                    }
            }
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
                Dune::GeometryType gt = it->type();
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
        IdExchange datahandle(grid,gridview,vertexmapper,myids,myprocs,owner,gid2index,index2gid,extraDOFs,doubled2Original);
        lc.template communicate<IdExchange>(datahandle,InteriorBorder_InteriorBorder_Interface,ForwardCommunication);

        OwnerExchange ownerDatahandle(grid,vertexmapper,myids,owner,index2gid);
        lc.template communicate<OwnerExchange>(ownerDatahandle,InteriorBorder_InteriorBorder_Interface,ForwardCommunication);

        //        for (typename std::map<IdType,int>::iterator i=owner.begin(); i!=owner.end(); ++i)
        //            std::cout << grid.comm().rank() << ": global = " << i->first << ", owner = " << i->second << std::endl;


        extraDOFs = 0;
        for (typename std::map<int,GIDSet>::iterator i=myids.begin(); i!=myids.end(); ++i)
            for (typename GIDSet::iterator j=(i->second).begin(); j!=(i->second).end(); ++j)
                gid2index[*j] = -1; // indicates "not assigned yet"
        VIterator vendit = gridview.template end<n>();
        for (VIterator it = gridview.template begin<n>(); it!=vendit; ++it)
            {
                IdType beta = grid.globalIdSet().id(*it);
                if (gid2index.find(beta)!=gid2index.end())
                    {
                        int alpha = vertexmapper.map(*it);
                        gid2index[beta] = alpha; // assign existing local index
                        index2gid[alpha] = beta;
                    }
            }
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
        std::set< tuple<IdType,int,int> > ownindices;
        for (typename std::map<int,GIDSet>::iterator i=myids.begin(); i!=myids.end(); ++i)
            for (typename GIDSet::iterator j=(i->second).begin(); j!=(i->second).end(); ++j)
                {
                    int a=slave;
                    if (owner[*j]==grid.comm().rank()) a=master;
                    int localIndex = gid2index[*j];
                    // do not add double indices
                    if (doubled2Original.find(localIndex) == doubled2Original.end())
                        info.addLocalIndex(tuple<IdType,int,int>(*j,localIndex,a));
                    else
                        info.addLocalIndex(tuple<IdType,int,int>(*j,doubled2Original[localIndex],a));
                }
        std::set< tuple<int,IdType,int> > remoteindices;
        for (typename std::map<int,ProcSet>::iterator i=myprocs.begin(); i!=myprocs.end(); ++i)
            {
                GIDSet& thisset = myids[i->first];
                for (typename GIDSet::iterator j=thisset.begin(); j!=thisset.end(); ++j)
                    for (typename ProcSet::iterator p=(i->second).begin(); p!=(i->second).end(); ++p)
                        {
                            int a=slave;
                            if (owner[*j]==(*p)) a=master;
                            if (*p!=grid.comm().rank()) info.addRemoteIndex(tuple<int,IdType,int>(*p,*j,a));
                        }
            }

        //        for (typename std::set< tuple<IdType,int,int> >::const_iterator it = info.localIndices().begin();
        //            it != info.localIndices().end(); ++it)
        //            std::cout << grid.comm().rank() << ": global = " << Element<0>::get(*it)
        //                << ", local = " << Element<1>::get(*it) << ", m = " << Element<2>::get(*it) << std::endl;
        //
        //        for (typename std::set< tuple<int,IdType,int> >::const_iterator it = info.remoteIndices().begin();
        //            it != info.remoteIndices().end(); ++it)
        //            std::cout << grid.comm().rank() << ": p = " << Element<0>::get(*it)
        //                << ", global = " << Element<1>::get(*it) << ", m = " << Element<2>::get(*it) << std::endl;

        // clear what is not needed anymore to save memory
        myids.clear();
        gid2index.clear();
        myprocs.clear();
        owner.clear();
        neighbors.clear();

        return;
    }

    //! fill data structures needed for extension
    void extend (const G& grid, const GV& gridView, const VM& vertexmapper,
                 std::map<int,GIDSet>& borderlinks, int& extraDOFs, std::map<IdType,int>& gid2index)
    {
        // initialize output parameters
        borderlinks.clear();
        extraDOFs = 0;
        gid2index.clear();

        // build local borderlinks from mesh
        Iterator eendit = gridView.template end<0>();
        for (Iterator it = gridView.template begin<0>(); it!=eendit; ++it)
            {
                Dune::GeometryType gt = it->type();
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
        VIterator vendit = gridView.template end<n>();
        for (VIterator it = gridView.template begin<n>(); it!=gridView.template end<n>(); ++it)
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

    //! fill data structures needed for extension
    void extend (const G& grid, const GV& gridView, const VM& vertexmapper,
                 std::map<int,GIDSet>& borderlinks, int& extraDOFs, std::map<IdType,int>& gid2index,
                 std::map<int,IdType>& index2gid)
    {
        // initialize output parameters
        borderlinks.clear();
        extraDOFs = 0;
        gid2index.clear();

        // build local borderlinks from mesh
        Iterator eendit = gridView.template end<0>();
        for (Iterator it = gridView.template begin<0>(); it!=eendit; ++it)
            {
                Dune::GeometryType gt = it->type();
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
        VIterator vendit = gridView.template end<n>();
        for (VIterator it = gridView.template begin<n>(); it!=gridView.template end<n>(); ++it)
            {
                IdType beta = grid.globalIdSet().id(*it);
                if (gid2index.find(beta)!=gid2index.end())
                    {
                        int alpha = vertexmapper.map(*it);
                        gid2index[beta] = alpha;
                        index2gid[alpha] = beta;
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

    P1ExtendOverlap (LC lcomm)
        : lc(lcomm)
    {}

private:
    LC lc;
};

/** \todo Please doc me! */

// forward declaration
template<class G, class RT> class P1FunctionExtendedManager;


//! class for P1 finite element functions on a grid
/*! This class implements the interface of a DifferentiableGridFunction
  with piecewise linear elements using a Lagrange basis. It is implemented
  using the general shape functions, thus it should work for all element types
  and dimensions.

  In addition to the DifferentiableGridFunction interface P1 functions can be initialized
  from a C0GridFunction via Lagrange interpolation. Dereferencing delivers
  the coefficient vector.

  \tparam GV The grid view
*/
template<class GV, class RT, class LC, int m=1>
class P1FunctionExtended
    : virtual public GridFunctionGlobalEvalDefault<GV,RT,m>
    , virtual public FunctionDefault<typename GV::Grid::ctype,RT,GV::Grid::dimension,m>
    , virtual public ElementwiseCInfinityFunction<GV,RT,m>
    , virtual public H1Function<typename GV::Grid::ctype,RT,GV::Grid::dimension,m>
    , virtual public C0GridFunction<typename GV::Grid,RT,m>
{
    //! get the Grid from the GridView
    typedef typename GV::Grid G;

    //! get domain field type from the grid
    typedef typename G::ctype DT;

    //! get domain dimension from the grid
    enum {n=G::dimension};

    //! get entity from the grid
    typedef typename G::template Codim<0>::Entity Entity;

    //! Parameter for mapper class
    template<int dim>
    struct P1Layout
    {
        bool contains (Dune::GeometryType gt)
        {
            return gt.dim() == 0;
        }
    };

    //! make copy constructor private
    P1FunctionExtended (const P1FunctionExtended&);

    // types
    typedef typename G::Traits::GlobalIdSet IDS;
    typedef typename IDS::IdType IdType;
    typedef std::set<IdType> GIDSet;
    typedef typename GV::IndexSet IS;

public:
    typedef FieldVector<RT,m> BlockType;
    typedef BlockVector<BlockType> RepresentationType;
    typedef MultipleCodimMultipleGeomTypeMapper<G,IS,P1Layout> VM;
    typedef typename P1ExtendOverlap<G,GV,VM,LC>::P1IndexInfoFromGrid P1IndexInfoFromGrid;

    //! allocate data
    P1FunctionExtended (const GV& gridView, LC lcomm, bool extendoverlap=false)
        : GridFunctionGlobalEvalDefault<GV,RT,m>(gridView)
        , gv(gridView), mapper_(gridView.grid(),gridView.indexSet()), lc(lcomm), oldcoeff(0)
    {
        // check if overlap extension is possible
        if (extendoverlap && gv.grid().overlapSize(0)>0)
            DUNE_THROW(GridError,"P1FunctionExtended: extending overlap requires nonoverlapping grid");

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
                P1ExtendOverlap<G,GV,VM,LC> extender(lc);
                extender.extend(gv.grid(),gv,mapper_,borderlinks,extraDOFs,gid2index);
            }

        // allocate the vector
        oldcoeff = 0;
        try {
            coeff = new RepresentationType(mapper_.size()+extraDOFs);
        }
        catch (std::bad_alloc) {
            std::cerr << "not enough memory in P1FunctionExtended" << std::endl;
            throw; // rethrow exception
        }
        dverb << "making FE function with " << mapper_.size()+extraDOFs << " components"
              << "(" << extraDOFs << " extra degrees of freedom)" << std::endl;
    }

    //! deallocate the vector
    ~P1FunctionExtended ()
    {
        delete coeff;
        if (oldcoeff!=0) delete oldcoeff;
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
        Dune::GeometryType gt = e.type(); // extract type of element
        for (int i=0; i<Dune::LagrangeShapeFunctions<DT,RT,n>::general(gt,1).size(); ++i)
            value += Dune::LagrangeShapeFunctions<DT,RT,n>::general(gt,1)[i].evaluateFunction(0,xi)*(*coeff)[mapper_.template map<n>(e,i)][comp];
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
        Dune::GeometryType gt = e.type(); // extract type of element
        y = 0;
        for (int i=0; i<Dune::LagrangeShapeFunctions<DT,RT,n>::general(gt,1).size(); ++i)
            {
                RT basefuncvalue=Dune::LagrangeShapeFunctions<DT,RT,n>::general(gt,1)[i].evaluateFunction(0,xi);
                int index = mapper_.template map<n>(e,i);
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
        Dune::GeometryType gt = e.type(); // extract type of element
        const typename Dune::LagrangeShapeFunctionSetContainer<DT,RT,n>::value_type&
            sfs=Dune::LagrangeShapeFunctions<DT,RT,n>::general(gt,1);
        Dune::FieldMatrix<DT,n,n> jac = e.geometry().jacobianInverseTransposed(xi);
        for (int i=0; i<sfs.size(); ++i)
            {
                Dune::FieldVector<DT,n> grad(0),temp;
                for (int l=0; l<n; l++)
                    temp[l] = sfs[i].evaluateDerivative(0,l,xi);
                jac.umv(temp,grad); // transform gradient to global ooordinates
                value += grad[dir] * (*coeff)[mapper_.template map<n>(e,i)][comp];
            }
        return value;
    }


    //! interpolate nodal values from a grid function
    /*! Lagrange interpolation of a P1 finite element function from given
      continuous grid function. Evaluation is done by visiting the vertices
      of each element and storing a bitvector of visited vertices.

      @param[in]  u    a continuous grid function
    */
    void interpolate (const C0GridFunction<G,RT,m>& u)
    {
        typedef typename GV::template Codim<0>::Iterator Iterator;
        std::vector<bool> visited(mapper_.size());
        for (int i=0; i<mapper_.size(); i++) visited[i] = false;

        Iterator eendit = gv.template end<0>();
        for (Iterator it = gv.template begin<0>(); it!=eendit; ++it)
            {
                Dune::GeometryType gt = it->type();
                for (int i=0; i<Dune::LagrangeShapeFunctions<DT,RT,n>::general(gt,1).size(); ++i)
                    if (!visited[mapper_.template map<n>(*it,i)])
                        {
                            for (int c=0; c<m; c++)
                                (*coeff)[mapper_.template map<n>(*it,i)][c] =
                                    u.evallocal(c,*it,Dune::LagrangeShapeFunctions<DT,RT,n>::general(gt,1)[i].position());
                            visited[mapper_.template map<n>(*it,i)] = true;
                        }
            }
    }

    /** \todo Please doc me! */
    void interpolate (const P0Function<GV,RT,m>& u)
    {
        typedef typename GV::template Codim<0>::Iterator Iterator;
        std::vector<char> counter(mapper_.size());
        for (int i=0; i<mapper_.size(); i++) counter[i] = 0;

        for (int i=0; i<(*coeff).size(); i++)
            (*coeff)[i] = 0;
        Iterator eendit = gv.template end<0>();
        for (Iterator it = gv.template begin<0>(); it!=eendit; ++it)
            {
                Dune::GeometryType gt = it->type();
                for (int i=0; i<Dune::LagrangeShapeFunctions<DT,RT,n>::general(gt,1).size(); ++i)
                    {
                        for (int c=0; c<m; c++)
                            (*coeff)[mapper_.template map<n>(*it,i)][c] +=
                                u.evallocal(c,*it,Dune::LagrangeShapeFunctions<DT,RT,n>::general(gt,1)[i].position());
                        counter[mapper_.template map<n>(*it,i)] += 1;
                    }
            }
        for (int i=0; i<counter.size(); i++)
            for (int c=0; c<m; c++)
                (*coeff)[i][c] /= counter[i];
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
    void fillIndexInfoFromGrid (P1IndexInfoFromGrid& info)
    {
        P1ExtendOverlap<G,GV,VM,LC> extender(lc);
        extender.fillIndexInfoFromGrid(gv.grid(),gv,mapper_,info);
    }


    /** empty method to maintain symmetry
        For vertex data nothing is required in preAdapt but for other
        finite element functions this method is necessary.
    */
    void preAdapt ()
    {
    }

    /** @brief Initiate update process

        Call this method after the grid has been adapted. The representation is
        now updated to the new grid and the finite element function can be used on
        the new grid. However the data is not initialized.
        The old representation (with respect to the old grid) can still be accessed if
        it has been saved. It is deleted in endUpdate().
    */
    void postAdapt (P1FunctionExtendedManager<G,RT>& manager)
    {
        typedef typename G::template Codim<n>::LeafIterator VLeafIterator;
        typedef typename G::template Codim<0>::LeafIterator ELeafIterator;
        typedef typename G::template Codim<n>::LevelIterator VLevelIterator;
        typedef typename G::template Codim<0>::LevelIterator ELevelIterator;
        typedef typename G::template Codim<0>::EntityPointer EEntityPointer;

        // \todo check that function is only called for data with respect to leafs
        // save the current representation
        oldcoeff = coeff;

        // allow mapper to recompute its internal sizes
        mapper_.update();

        // overlap extension, recompute extra DOFs
        if (extendOverlap)
            {
                // set of neighbors in global ids for border vertices
                std::map<int,GIDSet> borderlinks;
                std::map<IdType,int> gid2index;
                extraDOFs = 0;

                // compute extension
                P1ExtendOverlap<G,GV,VM,LC> extender(lc);
                extender.extend(gv.grid(),gv,mapper_,borderlinks,extraDOFs,gid2index);
            }

        // allocate data with new size (while keeping the old data ...)
        try {
            coeff = new RepresentationType(mapper_.size()+extraDOFs); // allocate new representation
        }
        catch (std::bad_alloc) {
            std::cerr << "not enough memory in P1FunctionExtended update" << std::endl;
            throw; // rethrow exception
        }
        std::cout << "P1  function enlarged to " << mapper_.size() << " components" << std::endl;

        // vector of flags to store which vertex has been handled already
        std::vector<bool> visited(mapper_.size());
        for (int i=0; i<mapper_.size(); i++) visited[i] = false;

        // now loop over the NEW mesh to copy the data that was already in the OLD mesh
        VLeafIterator veendit = gv.grid().template leafend<n>();
        for (VLeafIterator it = gv.grid().template leafbegin<n>(); it!=veendit; ++it)
            {
                // lookup in mapper
                int i;
                if (manager.savedMap().contains(*it,i))
                    {
                        //               std::cout << " found vertex=" << it->geometry().corner(0)
                        //                         << " at i=" << i
                        //                         << " oldindex=" << manager.oldIndex()[i]
                        //                         << " newindex=" << mapper_.map(*it)
                        //                         << std::endl;
                        // the vertex existed already in the old mesh, copy data
                        for (int c=0; c<m; c++)
                            (*coeff)[mapper_.map(*it)][c] = (*oldcoeff)[manager.oldIndex()[i]][c];
                        // ... and mark as visited
                        visited[mapper_.map(*it)] = true;
                    }
            }

        // now loop the second time to interpolate the new coefficients
        // new implementation using interpolation on codim 0
        for (int level=1; level<=gv.grid().maxLevel(); ++level)
            {
                ELevelIterator elendit = gv.grid().template lend<0>(level);
                for (ELevelIterator it = gv.grid().template lbegin<0>(level); it!=elendit; ++it)
                    {
                        GeometryType gte = it->type();
                        for (int i=0; i<Dune::LagrangeShapeFunctions<DT,RT,n>::general(gte,1).size(); ++i)
                            {
                                int index = mapper_.template map<n>(*it,i);
                                if (!visited[index])
                                    {
                                        // OK, this is a new vertex
                                        EEntityPointer father=it->father(); // the father element
                                        GeometryType gtf = father->type(); // fathers type
                                        const FieldVector<DT,n>& cpos=Dune::LagrangeShapeFunctions<DT,RT,n>::general(gte,1)[i].position();
                                        FieldVector<DT,n> pos = it->geometryInFather().global(cpos); // map corner to father element
                                        for (int c=0; c<m; c++)
                                            (*coeff)[index][c] = 0;
                                        for (int j=0; j<Dune::LagrangeShapeFunctions<DT,RT,n>::general(gtf,1).size(); ++j)
                                            {
                                                RT basefuncvalue = Dune::LagrangeShapeFunctions<DT,RT,n>::general(gtf,1)[j].evaluateFunction(0,pos);
                                                for (int c=0; c<m; c++)
                                                    (*coeff)[index][c] += basefuncvalue * (*coeff)[mapper_.template map<n>(*father,j)][c];
                                                //                       std::cout << "  corner=" << i
                                                //                                 << " cpos=" << father->geometry().corner(i)
                                                //                                 << " u=" << (*coeff)[mapper_.template map<n>(*father,i)]
                                                //                                 << std::endl;
                                            }
                                        //                   std::cout << "index=" << mapper_.map(*it) << " value=" << value << std::endl;
                                        visited[index] = true;
                                    }
                            }
                    }
            }

        // now really delete old representation
        if (oldcoeff!=0) delete oldcoeff;
        oldcoeff = 0;
    }

    /** @brief export the mapper for external use
     */
    const VM& mapper () const
    {
        return mapper_;
    }

    //! VTK output
    void vtkout (VTKWriter<GV>& vtkwriter, std::string s) const
    {
        typename VTKWriter<GV>::VTKFunction *p = new VTKGridFunctionWrapper<GV,RT,m>(*this,s);
        vtkwriter.addVertexData(p);
    }

private:
    // The GridView.  Don't use a reference here since Grid::levelView
    // returns a temporary
    const GV gv;

    // we need a mapper
    VM mapper_;

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


/** \brief P1 finite element function on the leaf grid

    \param G The grid
    \param RT The type used for the component values of the function
    \param m Vector-valued functions: number of components
*/
template<class G, class RT, int m=1>
class LeafP1FunctionExtended : public P1FunctionExtended<typename G::LeafGridView,RT,LeafCommunicate<G>,m>
{
public:
    /** \brief Constructor for a given grid
        \todo Please doc the second argument
    */
    LeafP1FunctionExtended (const G& grid, bool extendoverlap=false)
        : GridFunctionGlobalEvalDefault<typename G::LeafGridView,RT,m>(grid.leafView())
        , P1FunctionExtended<typename G::LeafGridView,RT,LeafCommunicate<G>,m>(grid.leafView(),LeafCommunicate<G>(grid),extendoverlap)
    {}
};


/** \brief P1 finite element function on a given level grid

    \param G The grid
    \param RT The type used for the component values of the function
    \param m Vector-valued functions: number of components
*/
template<class G, class RT, int m=1>
class LevelP1FunctionExtended : public P1FunctionExtended<typename G::LevelGridView,RT,LevelCommunicate<G>,m>
{
public:
    /** \brief Constructor for a given grid
        \todo Please doc the third argument
    */
    LevelP1FunctionExtended (const G& grid, int level, bool extendoverlap=false)
        : GridFunctionGlobalEvalDefault<typename G::LevelGridView,RT,m>(grid.levelView(level))
        , P1FunctionExtended<typename G::LevelGridView,RT,LevelCommunicate<G>,m>(grid.levelView(level),LevelCommunicate<G>(grid,level),extendoverlap)
    {}
};


/** \brief Manage mesh adaptation and load balancing for several P1 finite element functions

    Adaptivity management is only required for the leaf finite element functions,
    therefore we do only allow those to be registered.

    \todo If we would have P1 vector valued functions with different
    numbers of components we still would need a separate manager for every size, although
    this is actually not necessary.
*/
template<class G, class RT>
class P1FunctionExtendedManager {
    enum {dim=G::dimension};
    typedef typename G::ctype DT;
    typedef typename G::template Codim<dim>::LeafIterator VLeafIterator;
    typedef typename G::template Codim<0>::EntityPointer EEntityPointer;

    template<int dim>
    struct P1Layout
    {
        bool contains (Dune::GeometryType gt)
        {
            return gt.dim() == 0;
        }
    };

public:

    //! manages nothing
    P1FunctionExtendedManager (const G& g) : mapper(g,g.leafIndexSet()), grid(g), savedmap(g)
    {
        // allocate index array to correct size (this possible for vertex data)
        oldindex.resize(mapper.size());

        // and allocate the universal mapper to acces the old indices
        savedmap.clear(); //should be empty already

        // now loop over all vertices and copy the index provided by the mapper
        VLeafIterator veendit = grid.template leafend<dim>();
        for (VLeafIterator it = grid.template leafbegin<dim>(); it!=veendit; ++it)
            {
                oldindex[savedmap.map(*it)] = mapper.map(*it);
            }
    }

    const GlobalUniversalMapper<G>& savedMap ()
    {
        return savedmap;
    }

    const std::vector<int>& oldIndex ()
    {
        return oldindex;
    }

private:
    // we need a mapper
    MultipleCodimMultipleGeomTypeMapper<G,typename G::template Codim<0>::LeafIndexSet,P1Layout> mapper;

    // store a reference to the grid that is managed
    const G& grid;

    // We need a persistent consecutive enumeration
    GlobalUniversalMapper<G> savedmap;

    // The old leaf indices are stored in a dynamically allocated vector
    std::vector<int> oldindex;
};


/** @} */

}
#endif
