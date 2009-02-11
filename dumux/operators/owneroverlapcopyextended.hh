// $Id$
#ifndef DUNE_OWNEROVERLAPCOPYEXTENDED_HH
#define DUNE_OWNEROVERLAPCOPYEXTENDED_HH

#include<dune/istl/owneroverlapcopy.hh>

namespace Dune {

#if HAVE_MPI

  /**
   * @brief A class setting up standard communication for a two-valued
   * attribute set with owner/overlap/copy semantics.
   *
   * set up communication from known distribution with owner/overlap/copy semantics
   */
  template <class GlobalIdType, class LocalIdType>
  class OwnerOverlapCopyExtendedCommunication : public OwnerOverlapCopyCommunication<GlobalIdType,LocalIdType>
  {
    typedef typename IndexInfoFromGrid<GlobalIdType,LocalIdType>::IndexTripel IndexTripel;
    typedef typename IndexInfoFromGrid<GlobalIdType,LocalIdType>::RemoteIndexTripel RemoteIndexTripel;
    typedef typename std::set<IndexTripel>::const_iterator localindex_iterator;
    typedef typename std::set<RemoteIndexTripel>::const_iterator remoteindex_iterator;
    typedef typename OwnerOverlapCopyAttributeSet::AttributeSet AttributeSet;
    typedef Dune::ParallelLocalIndex<AttributeSet> LI;
    typedef Dune::ParallelIndexSet<GlobalIdType,LI,512> PIS;
    typedef Dune::RemoteIndices<PIS> RI;
    typedef Dune::RemoteIndexListModifier<PIS,false> RILM;
    typedef typename RI::RemoteIndex RX;
    typedef Dune::BufferedCommunicator<PIS> BC;
    typedef Dune::Interface<PIS> IF;

    template<typename T>
    struct AddGatherScatter
    {
      typedef typename CommPolicy<T>::IndexedType V;

      static V gather(const T& a, int i)
      {
        return a[i];
      }

      static void scatter(T& a, V v, int i)
      {
        a[i] += v;
      }
    };

    void buildAllToAllInterface () const
    {
      if (AllToAllInterfaceBuilt)
        AllToAllInterface.free();
      typedef Combine<EnumItem<AttributeSet,OwnerOverlapCopyAttributeSet::owner>,EnumItem<AttributeSet,OwnerOverlapCopyAttributeSet::copy>,AttributeSet> OwnerCopySet;
      typedef Combine<EnumItem<AttributeSet,OwnerOverlapCopyAttributeSet::owner>,EnumItem<AttributeSet,OwnerOverlapCopyAttributeSet::overlap>,AttributeSet> OwnerOverlapSet;
      typedef Combine<OwnerOverlapSet,EnumItem<AttributeSet,OwnerOverlapCopyAttributeSet::copy>,AttributeSet> AllSet;
      AllSet sourceFlags;
      AllSet destFlags;
      AllToAllInterface.build(this->remoteIndices(),sourceFlags,destFlags);
      AllToAllInterfaceBuilt = true;
    }

  public:
    /**
     * @brief Communicate values from all data points to all other data points and add them to those values.
     *
     * @brief source The data to send from.
     * @brief dest The data to add them communicated values to.
     */
    template<class T>
    void addAllToAll (const T& source, T& dest) const
    {
      if (!AllToAllInterfaceBuilt)
        buildAllToAllInterface ();
      BC communicator;
      communicator.template build<T>(AllToAllInterface);
      communicator.template forward<AddGatherScatter<T> >(source,dest);
      communicator.free();
    }

    typedef Dune::EnumItem<AttributeSet,OwnerOverlapCopyAttributeSet::copy> CopyFlags;

    /** @brief The type of the parallel index set. */
    typedef Dune::ParallelIndexSet<GlobalIdType,LI,512> ParallelIndexSet;

    /** @brief The type of the remote indices. */
    typedef Dune::RemoteIndices<PIS> RemoteIndices;

    /**
     * @brief The type of the reverse lookup of indices. */
    typedef Dune::GlobalLookupIndexSet<ParallelIndexSet> GlobalLookupIndexSet;


    /**
     * @brief Construct the communication without any indices.
     *
     * The local index set and the remote indices have to be set up
     * later on.
     * @param comm_ The MPI Communicator to use, e. g. MPI_COMM_WORLD
     */
    OwnerOverlapCopyExtendedCommunication (MPI_Comm comm_)
      : OwnerOverlapCopyCommunication<GlobalIdType,LocalIdType> (comm_), AllToAllInterfaceBuilt(false)
    {}

    /**
     * @brief Constructor
     * @param indexinfo The set of IndexTripels describing the local and remote indices.
     * @param comm_ The communicator to use in the communication.
     */
    OwnerOverlapCopyExtendedCommunication (const IndexInfoFromGrid<GlobalIdType,LocalIdType>& indexinfo, MPI_Comm comm_)
      : OwnerOverlapCopyCommunication<GlobalIdType,LocalIdType>(indexinfo, comm_), AllToAllInterfaceBuilt(false)
    {}

    // destructor: free memory in some objects
    ~OwnerOverlapCopyExtendedCommunication ()
    {
      if (AllToAllInterfaceBuilt) AllToAllInterface.free();
    }

  private:
    mutable IF AllToAllInterface;
    mutable bool AllToAllInterfaceBuilt;
  };

#endif


  /** @} end documentation */

} // end namespace

#endif
