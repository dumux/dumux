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
      typedef OwnerOverlapCopyCommunication<GlobalIdType,LocalIdType> ParentType;

      typedef typename OwnerOverlapCopyAttributeSet::AttributeSet AttributeSet;
      typedef typename ParentType::BC   BC; // BufferedCommunicator
      typedef typename ParentType::IF   IF; // Interface

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
      communicator.template forward<ParentType::template AddGatherScatter<T> >(source,dest);
      communicator.free();
    }
  public:

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
