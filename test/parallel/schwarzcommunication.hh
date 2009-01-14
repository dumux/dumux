#ifndef SCHWARZCOMMUNICATION_HH
#define SCHWARZCOMMUNICATION_HH

#include<dune/common/exceptions.hh>
#include<dune/grid/common/grid.hh>
#include<dune/grid/common/mcmgmapper.hh>

/**
 * @brief A class setting up communication needed for overlapping Schwarz methods
 * based on the grid communication.
 *
 * It works on leaf and level-wise grids for P1 finite elements.
 */
template <class Grid, class IndexSet, class CommWrapper>
class P1OverlappingSchwarzCommunication
{
  // extract constants
  enum {dim=Grid::dimension};

  // Parameter for mapper class
  template<int dim>
  struct P1Layout
  {
    bool contains (Dune::GeometryType gt)
    {
      return gt.isVertex();
    }
  };

  // extract types
  typedef typename Grid::ctype ct;
  typedef typename IndexSet::template Codim<dim>::template Partition<Dune::All_Partition>::Iterator VIterator;
  typedef Dune::MultipleCodimMultipleGeomTypeMapper<Grid,IndexSet,P1Layout> Mapper;

  // A DataHandle class that computes owner for each vertex
  class OwnerExchange
    : public Dune::CommDataHandleIF<OwnerExchange,int>
  {
  public:
    //! export type of data for message buffer
    typedef int DataType;

    //! returns true if data for this codim should be communicated
    bool contains (int d, int c) const
    {
      return (c==d);
    }

    //! returns true if size per entity of given dim and codim is a constant
    bool fixedsize (int d, int c) const
    {
      return true;
    }

    /*! how many objects of type DataType have to be sent for a given entity

    Note: Only the sender side needs to know this size.
    */
    template<class EntityType>
    size_t size (EntityType& e) const
    {
      return 1;
    }

    //! pack data from user to message buffer
    template<class MessageBuffer, class EntityType>
    void gather (MessageBuffer& buff, const EntityType& e) const
    {
      // just send your rank
      buff.write(grid.comm().rank());
    }

    /*! unpack data from message buffer to user

    n is the number of objects sent by the sender
    */
    template<class MessageBuffer, class EntityType>
    void scatter (MessageBuffer& buff, const EntityType& e, size_t n)
    {
      int q;
      buff.read(q); // get rank of sender

      // smallest rank is the owner
      if (grid.comm().rank()>q)
        {
          int index=mapper.map(e);
          uniquemask[index] = false;
        }
    }

    //! constructor
    OwnerExchange (const Grid& grid_, const Mapper& mapper_, std::vector<bool>& uniquemask_)
      : grid(grid_), mapper(mapper_), uniquemask(uniquemask_)
    {}

  private:
    const Grid& grid;
    const Mapper& mapper;
    std::vector<bool>& uniquemask;
  };


  // A DataHandle class that sends and adds
  template<class T>
  class AddingDataHandle
    : public Dune::CommDataHandleIF<AddingDataHandle<T>,
                    typename T::value_type>
  {
  public:
    //! export type of data for message buffer
    typedef typename T::value_type DataType;

    //! returns true if data for this codim should be communicated
    bool contains (int d, int c) const
    {
      return (c==d);
    }

    //! returns true if size per entity of given dim and codim is a constant
    bool fixedsize (int d, int c) const
    {
      return true;
    }

    /*! how many objects of type DataType have to be sent for a given entity

    Note: Only the sender side needs to know this size.
    */
    template<class EntityType>
    size_t size (EntityType& e) const
    {
      return 1;
    }

    //! pack data from user to message buffer
    template<class MessageBuffer, class EntityType>
    void gather (MessageBuffer& buff, const EntityType& e) const
    {
      // get index of entity
      int index=mapper.map(e);
      buff.write(source[index]);
    }

    /*! unpack data from message buffer to user

    n is the number of objects sent by the sender
    */
    template<class MessageBuffer, class EntityType>
    void scatter (MessageBuffer& buff, const EntityType& e, size_t n)
    {
      DataType x;
      buff.read(x);
      int index=mapper.map(e);
      dest[index] += x;
    }

    //! constructor
    AddingDataHandle (const Grid& grid_, const Mapper& mapper_,
                              const T& source_, T& dest_)
      : grid(grid_), mapper(mapper_), source(source_), dest(dest_)
    {}

  private:
    const Grid& grid;
    const Mapper& mapper;
    const T& source;
    T& dest;
  };

  // A DataHandle class that sends and adds
  template<class T>
  class CopyOwnerDataHandle
    : public Dune::CommDataHandleIF<CopyOwnerDataHandle<T>,
                    typename T::value_type>
  {
  public:
    //! export type of data for message buffer
    typedef typename T::value_type DataType;

    //! returns true if data for this codim should be communicated
    bool contains (int d, int c) const
    {
      return (c==d);
    }

    //! returns true if size per entity of given dim and codim is a constant
    bool fixedsize (int d, int c) const
    {
      return true;
    }

    /*! how many objects of type DataType have to be sent for a given entity

    Note: Only the sender side needs to know this size.
    */
    template<class EntityType>
    size_t size (EntityType& e) const
    {
      return 1;
    }

    //! pack data from user to message buffer
    template<class MessageBuffer, class EntityType>
    void gather (MessageBuffer& buff, const EntityType& e) const
    {
      // get index of entity
      int index=mapper.map(e);
      if (uniquemask[index])
        buff.write(source[index]);
      else
        {
          dest[index] = 0;
          buff.write(dest[index]);
        }
    }

    /*! unpack data from message buffer to user

    n is the number of objects sent by the sender
    */
    template<class MessageBuffer, class EntityType>
    void scatter (MessageBuffer& buff, const EntityType& e, size_t n)
    {
      DataType x;
      buff.read(x);
      int index=mapper.map(e);
      dest[index] += x;
    }

    //! constructor
    CopyOwnerDataHandle (const Grid& grid_, const Mapper& mapper_,
                         const T& source_, T& dest_,
                         const std::vector<bool>& uniquemask_)
      : grid(grid_), mapper(mapper_), source(source_), dest(dest_), uniquemask(uniquemask_)
    {}

  private:
    const Grid& grid;
    const Mapper& mapper;
    const T& source;
    T& dest;
    const std::vector<bool>& uniquemask;
  };

public:
  enum{
    category = Dune::SolverCategory::overlapping
      };

  const typename Grid::template Codim<0>::CollectiveCommunication& communicator() const
  {
    return grid.comm();
  }

  /**
   * @brief Communicate values from owner data points to all other data points.
   *
   * @brief source The data to send from.
   * @brief dest The data to send to.
   */
  template<class T>
  void copyOwnerToAll (const T& source, T& dest) const
  {
    CopyOwnerDataHandle<T> datahandle(grid,mapper,source,dest,uniquemask);
    commwrapper.template communicate<CopyOwnerDataHandle<T> >(datahandle,Dune::All_All_Interface,
                                     Dune::ForwardCommunication);
  }


  /**
   * @brief Communicate values from owner data points to all other data points and add them to those values.
   *
   * @brief source The data to send from.
   * @brief dest The data to add them communicated values to.
   */
  template<class T>
  void addOwnerOverlapToAll (const T& source, T& dest) const
  {
    AddingDataHandle<T> datahandle(grid,mapper,source,dest);
    commwrapper.template communicate<AddingDataHandle<T> >(datahandle,Dune::All_All_Interface,
                                                    Dune::ForwardCommunication);
  }

  /**
   * @brief Compute a global dot product of two vectors.
   *
   * @param x The first vector of the product.
   * @param y The second vector of the product.
   * @param res Reference to store the result in.
   */
  template<class T1, class T2>
  void dot (const T1& x, const T1& y, T2& result) const
  {
    if (uniquemask.size()!=static_cast<typename std::vector<bool>::size_type>(x.size()))
      DUNE_THROW(Dune::RangeError,"size mismatch in dot");

    result = 0;
    for (int i=0; i<x.size(); i++)
      if (uniquemask[i])
        result += x[i]*y[i];
    result = grid.comm().sum(result);
    return;
  }

  /**
   * @brief Compute the global euclidian norm of a vector.
   *
   * @param x The vector to compute the norm of.
   * @return The global euclidian norm of that vector.
   */
  template<class T1>
  double norm (const T1& x) const
  {
    if (uniquemask.size()!=static_cast<typename std::vector<bool>::size_type>(x.size()))
      DUNE_THROW(Dune::RangeError,"size mismatch in norm");

    double result = 0;
    for (int i=0; i<x.size(); i++)
      if (uniquemask[i])
        result += x[i].two_norm2();
    return sqrt(grid.comm().sum(result));
  }


  /**
   * @brief Set vector to zero at front nodes
   *
   * @param x The vector to project.
   */
  template<class T1>
  void project (T1& x) const
  {
     if (nonfrontmask.size()!=static_cast<typename std::vector<bool>::size_type>(x.size()))
      DUNE_THROW(Dune::RangeError,"size mismatch in project");
    for (int i=0; i<x.size(); i++)
      if (nonfrontmask[i]==false)
        x[i] = 0;
  }


  /**
   * @brief Constructor
   * @param indexinfo The set of IndexTripels describing the local and remote indices.
   * @param comm_ The communicator to use in the communication.
   */
  P1OverlappingSchwarzCommunication (const Grid& grid_, const IndexSet& indexset_,
                                     const CommWrapper& commwrapper_)
    : grid(grid_), indexset(indexset_), commwrapper(commwrapper_), mapper(grid_,indexset_)
  {
    // construct mask vector holding 1 for all vertices that are neither front nor ghost
    nonfrontmask.resize(mapper.size());
    VIterator vendit = indexset.template end<dim,Dune::All_Partition>();
    for (VIterator it = indexset.template begin<dim,Dune::All_Partition>(); it!=vendit; ++it)
      if (it->partitionType()!=Dune::FrontEntity && it->partitionType()!=Dune::GhostEntity)
        nonfrontmask[mapper.map(*it)] = true;
      else
        nonfrontmask[mapper.map(*it)] = false;

    // construct mask vector holding 1 for all vertices assigned uniquely to this process
    uniquemask.resize(mapper.size());
    for (VIterator it = indexset.template begin<dim,Dune::All_Partition>(); it!=vendit; ++it)
      if (it->partitionType()==Dune::InteriorEntity || it->partitionType()==Dune::BorderEntity)
        uniquemask[mapper.map(*it)] = true;
      else
        uniquemask[mapper.map(*it)] = false;
    OwnerExchange datahandle(grid,mapper,uniquemask);
    commwrapper.template communicate<OwnerExchange>(datahandle,Dune::InteriorBorder_InteriorBorder_Interface,
                                                    Dune::ForwardCommunication);

//     for (VIterator it = indexset.template begin<dim,Dune::All_Partition>(); it!=vendit; ++it)
//       std::cout << "rank=" << grid.comm().rank()
//                 << " index=" << mapper.map(*it)
//                 << " pos=" << it->geometry().corner(0)
//                 << " nonfront=" << nonfrontmask[mapper.map(*it)]
//                 << " unique=" << uniquemask[mapper.map(*it)]
//                 << std::endl;

  }

private:
  const Grid& grid;
  const IndexSet& indexset;
  CommWrapper commwrapper;
  Mapper mapper;
  std::vector<bool> nonfrontmask;
  std::vector<bool> uniquemask;
};


/** \brief Leaf version of grid based communication in the Schwarz method for overlapping grids

\tparam Grid The grid
*/
template<class Grid>
class LeafP1OverlappingSchwarzCommunication
  : public P1OverlappingSchwarzCommunication<Grid,typename Grid::template Codim<0>::LeafIndexSet,
                                             Dune::LeafCommunicate<Grid> >
{
public:
  /** \brief Constructor for a given grid
  */
  LeafP1OverlappingSchwarzCommunication (const Grid& grid)
    : P1OverlappingSchwarzCommunication<Grid,typename Grid::template Codim<0>::LeafIndexSet,
                                        Dune::LeafCommunicate<Grid> >
  (grid,grid.leafIndexSet(),Dune::LeafCommunicate<Grid>(grid))
  {}
};


/** \brief Level version of grid based communication in the Schwarz method for overlapping grids

\tparam Grid The grid
*/
template<class Grid>
class LevelP1OverlappingSchwarzCommunication
  : public P1OverlappingSchwarzCommunication<Grid,typename Grid::template Codim<0>::LevelIndexSet,
                                             Dune::LevelCommunicate<Grid> >
{
public:
  /** \brief Constructor for a given grid
  */
  LevelP1OverlappingSchwarzCommunication (const Grid& grid, int level)
    : P1OverlappingSchwarzCommunication<Grid,typename Grid::template Codim<0>::LevelIndexSet,
                                        Dune::LevelCommunicate<Grid> >
  (grid,grid.levelIndexSet(level),Dune::LevelCommunicate<Grid>(grid,level))
  {}
};





#endif
