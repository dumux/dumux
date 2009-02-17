// A DataHandle class to exchange entries of a vector
template<class M, class V> // mapper type and vector type
class VectorExchange 
  : public Dune::CommDataHandleIF<VectorExchange<M,V>,
				  typename V::value_type>
{
public:
  //! export type of data for message buffer
  typedef typename V::value_type DataType;

  //! returns true if data for this codim should be communicated
  bool contains (int dim, int codim) const
  {
	return (codim==0);
  }

  //! returns true if size per entity of given dim and codim is a constant
  bool fixedsize (int dim, int codim) const
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
	buff.write(c[mapper.map(e)]);
  }

  /*! unpack data from message buffer to user

  n is the number of objects sent by the sender
  */
  template<class MessageBuffer, class EntityType>
  void scatter (MessageBuffer& buff, const EntityType& e, size_t n)
  {
	DataType x;
	buff.read(x);
	c[mapper.map(e)]=x;
  }

  //! constructor
  VectorExchange (const M& mapper_, V& c_) 
	: mapper(mapper_), c(c_)
  {}
 
private:
  const M& mapper;
  V& c;
};
