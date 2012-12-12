// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*****************************************************************************
 *   This file is based on the file of DUNE-PDELab with the same name,       *
 *   http://www.dune-project.org/pdelab/.                                    *
 *   Some changes are made such that the backends work for the vectors and   *
 *   matrices used in Dumux.                                                 *
 *****************************************************************************/
#ifndef DUNE_ISTLVECTORBACKEND_HH
#define DUNE_ISTLVECTORBACKEND_HH

#include<vector>

#include <dune/common/fvector.hh>
#include <dune/istl/bvector.hh>
#include <dune/pdelab/backend/backendselector.hh>
#include <dune/pdelab/backend/istlmatrixbackend.hh>

namespace Dumux {

    template<int> class ISTLVectorBackend;

    template<typename T, typename E, int BLOCKSIZE>
    class ISTLBlockVectorContainer
    {
    public:
      typedef E ElementType;
      typedef Dune::BlockVector< Dune::FieldVector<E,BLOCKSIZE> > ContainerType;
      typedef ContainerType BaseT;
      typedef typename ContainerType::field_type field_type;
      typedef typename ContainerType::iterator iterator;
      typedef typename ContainerType::const_iterator const_iterator;
      typedef typename ContainerType::block_type block_type;
      typedef typename ContainerType::size_type size_type;
      typedef ISTLVectorBackend<BLOCKSIZE> Backend;

      ISTLBlockVectorContainer ()
      {}
      ISTLBlockVectorContainer (const T& t_) : container(t_.globalSize()/BLOCKSIZE)
      {}
      ISTLBlockVectorContainer (const T& t_, const E& e) : container(t_.globalSize()/BLOCKSIZE)
      {
        container=e;
      }

      size_type N() const
      {
        return container.N();
      }


      ISTLBlockVectorContainer& operator= (const E& e)
      {
        container=e;
        return *this;
      }

      ISTLBlockVectorContainer& operator*= (const E& e)
      {
        container*=e;
        return *this;
      }


      ISTLBlockVectorContainer& operator+= (const E& e)
      {
        container+=e;
        return *this;
      }

      ISTLBlockVectorContainer& operator+= (const ISTLBlockVectorContainer& e)
      {
        container+=e;
        return *this;
      }

      ISTLBlockVectorContainer& operator-= (const ISTLBlockVectorContainer& e)
      {
        container-=e;
        return *this;
      }

      block_type& operator[](std::size_t i)
      {
        return container[i];
      }

      const block_type& operator[](std::size_t i) const
      {
        return container[i];
      }

      typename Dune::template FieldTraits<E>::real_type two_norm() const
      {
        return container.two_norm();
      }

      typename Dune::template FieldTraits<E>::real_type two_norm2() const
      {
        return container.two_norm2();
      }

      typename Dune::template FieldTraits<E>::real_type one_norm() const
      {
        return container.one_norm();
      }

      typename Dune::template FieldTraits<E>::real_type infinity_norm() const
      {
        return container.infinity_norm();
      }

      E operator*(const ISTLBlockVectorContainer& y) const
      {
        return container*y.base();
      }

      E dot(const ISTLBlockVectorContainer& y) const
      {
        return container.dot(y.base());
      }

      ISTLBlockVectorContainer& axpy(const E& a, const ISTLBlockVectorContainer& y)
      {
        container.axpy(a, y);
        return *this;
      }

      // for debugging and AMG access
      ContainerType& base ()
      {
        return container;
      }

      const ContainerType& base () const
      {
        return container;
      }

      operator ContainerType&()
      {
        return container;
      }

      operator const ContainerType&() const
      {
        return container;
      }

      iterator begin()
      {
        return container.begin();
      }


      const_iterator begin() const
      {
        return container.begin();
      }

      iterator end()
      {
        return container.end();
      }


      const_iterator end() const
      {
        return container.end();
      }

      size_t flatsize() const
      {
        return container.size()*BLOCKSIZE;
      }

      size_t size() const
      {
        return container.size();
      }

      void resize(size_t n) 
      {
        container.resize(n);
      }

      template<typename X>
      void std_copy_to (std::vector<X>& x) const
      {
        size_t n = flatsize();
        x.resize(n);
        for (size_t i=0; i<n; i++)
          x[i] = container[i/BLOCKSIZE][i%BLOCKSIZE];
      }

      template<typename X>
      void std_copy_from (const std::vector<X>& x)
      {
        //test if x has the same size as the container
        assert (x.size() == flatsize());
        for (size_t i=0; i<flatsize(); i++)
          container[i/BLOCKSIZE][i%BLOCKSIZE] = x[i];
      }

    private:
      Dune::BlockVector< Dune::FieldVector<E,BLOCKSIZE> > container;
    };


    //! ISTL backend for FunctionSpace
    template<int BLOCKSIZE=1>
    class ISTLVectorBackend
    {
    public:
      enum{
        //! \brief export the block size
        BlockSize = BLOCKSIZE
      };

      //export Matrix Backend Type
      typedef ISTLBCRSMatrixBackend<BLOCKSIZE,BLOCKSIZE> MatrixBackend;

      //! container construction

      // extract type of container element
      template<class C>
      struct Value
      {
        typedef typename C::field_type Type;
      };

      //! The size type
      typedef typename Dune::BlockVector< Dune::FieldVector<float,BLOCKSIZE> >::size_type size_type;

      // get const_reference to container element
      // note: this method does not depend on T!
      template<typename C>
      static const typename C::field_type& access (const C& c, size_type i)
      {
        return c.base()[i/BLOCKSIZE][i%BLOCKSIZE];
      }

      // get non const_reference to container element
      // note: this method does not depend on T!
      template<typename C>
      static typename C::field_type& access (C& c, size_type i)
      {
        return c.base()[i/BLOCKSIZE][i%BLOCKSIZE];
      }
    };

    template<int BLOCKSIZE,typename T, typename E>
    struct BackendVectorSelectorHelper<ISTLVectorBackend<BLOCKSIZE>, T, E>
    {
      typedef ISTLBlockVectorContainer<T,E,BLOCKSIZE> Type;
    };



} // namespace Dumux

#endif
