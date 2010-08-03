/*
 * \file
 *
 * \brief This is basically a copy of PDELab's ISTLVectorBackend which
 *        allows to set the global function space later
 */
#ifndef DUMUX_BOX_SOLUTION_VECTOR_HH
#define DUMUX_BOX_SOLUTION_VECTOR_HH

#include<vector>

#include<dune/common/fvector.hh>
#include<dune/istl/bvector.hh>

namespace Dumux {
namespace Properties {
NEW_PROP_TAG(NumEq);
NEW_PROP_TAG(Model);
};

namespace PDELab {

//! ISTL backend for FunctionSpace
template<class TypeTag>
class BoxISTLVectorBackend
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Model)) Model;

public:
    //! \brief export the block size
    static const int BlockSize = GET_PROP_VALUE(TypeTag, PTAG(NumEq));

    //! container construction
    template<typename T, typename E>
    class VectorContainer : public Dune::BlockVector< Dune::FieldVector<E,BlockSize> >
    {
    public:
        typedef E ElementType;
        typedef Dune::BlockVector< Dune::FieldVector<E,BlockSize> > BaseT;
        typedef BoxISTLVectorBackend<TypeTag> Backend;

        VectorContainer()
        {}
        VectorContainer(const Model &model) : BaseT(model.numDofs()) {}
        VectorContainer(const Model &model, E val) : BaseT(model.numDofs(), val) {}
        VectorContainer (const T& t_) : BaseT(t_.globalSize()/BlockSize)
        {}
        VectorContainer (const T& t_, const E& e) : BaseT(t_.globalSize()/BlockSize)
        {
            BaseT::operator=(e);
        }
        VectorContainer& operator= (const E& e)
        {
            BaseT::operator=(e);
            return *this;
        }

        // for debugging and AMG access
        BaseT& base ()
        {
            return *this;
        }

        const BaseT& base () const
        {
            return *this;
        }

        template<typename X>
        void std_copy_to (std::vector<X>& x) const
        {
            size_t n = this->size()*BlockSize;
            x.resize(n);
            for (size_t i=0; i<n; i++)
                x[i] = (*this)[i/BlockSize][i%BlockSize];
        }

        template<typename X>
        void std_copy_from (const std::vector<X>& x)
        {
            size_t n = this->size()*BlockSize;
            x.resize(n);
            for (size_t i=0; i<n; i++)
                (*this)[i/BlockSize][i%BlockSize] = x[i];
        }
    };

    // extract type of container element
    template<class C>
    struct Value
    {
        typedef typename C::field_type Type;
    };

    //! The size type
    typedef typename Dune::BlockVector< Dune::FieldVector<float,BlockSize> >::size_type size_type;

    // get const_reference to container element
    // note: this method does not depend on T!
    template<typename C>
    static const typename C::field_type& access (const C& c, size_type i)
    {
        return c[i/BlockSize][i%BlockSize];
    }

    // get non const_reference to container element
    // note: this method does not depend on T!
    template<typename C>
    static typename C::field_type& access (C& c, size_type i)
    {
        return c[i/BlockSize][i%BlockSize];
    }
};

} // namespace PDELab
} // namespace Dune

#endif
