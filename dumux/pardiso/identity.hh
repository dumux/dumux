// $Id:$
/*****************************************************************************
 *   Copyright (C) <YEARS> by <ADD_AUTHOR_HERE>                              *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/

#ifndef DUNE_IDENTITY_HH
#define DUNE_IDENTITY_HH

#include <dune/istl/preconditioners.hh>


namespace Dune {


/*! \brief The sequential Pardiso preconditioner.

  Put the Pardiso direct solver into the preconditioner framework.
*/
template<class M, class X, class Y>
class SeqIdentity : public Preconditioner<X,Y> {
public:
    //! \brief The matrix type the preconditioner is for.
    typedef M matrix_type;
    //! \brief The domain type of the preconditioner.
    typedef X domain_type;
    //! \brief The range type of the preconditioner.
    typedef Y range_type;
    //! \brief The field type of the preconditioner.
    typedef typename X::field_type field_type;

    typedef typename M::RowIterator RowIterator;
    typedef typename M::ColIterator ColIterator;

    // define the category
    enum {
        //! \brief The category the preconditioner is part of
        category=SolverCategory::sequential
    };

    /*! \brief Constructor.

      Constructor gets all parameters to operate the prec.
      \param A The matrix to operate on.
      \param n The number of iterations to perform.
      \param w The relaxation factor.
    */
    SeqIdentity (const M& A)
        : A_(A)
    {    }

    /*!
      \brief Prepare the preconditioner.

      \copydoc Preconditioner::pre(X&,Y&)
    */
    virtual void pre (X& x, Y& b) {}

    /*!
      \brief Apply the preconditioner.

      \copydoc Preconditioner::apply(X&,const Y&)
    */
    virtual void apply (X& v, const Y& d)
    {
        v = d;
        v *= 1;
    }

    /*!
      \brief Clean up.

      \copydoc Preconditioner::post(X&)
    */
    virtual void post (X& x) {}

    ~SeqIdentity()
    { }

private:
    M A_; //!< The matrix we operate on.
};

}






#endif

