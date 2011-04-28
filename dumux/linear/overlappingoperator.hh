// $Id$
/*****************************************************************************
 *   Copyright (C) 2011 by Andreas Lauser                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief A overlapping linear operator for use with ISTL
 */
#ifndef DUMUX_OVERLAPPING_OPERATOR_HH
#define DUMUX_OVERLAPPING_OPERATOR_HH

namespace Dumux {

// operator that resets result to zero at constrained DOFS
template<class OverlappingMatrix, class DomainVector, class RangeVector>
class OverlappingOperator : 
    public Dune::AssembledLinearOperator<OverlappingMatrix,
                                         DomainVector,
                                         RangeVector>
{
public:
    //! export types
    typedef OverlappingMatrix matrix_type;
    typedef DomainVector domain_type;
    typedef RangeVector range_type;
    typedef typename domain_type::field_type field_type;

    //redefine the category, that is the only difference
    enum {category=Dune::SolverCategory::overlapping};

    OverlappingOperator (const OverlappingMatrix& A)
        : A_(A)
    {}

    //! apply operator to x:  \f$ y = A(x) \f$
    virtual void apply (const DomainVector& x, RangeVector& y) const
    {
        A_.mv(x,y);

    }

    //! apply operator to x, scale and add:  \f$ y = y + \alpha A(x) \f$
    virtual void applyscaleadd (field_type alpha, const DomainVector& x, RangeVector& y) const
    {
        A_.usmv(alpha,x,y);
    }

    //! returns the matrix
    virtual const OverlappingMatrix& getmat () const
    {
        return A_;
    }

private:
    const OverlappingMatrix &A_;
};

} // namespace Dumux

#endif
