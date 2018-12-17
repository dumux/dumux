// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup InputOutput
 * \brief Provides a vector-valued function using Dune::FieldVectors
 *        as elements.
 */
#ifndef VTK_NESTED_FUNCTION_HH
#define VTK_NESTED_FUNCTION_HH

#include <string>

#include <dune/common/fvector.hh>
#include <dune/istl/bvector.hh>
#include <dune/grid/io/file/vtk/function.hh>

namespace Dumux {

/*!
 * \ingroup InputOutput
 * \brief Provides a vector-valued function using Dune::FieldVectors
 *        as elements.
 */
template <class GridView, class Mapper, class Buffer>
class VtkNestedFunction : public Dune::VTKFunction<GridView>
{
    enum { dim = GridView::dimension };
    using ctype = typename GridView::ctype;
    using Element = typename GridView::template Codim<0>::Entity;
    using ReferenceElements = Dune::ReferenceElements<ctype, dim>;

public:
    VtkNestedFunction(std::string name,
                      const GridView &gridView,
                      const Mapper &mapper,
                      const Buffer &buf,
                      int codim,
                      int numComp)
        : name_(name)
        , gridView_(gridView)
        , mapper_(mapper)
        , buf_(buf)
        , codim_(codim)
        , numComp_(numComp)
    {
        assert(std::size_t(buf_.size()) == std::size_t(mapper_.size()));
    }

    virtual std::string name () const
    { return name_; }

    virtual int ncomps() const
    { return numComp_; }

    virtual double evaluate(int mycomp,
                            const Element &element,
                            const Dune::FieldVector< ctype, dim > &xi) const
    {
        int idx;
        if (codim_ == 0) {
            // cells. map element to the index
            idx = mapper_.index(element);
        }
        else if (codim_ == dim) {
            // find vertex which is closest to xi in local
            // coordinates. This code is based on Dune::P1VTKFunction
            double min=1e100;
            int imin=-1;
            Dune::GeometryType geomType = element.type();
            int n = element.subEntities(dim);

            for (int i=0; i < n; ++i)
            {
                Dune::FieldVector<ctype,dim> local =
                    ReferenceElements::general(geomType).position(i,dim);
                local -= xi;
                if (local.infinity_norm()<min)
                {
                    min = local.infinity_norm();
                    imin = i;
                }
            }

            // map vertex to an index
            idx = mapper_.subIndex(element, imin, codim_);
        }
        else
            DUNE_THROW(Dune::InvalidStateException,
                       "Only element and vertex based vector "
                       " fields are supported so far.");

        double val = buf_[idx][mycomp];
        using std::abs;
        if (abs(val) < std::numeric_limits<float>::min())
            val = 0;

        return val;
    }

private:
    const std::string name_;
    const GridView gridView_;
    const Mapper &mapper_;
    const Buffer &buf_;
    int codim_;
    int numComp_;
};

} // end namespace Dumux

#endif
