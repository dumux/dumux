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
 * \brief Dune style VTK functions
 */
#ifndef DUMUX_IO_VTK_FUNCTION_HH
#define DUMUX_IO_VTK_FUNCTION_HH

#include <string>
#include <memory>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/common.hh>
#include <dune/grid/io/file/vtk/function.hh>

#include <dumux/io/vtkprecision.hh>
#include <dumux/common/typetraits/typetraits.hh>

namespace Dumux::Vtk {

namespace Detail {

template<class Field>
double accessEntry(const Field& f, [[maybe_unused]] int mycomp, [[maybe_unused]] int i)
{
    if constexpr (IsIndexable<std::decay_t<decltype(std::declval<Field>()[0])>>{})
    {
        if constexpr (IsIndexable<std::decay_t<decltype(std::declval<Field>()[0][0])>>{})
            DUNE_THROW(Dune::InvalidStateException, "Invalid field type");
        else
            return f[i][mycomp];
    }
    else
        return f[i];
}

} // end namespace Detail

/*!
 * \ingroup InputOutput
 * \brief a VTK function that supports both scalar and vector values for each element
 *
 * \tparam GridView The Dune grid view type
 * \tparam Mapper The type used for mapping elements to indices in the field
 * \tparam F The field type (either vector of scalars or vectors)
 */
template <typename GridView, typename Mapper, typename F>
struct VectorP0VTKFunction : Dune::VTKFunction<GridView>
{
    enum { dim = GridView::dimension };
    using ctype = typename GridView::ctype;
    using Element = typename GridView::template Codim<0>::Entity;

public:

    //! return number of components
    int ncomps() const final { return nComps_; }

    //! get name
    std::string name() const final { return name_; }

    //! evaluate
    double evaluate(int mycomp, const Element& e, const Dune::FieldVector<ctype, dim>&) const final
    { return Detail::accessEntry(field_, mycomp, mapper_.index(e)); }

    //! get output precision for the field
    Dumux::Vtk::Precision precision() const final
    { return precision_; }

    //! Constructor
    VectorP0VTKFunction(const GridView& gridView,
                        const Mapper& mapper,
                        const F& field,
                        const std::string& name,
                        int nComps,
                        Dumux::Vtk::Precision precision = Dumux::Vtk::Precision::float32)
    : field_(field), name_(name), nComps_(nComps), mapper_(mapper), precision_(precision)
    {
        if (field.size()!=(unsigned int)(mapper.size()))
            DUNE_THROW(Dune::IOError, "VectorP0VTKFunction: size mismatch between field "
                                       << name << " (" << field.size() << ") and mapper (" << mapper.size() << ")");
    }

private:
    const F& field_;
    const std::string name_;
    int nComps_;
    const Mapper& mapper_;
    Dumux::Vtk::Precision precision_;
};

/*!
 * \ingroup InputOutput
 * \brief a VTK function that supports both scalar and vector values for each vertex
 *
 * \tparam GridView The Dune grid view type
 * \tparam Mapper The type used for mapping vertices to indices in the field
 * \tparam F The field type (either vector of scalars or vectors)
 */
template <typename GridView, typename Mapper, typename F>
struct VectorP1VTKFunction : Dune::VTKFunction<GridView>
{
    enum { dim = GridView::dimension };
    using ctype = typename GridView::ctype;
    using Element = typename GridView::template Codim<0>::Entity;

public:

    //! return number of components
    int ncomps() const final { return nComps_; }

    //! get name
    std::string name() const final { return name_; }

    //! evaluate
    double evaluate(int mycomp, const Element& e, const Dune::FieldVector<ctype, dim>& xi) const final
    {
        const unsigned int dim = Element::mydimension;
        const unsigned int nVertices = e.subEntities(dim);

        std::vector<Dune::FieldVector<ctype, 1>> cornerValues(nVertices);
        for (unsigned i = 0; i < nVertices; ++i)
            cornerValues[i] = Detail::accessEntry(field_, mycomp, mapper_.subIndex(e, i, dim));

        // (Ab)use the MultiLinearGeometry class to do multi-linear interpolation between scalars
        const Dune::MultiLinearGeometry<ctype, dim, 1> interpolation(e.type(), std::move(cornerValues));
        return interpolation.global(xi);
    }

    //! get output precision for the field
    Dumux::Vtk::Precision precision() const final
    { return precision_; }

    //! Constructor
    VectorP1VTKFunction(const GridView& gridView,
                        const Mapper& mapper,
                        const F& field,
                        const std::string& name,
                        int nComps,
                        Dumux::Vtk::Precision precision = Dumux::Vtk::Precision::float32)
    : field_(field), name_(name), nComps_(nComps), mapper_(mapper), precision_(precision)
    {
        if (field.size()!=(unsigned int)( mapper.size() ))
            DUNE_THROW(Dune::IOError, "VectorP1VTKFunction: size mismatch between field "
                                      << name << " (" << field.size() << ") and mapper (" << mapper.size() << ")");
    }
private:
    const F& field_;
    const std::string name_;
    int nComps_;
    const Mapper& mapper_;
    Dumux::Vtk::Precision precision_;
};

/*!
 * \ingroup InputOutput
 * \brief A VTK function that supports both scalar and vector values for each vertex.
 *        This expects the data to be organized by a two-dimensional field storing for
 *        each element the element-local nodal values. This can be used for the output
 *        of fields that are non-conforming due to e.g. constitutive relationships and
 *        where no extra degrees of freedom exist to disply the discontinuities.
 *
 * \tparam GridView The Dune grid view type
 * \tparam Mapper The type used for mapping elements to indices in the field
 * \tparam F The field type (either vector of scalars or vectors)
 */
template <typename GridView, typename Mapper, typename F>
struct VectorP1NonConformingVTKFunction : Dune::VTKFunction<GridView>
{
    enum { dim = GridView::dimension };
    using ctype = typename GridView::ctype;
    using Element = typename GridView::template Codim<0>::Entity;

public:

    //! return number of components
    int ncomps() const final { return nComps_; }

    //! get name
    std::string name() const final { return name_; }

    //! evaluate
    double evaluate(int mycomp, const Element& e, const Dune::FieldVector<ctype, dim>& xi) const final
    {
        const unsigned int dim = Element::mydimension;
        const unsigned int nVertices = e.subEntities(dim);

        std::vector<Dune::FieldVector<ctype, 1>> cornerValues(nVertices);
        for (unsigned i = 0; i < nVertices; ++i)
            cornerValues[i] = accessEntry_(mycomp, mapper_.index(e), i);

        // (Ab)use the MultiLinearGeometry class to do multi-linear interpolation between scalars
        const Dune::MultiLinearGeometry<ctype, dim, 1> interpolation(e.type(), std::move(cornerValues));
        return interpolation.global(xi);
    }

    //! get output precision for the field
    Dumux::Vtk::Precision precision() const final
    { return precision_; }

    //! Constructor
    VectorP1NonConformingVTKFunction(const GridView& gridView,
                                     const Mapper& mapper,
                                     const F& field,
                                     const std::string& name,
                                     int nComps,
                                     Dumux::Vtk::Precision precision = Dumux::Vtk::Precision::float32)
    : field_(field), name_(name), nComps_(nComps), mapper_(mapper), precision_(precision)
    {
        if (field.size()!=(unsigned int)(mapper.size()))
            DUNE_THROW(Dune::IOError, "VectorP1NonConformingVTKFunction: size mismatch between field "
                                      << name << " (" << field.size() << ") and mapper (" << mapper.size() << ")");
    }
private:
    //! access to the field
    double accessEntry_([[maybe_unused]] int mycomp, [[maybe_unused]] int eIdx, [[maybe_unused]] int cornerIdx) const
    {
        if constexpr (IsIndexable<decltype(std::declval<F>()[0])>{})
        {
            if constexpr (IsIndexable<decltype(std::declval<F>()[0][0])>{})
                return field_[eIdx][cornerIdx][mycomp];
            else
                return field_[eIdx][cornerIdx];
        }
        else
            DUNE_THROW(Dune::InvalidStateException, "Invalid field type");
    }

    const F& field_;
    const std::string name_;
    int nComps_;
    const Mapper& mapper_;
    Dumux::Vtk::Precision precision_;

};

/*!
 * \ingroup InputOutput
 * \brief struct that can hold any field that fulfills the VTKFunction interface
 *
 * \tparam GridView The Dune grid view type
 */
template<class GridView>
class Field
{
    enum { dim = GridView::dimension };
    using ctype = typename GridView::ctype;
    using Element = typename GridView::template Codim<0>::Entity;

public:
    // template constructor selects the right VTKFunction implementation
    template <typename F, class Mapper>
    Field(const GridView& gridView, const Mapper& mapper, F const& f,
          const std::string& name, int numComp = 1, int codim = 0,
          Dune::VTK::DataMode dm = Dune::VTK::conforming,
          Dumux::Vtk::Precision precision = Dumux::Vtk::Precision::float32)
    : codim_(codim)
    {
        if (codim == GridView::dimension)
        {
            if (dm == Dune::VTK::conforming)
                field_ = std::make_shared< VectorP1VTKFunction<GridView, Mapper, F> >(gridView, mapper, f, name, numComp, precision);
            else
                field_ = std::make_shared< VectorP1NonConformingVTKFunction<GridView, Mapper, F> >(gridView, mapper, f, name, numComp, precision);
        }
        else if (codim == 0)
            field_ = std::make_shared< VectorP0VTKFunction<GridView, Mapper, F> >(gridView, mapper, f, name, numComp, precision);
        else
            DUNE_THROW(Dune::NotImplemented, "Only element or vertex quantities allowed.");
    }

    //! return the name of this field
    std::string name () const { return field_->name(); }

    //! return the number of components of this field
    int ncomps() const { return field_->ncomps(); }

    //! return the precision of this field
    Dumux::Vtk::Precision precision() const
    { return field_->precision(); }

    //! codimension of the entities on which the field values live
    int codim() const { return codim_; }

    //! element-local evaluation of the field
    double evaluate(int mycomp,
                    const Element &element,
                    const Dune::FieldVector< ctype, dim > &xi) const
    { return field_->evaluate(mycomp, element, xi); }

    //! returns the underlying vtk function
    std::shared_ptr<const Dune::VTKFunction<GridView>> get() const
    { return field_; }

private:
    int codim_;
    // can point to anything fulfilling the VTKFunction interface
    std::shared_ptr<Dune::VTKFunction<GridView>> field_;
};

} // end namespace Dumux::Vtk

#endif
