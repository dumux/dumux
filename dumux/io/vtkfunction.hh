// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
#ifndef VTK_FUNCTION_HH
#define VTK_FUNCTION_HH

#include <string>
#include <memory>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/common.hh>
#include <dune/grid/io/file/vtk/function.hh>

#include <dumux/common/typetraits/typetraits.hh>

namespace Dumux {
namespace Vtk {

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
    virtual int ncomps() const { return nComps_; }

    //! get name
    virtual std::string name() const { return name_; }

    //! evaluate
    virtual double evaluate(int mycomp, const Element& e, const Dune::FieldVector<ctype, dim>&) const
    { return accessChooser_(mycomp, mapper_.index(e), IsIndexable<decltype(field_[0])>()); }

    //! Constructor
    VectorP0VTKFunction(const GridView& gridView, const Mapper& mapper, const F& field, const std::string& name, int nComps)
    : field_(field), name_(name), nComps_(nComps), mapper_(mapper)
    {
        if (field.size()!=(unsigned int)(mapper.size()))
            DUNE_THROW(Dune::IOError, "VectorP0VTKFunction: size mismatch between field "
                                       << name << " (" << field.size() << ") and mapper (" << mapper.size() << ")");
    }

private:

    //! access for vectorial fields
    double accessChooser_(int mycomp, int i, std::true_type) const
    { return vectorFieldAccess_(mycomp, i, IsIndexable<decltype(field_[0][0])>()); }

    //! access for scalar fields
    double accessChooser_(int mycomp, int i, std::false_type) const { return field_[i]; }

    //! access to permissive vectorial fields
    double vectorFieldAccess_(int mycomp, int i, std::false_type) const { return field_[i][mycomp]; }

    //! if the field is indexable more than two times, throw error
    double vectorFieldAccess_(int mycomp, int i, std::true_type) const { DUNE_THROW(Dune::InvalidStateException, "Invalid field type"); }

    const F& field_;
    const std::string name_;
    int nComps_;
    const Mapper& mapper_;
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
    virtual int ncomps() const { return nComps_; }

    //! get name
    virtual std::string name() const { return name_; }

    //! evaluate
    virtual double evaluate(int mycomp, const Element& e, const Dune::FieldVector<ctype, dim>& xi) const
    {
        const unsigned int dim = Element::mydimension;
        const unsigned int nVertices = e.subEntities(dim);

        std::vector<Dune::FieldVector<ctype, 1>> cornerValues(nVertices);
        for (unsigned i = 0; i < nVertices; ++i)
            cornerValues[i] = accessChooser_(mycomp, mapper_.subIndex(e, i, dim), IsIndexable<decltype(field_[0])>());

        // (Ab)use the MultiLinearGeometry class to do multi-linear interpolation between scalars
        const Dune::MultiLinearGeometry<ctype, dim, 1> interpolation(e.type(), std::move(cornerValues));
        return interpolation.global(xi);
    }

    //! Constructor
    VectorP1VTKFunction(const GridView& gridView, const Mapper& mapper, const F& field, const std::string& name, int nComps)
    : field_(field), name_(name), nComps_(nComps), mapper_(mapper)
    {
        if (field.size()!=(unsigned int)( mapper.size() ))
            DUNE_THROW(Dune::IOError, "VectorP1VTKFunction: size mismatch between field "
                                      << name << " (" << field.size() << ") and mapper (" << mapper.size() << ")");
    }
private:

    //! access for vectorial fields
    double accessChooser_(int mycomp, int i, std::true_type) const
    { return vectorFieldAccess_(mycomp, i, IsIndexable<decltype(field_[0][0])>()); }

    //! access for scalar fields
    double accessChooser_(int mycomp, int i, std::false_type) const { return field_[i]; }

    //! access to permissive vectorial fields
    double vectorFieldAccess_(int mycomp, int i, std::false_type) const { return field_[i][mycomp]; }

    //! if the field is indexable more than two times, throw error
    double vectorFieldAccess_(int mycomp, int i, std::true_type) const { DUNE_THROW(Dune::InvalidStateException, "Invalid field type"); }

    const F& field_;
    const std::string name_;
    int nComps_;
    const Mapper& mapper_;
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
    virtual int ncomps() const { return nComps_; }

    //! get name
    virtual std::string name() const { return name_; }

    //! evaluate
    virtual double evaluate(int mycomp, const Element& e, const Dune::FieldVector<ctype, dim>& xi) const
    {
        const unsigned int dim = Element::mydimension;
        const unsigned int nVertices = e.subEntities(dim);

        std::vector<Dune::FieldVector<ctype, 1>> cornerValues(nVertices);
        for (unsigned i = 0; i < nVertices; ++i)
            cornerValues[i] = accessChooser_(mycomp, mapper_.index(e), i, IsIndexable<decltype(field_[0])>());

        // (Ab)use the MultiLinearGeometry class to do multi-linear interpolation between scalars
        const Dune::MultiLinearGeometry<ctype, dim, 1> interpolation(e.type(), std::move(cornerValues));
        return interpolation.global(xi);
    }

    //! Constructor
    VectorP1NonConformingVTKFunction(const GridView& gridView, const Mapper& mapper, const F& field, const std::string& name, int nComps)
    : field_(field), name_(name), nComps_(nComps), mapper_(mapper)
    {
        if (field.size()!=(unsigned int)(mapper.size()))
            DUNE_THROW(Dune::IOError, "VectorP1NonConformingVTKFunction: size mismatch between field "
                                      << name << " (" << field.size() << ") and mapper (" << mapper.size() << ")");
    }
private:

    //! access to the field
    double accessChooser_(int mycomp, int eIdx, int cornerIdx, std::true_type) const
    { return fieldAccess_(mycomp, eIdx, cornerIdx, IsIndexable<decltype(field_[0][0])>()); }

    //! fields have to be indexable at least twice
    double accessChooser_(int mycomp, int eIdx, int cornerIdx, std::false_type) const
    { DUNE_THROW(Dune::InvalidStateException, "Invalid field type"); }

    //! scalar field access
    double fieldAccess_(int mycomp, int eIdx, int cornerIdx, std::false_type) const
    { return field_[eIdx][cornerIdx]; }

    //! vectorial field access
    double fieldAccess_(int mycomp, int eIdx, int cornerIdx, std::true_type) const
    { return field_[eIdx][cornerIdx][mycomp]; }

    const F& field_;
    const std::string name_;
    int nComps_;
    const Mapper& mapper_;
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
          Dune::VTK::DataMode dm = Dune::VTK::conforming)
    : codim_(codim)
    {
        if (codim == GridView::dimension)
        {
            if (dm == Dune::VTK::conforming)
                field_ = std::make_shared< VectorP1VTKFunction<GridView, Mapper, F> >(gridView, mapper, f, name, numComp);
            else
                field_ = std::make_shared< VectorP1NonConformingVTKFunction<GridView, Mapper, F> >(gridView, mapper, f, name, numComp);
        }
        else if (codim == 0)
            field_ = std::make_shared< VectorP0VTKFunction<GridView, Mapper, F> >(gridView, mapper, f, name, numComp);
        else
            DUNE_THROW(Dune::NotImplemented, "Only element or vertex quantities allowed.");
    }

    virtual ~Field() {}

    //! return the name of this field
    virtual std::string name () const { return field_->name(); }

    //! return the number of components of this field
    virtual int ncomps() const { return field_->ncomps(); }

    //! codimension of the entities on which the field values live
    int codim() const { return codim_; }

    //! element-local evaluation of the field
    virtual double evaluate(int mycomp,
                            const Element &element,
                            const Dune::FieldVector< ctype, dim > &xi) const
    { return field_->evaluate(mycomp, element, xi); }

    //! returns the underlying vtk function
    const std::shared_ptr<Dune::VTKFunction<GridView>>& get() const
    { return field_; }

private:
    int codim_;
    // can point to anything fulfilling the VTKFunction interface
    std::shared_ptr<Dune::VTKFunction<GridView>> field_;
};

} // end namespace Vtk

} // end namespace Dumux

#endif
