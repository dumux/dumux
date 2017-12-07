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
 * \ingroup Properties
 * \ingroup ImplicitProperties
 * \ingroup MultidomainModel
 * \brief Specify default properties required in the subdomains of dune-multidomain
 */
#ifndef DUMUX_SUBDOMAIN_PROPERTY_DEFAULTS_HH
#define DUMUX_SUBDOMAIN_PROPERTY_DEFAULTS_HH

#include <numeric>
#include <dune/grid/multidomaingrid.hh>
#include <dune/pdelab/backend/istlvectorbackend.hh>
#include <dune/pdelab/backend/istlmatrixbackend.hh>
#include <dune/pdelab/finiteelementmap/qkfem.hh>
#include <dune/pdelab/constraints/conforming.hh>

#include "subdomainproperties.hh"
#include "properties.hh"
#include "localoperator.hh"
#include "boxcouplinglocalresidual.hh"

namespace Dumux
{

namespace Properties
{

// Specifies the grid type for the subdomains
SET_PROP(SubDomain, Grid)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, MultiDomainTypeTag) MultiDomain;
    typedef typename GET_PROP_TYPE(MultiDomain, Grid) HostGrid;
    enum { maxSubDomains = GET_PROP_VALUE(MultiDomain, MaxSubDomains) };
    typedef typename Dune::mdgrid::FewSubDomainsTraits<HostGrid::dimension,maxSubDomains> MDGridTraits;
    typedef typename Dune::MultiDomainGrid<HostGrid, MDGridTraits> Grid;
public:
    typedef typename Grid::SubDomainGrid type;
};

// set the default BaseLocalResidual to BoxCouplingLocalResidual
SET_TYPE_PROP(SubDomain, BaseLocalResidual, BoxCouplingLocalResidual<TypeTag>);

// set the local operator used for submodels
SET_TYPE_PROP(SubDomain, LocalOperator,
              PDELab::MultiDomainLocalOperator<TypeTag>);

// use the time manager for the coupled problem in the sub problems
SET_PROP(SubDomain, TimeManager)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, MultiDomainTypeTag) MultiDomainTypeTag;
public:
    typedef typename GET_PROP_TYPE(MultiDomainTypeTag, TimeManager) type;
};

// set the constraints for the sub-models
SET_TYPE_PROP(SubDomain, Constraints, Dune::PDELab::NoConstraints);

// set the grid functions space for the sub-models
SET_PROP(SubDomain, ScalarGridFunctionSpace)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, LocalFEMSpace) FEM;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Constraints) Constraints;
    enum{numEq = GET_PROP_VALUE(TypeTag, NumEq)};
public:
    typedef Dune::PDELab::GridFunctionSpace<GridView, FEM, Constraints,
        Dune::PDELab::ISTLVectorBackend<> > type;
};

// set the grid functions space for the sub-models
SET_PROP(SubDomain, GridFunctionSpace)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, ScalarGridFunctionSpace) ScalarGridFunctionSpace;
    enum{numEq = GET_PROP_VALUE(TypeTag, NumEq)};
    typedef typename Dune::PDELab::EntityBlockedOrderingTag OrderingTag;
    typedef typename Dune::PDELab::ISTLVectorBackend<> VBE;
public:
    typedef Dune::PDELab::PowerGridFunctionSpace<ScalarGridFunctionSpace, numEq, VBE, OrderingTag> type;
};

// use the local FEM space associated with cubes by default
SET_PROP(SubDomain, LocalFEMSpace)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    enum{dim = GridView::dimension};
public:
    typedef Dune::PDELab::QkLocalFiniteElementMap<GridView,Scalar,Scalar,1>  type;
};

SET_PROP(SubDomain, ParameterTree)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, MultiDomainTypeTag) MultiDomainTypeTag;
    typedef typename GET_PROP(MultiDomainTypeTag, ParameterTree) ParameterTree;
public:
    typedef typename ParameterTree::type type;

    static type &tree()
    { return ParameterTree::tree(); }

    static type &compileTimeParams()
    { return ParameterTree::compileTimeParams(); }

    static type &runTimeParams()
    { return ParameterTree::runTimeParams(); }

    static type &deprecatedRunTimeParams()
    { return ParameterTree::deprecatedRunTimeParams(); }

    static type &unusedNewRunTimeParams()
    { return ParameterTree::unusedNewRunTimeParams(); }
};

} // namespace Properties
} // namespace Dumux

#endif // DUMUX_SUBDOMAIN_PROPERTY_DEFAULTS_HH
