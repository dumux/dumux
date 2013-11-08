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
#ifndef DUMUX_SUBDOMAIN_PROPERTY_DEFAULTS_HH
#define DUMUX_SUBDOMAIN_PROPERTY_DEFAULTS_HH

#include <dune/grid/multidomaingrid.hh>
#include <dumux/implicit/pdelab/pdelabadapter.hh>
#include <dumux/implicit/common/boxcouplinglocalresidual.hh>
#include <dumux/modelcoupling/common/multidomainboxlocaloperator.hh>
#include "coupledproperties.hh"

/*!
 * \file
 * \brief Specify default properties required in the subdomains of dune-multidomain
 */
namespace Dumux
{
namespace Properties
{
//! The type tag for problems which use dune-multidomain
NEW_TYPE_TAG(SubDomain, INHERITS_FROM(BoxPDELab, CoupledSubProblem));

//////////////////////////////////////
// Set property values
//////////////////////////////////////

// Specifies the grid type for the subdomains
SET_PROP(SubDomain, Grid)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, CoupledProblemTypeTag) CoupledTypeTag;
    typedef typename GET_PROP_TYPE(CoupledTypeTag, Grid) HostGrid;
    typedef typename Dune::mdgrid::FewSubDomainsTraits<HostGrid::dimension,4> MDGridTraits;
    typedef typename Dune::MultiDomainGrid<HostGrid, MDGridTraits> Grid;

public:
    typedef typename Grid::SubDomainGrid type;
};

// Set the default BaseLocalResidual to BoxCouplingLocalResidual
SET_TYPE_PROP(SubDomain, BaseLocalResidual, BoxCouplingLocalResidual<TypeTag>);

// set the local operator used for submodels
SET_TYPE_PROP(SubDomain, LocalOperator,
              Dumux::PDELab::MultiDomainBoxLocalOperator<TypeTag>);


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
        Dune::PDELab::ISTLVectorBackend<1> > type;
};

// \}

}
}

#endif
