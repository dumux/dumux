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
 * \brief The stress variables
 */
#ifndef DUMUX_GEOMECHANICS_IMPLICIT_STRESSVARIABLES_HH
#define DUMUX_GEOMECHANICS_IMPLICIT_STRESSVARIABLES_HH

#include <dumux/implicit/properties.hh>
#include <dumux/implicit/fluxvariablesbase.hh>

namespace Dumux
{

namespace Properties
{
NEW_PROP_TAG(MechanicalLawType);
}

/*!
 * \ingroup ImplicitModel
 * \brief the flux variables class
 *        specializations are provided for combinations of physical processes
 */
template<class TypeTag, bool enableAdvection, bool enableMolecularDiffusion, bool enableEnergyBalance>
class StressVariables {};

/*!
 * \ingroup ImplicitModel
 * \brief Base class for the flux variables
 *        Actual flux variables inherit from this class
 */
template<class TypeTag>
class StressVariables<TypeTag, false, false, false> : public FluxVariablesBase<TypeTag, StressVariables<TypeTag, false, false, false>>
{
    using ParentType = FluxVariablesBase<TypeTag, StressVariables<TypeTag, false, false, false>>;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using IndexType = typename GridView::IndexSet::IndexType;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Stencil = std::vector<IndexType>;
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using MechanicalLawType = typename GET_PROP_TYPE(TypeTag, MechanicalLawType);

    enum{ enableFluxVarsCache = GET_PROP_VALUE(TypeTag, EnableFluxVariablesCache) };
    enum { dim = GridView::dimension} ;

    typedef Dune::FieldMatrix<Scalar, dim, dim> DimMatrix;
    typedef Dune::FieldVector<Scalar, dim> DimVector;

public:

    void initAndComputeFluxes(const Problem& problem,
                              const Element& element,
                              const SubControlVolumeFace &scvFace)
    {
        ParentType::init(problem, element, scvFace);
    }

    Stencil computeStencil(const Problem& problem, const SubControlVolumeFace& scvFace)
    { return MechanicalLawType::stencil(problem, scvFace); }

    DimVector stressVector()
    { return MechanicalLawType::stressVector(this->problem(), this->scvFace()); }

    DimMatrix stressTensor()
    { return MechanicalLawType::stressTensor(this->problem(), this->scvFace()); }
};

} // end namespace

#endif
