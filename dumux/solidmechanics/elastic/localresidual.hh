// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Elastic
 * \brief Element-wise calculation of the local residual for problems
 *        using the elastic model considering linear elasticity.
 */
#ifndef DUMUX_SOLIDMECHANICS_ELASTIC_LOCAL_RESIDUAL_HH
#define DUMUX_SOLIDMECHANICS_ELASTIC_LOCAL_RESIDUAL_HH

#include <ranges>

#include <dune/common/fmatrix.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/typetraits/localdofs_.hh>
#include <dumux/discretization/defaultlocaloperator.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>

namespace Dumux {

/*!
 * \ingroup Elastic
 * \brief Element-wise calculation of the local residual for problems
 *        using the elastic model considering linear elasticity.
 *
 * Two assembly interfaces are supported, chosen by the grid variables. With
 * `FVGridVariables` (`FVAssembler`), the residual uses the volume variables, the stress type and
 * the problem's `neumann`/`dirichlet` interface. With `Experimental::GridVariables`
 * (`Experimental::Assembler`), it is written in terms of local degrees of freedom: degrees of
 * freedom with a control volume are tested with its indicator function, all others with their
 * basis function. The problem then provides the body force as `source(elemDisc, elemVars, ipData)`,
 * boundary terms as `boundaryFlux(elemDisc, elemVars, ipData)` and Dirichlet conditions as
 * constraints, and the spatial parameters provide `lameParams(elemDisc, ipData)`.
 */
template<class TypeTag>
class ElasticLocalResidual
: public DiscretizationDefaultLocalOperator<TypeTag>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using ParentType = DiscretizationDefaultLocalOperator<TypeTag>;

    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;

    // class assembling the stress tensor
    using StressType = GetPropType<TypeTag, Properties::StressType>;

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Tensor = Dune::FieldMatrix<Scalar, dim, dimWorld>;

public:
    using ParentType::ParentType;

    /*!
     * \brief For the elastic model the storage term is zero since
     *        we neglect inertia forces.
     *
     * \param problem The problem
     * \param scv The sub control volume
     * \param volVars The current or previous volVars
     */
    template<class SubControlVolume, class VolumeVariables>
    NumEqVector computeStorage(const Problem& problem,
                               const SubControlVolume& scv,
                               const VolumeVariables& volVars) const
    {
        return NumEqVector(0.0);
    }

    /*!
     * \brief Evaluates the force in all grid directions acting
     *        on a face of a sub-control volume.
     *
     * \param problem The problem
     * \param element The current element.
     * \param fvGeometry The finite-volume geometry
     * \param elemVolVars The volume variables of the current element
     * \param scvf The sub control volume face to compute the flux on
     * \param elemFluxVarsCache The cache related to flux computation
     */
    template<class ElementVolumeVariables, class SubControlVolumeFace, class ElementFluxVariablesCache>
    NumEqVector computeFlux(const Problem& problem,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolumeFace& scvf,
                            const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        // obtain force on the face from stress type
        return StressType::force(problem, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);
    }

    /*!
     * \brief Calculate the source term of the equation
     *
     * \param problem The problem to solve
     * \param element The DUNE Codim<0> entity for which the residual
     *                ought to be calculated
     * \param fvGeometry The finite-volume geometry of the element
     * \param elemVolVars The volume variables associated with the element stencil
     * \param scv The sub-control volume over which we integrate the source term
     * \note This is the default implementation for geomechanical models adding to
     *       the user defined sources the source stemming from the gravitational acceleration.
     *
     */
    template<class ElementVolumeVariables, class SubControlVolume>
    NumEqVector computeSource(const Problem& problem,
                              const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const ElementVolumeVariables& elemVolVars,
                              const SubControlVolume &scv) const
    {
        NumEqVector source(0.0);

        // add contributions from volume flux sources
        source += problem.source(element, fvGeometry, elemVolVars, scv);

        // add contribution from possible point sources
        source += problem.scvPointSources(element, fvGeometry, elemVolVars, scv);

        // maybe add gravitational acceleration
        static const bool gravity = getParamFromGroup<bool>(problem.paramGroup(), "Problem.EnableGravity");
        if (gravity)
        {
            const auto& g = problem.spatialParams().gravity(scv.center());
            for (int dir = 0; dir < GridView::dimensionworld; ++dir)
                source[Indices::momentum(dir)] += elemVolVars[scv].solidDensity()*g[dir];
        }

        return source;
    }
    //! The storage term is zero for quasi-static elasticity (assembly with local degrees of freedom)
    template<class ElementDiscretization, class ElementVariables, class Scv>
    NumEqVector storageIntegral(const ElementDiscretization&, const ElementVariables&, const Scv&, bool) const
    { return NumEqVector(0.0); }

    //! The body force \f$ \int \mathbf{f} \f$ over a control volume (assembly with local degrees of freedom)
    template<class ElementDiscretization, class ElementVariables, class Scv>
    NumEqVector sourceIntegral(const ElementDiscretization& elemDisc,
                               const ElementVariables& elemVars,
                               const Scv& scv) const
    {
        const auto& problem = this->asImp().problem();
        NumEqVector source(0.0);
        for (const auto& qpData : CVFE::quadratureRule(elemDisc, scv))
            source.axpy(qpData.weight(), problem.source(elemDisc, elemVars, qpData.ipData()));
        source *= elemVars[scv].extrusionFactor();
        return source;
    }

    //! The stress flux \f$ -\int \boldsymbol{\sigma}\mathbf{n} \f$ over a control volume face (assembly with local degrees of freedom)
    template<class ElementDiscretization, class ElementVariables, class Scvf>
    NumEqVector fluxIntegral(const ElementDiscretization& elemDisc,
                             const ElementVariables& elemVars,
                             const Scvf& scvf) const
    {
        const auto& problem = this->asImp().problem();
        NumEqVector flux(0.0);
        for (const auto& qpData : CVFE::quadratureRule(elemDisc, scvf))
        {
            const auto& ipData = qpData.ipData();
            const auto sigma = stressTensor_(problem, elemDisc, elemVars, ipData);
            for (int dir = 0; dir < dim; ++dir)
                flux[Indices::momentum(dir)] -= qpData.weight()*(sigma[dir]*scvf.unitOuterNormal());
        }
        flux *= elemVars[elemDisc.scv(scvf.insideScvIdx())].extrusionFactor();
        return flux;
    }

    /*!
     * \brief The Galerkin terms \f$ \int \boldsymbol{\sigma} : \nabla \mathbf{w} - \mathbf{f} \cdot \mathbf{w} \f$
     *        of the degrees of freedom without control volume (assembly with local degrees of freedom)
     */
    template<class ElementResidualVector, class ElementDiscretization, class ElementVariables>
    void addToElementFluxAndSourceResidual(ElementResidualVector& residual,
                                           const Problem& problem,
                                           const Element& element,
                                           const ElementDiscretization& elemDisc,
                                           const ElementVariables& elemVars) const
    {
        if constexpr (Detail::LocalDofs::hasNonCVLocalDofsInterface<ElementDiscretization>())
        {
            if (std::ranges::empty(nonCVLocalDofs(elemDisc)))
                return;

            for (const auto& qpData : CVFE::quadratureRule(elemDisc, element))
            {
                const auto& ipData = qpData.ipData();
                const auto& ipCache = cache(elemVars, ipData);
                const auto sigma = stressTensor_(problem, elemDisc, elemVars, ipData);
                const auto bodyForce = problem.source(elemDisc, elemVars, ipData);
                const auto& shapeValues = ipCache.shapeValues();
                for (const auto& localDof : nonCVLocalDofs(elemDisc))
                {
                    const auto idx = localDof.index();
                    for (int dir = 0; dir < dim; ++dir)
                    {
                        const auto eqIdx = Indices::momentum(dir);
                        residual[idx][eqIdx] += qpData.weight()*(sigma[dir]*ipCache.gradN(idx)
                                                                 - shapeValues[idx][0]*bodyForce[eqIdx]);
                    }
                }
            }
        }
    }

private:
    //! The stress \f$ \boldsymbol{\sigma} = \lambda \operatorname{tr}(\boldsymbol{\varepsilon}) \mathbf{I} + 2\mu \boldsymbol{\varepsilon} \f$ at an interpolation point
    template<class ElementDiscretization, class ElementVariables, class IpData>
    Tensor stressTensor_(const Problem& problem,
                         const ElementDiscretization& elemDisc,
                         const ElementVariables& elemVars,
                         const IpData& ipData) const
    {
        const auto& ipCache = cache(elemVars, ipData);
        Tensor gradU(0.0);
        for (const auto& localDof : localDofs(elemDisc))
            for (int dir = 0; dir < dim; ++dir)
                gradU[dir].axpy(elemVars[localDof].displacement(dir), ipCache.gradN(localDof.index()));

        const auto& lame = problem.spatialParams().lameParams(elemDisc, ipData);
        Tensor sigma(0.0);
        Scalar trace = 0.0;
        for (int i = 0; i < dim; ++i)
            trace += gradU[i][i];
        for (int i = 0; i < dim; ++i)
        {
            for (int j = 0; j < dimWorld; ++j)
                sigma[i][j] = lame.mu()*(gradU[i][j] + gradU[j][i]);
            sigma[i][i] += lame.lambda()*trace;
        }
        return sigma;
    }
};

} // end namespace Dumux

#endif
