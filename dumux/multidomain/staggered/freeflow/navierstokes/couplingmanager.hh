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
 * \ingroup MultiDomain
 * \brief The interface of the coupling manager for multi domain problems
 */

#ifndef DUMUX_MD_STAGGERED_COUPLING_MANAGER_HH
#define DUMUX_MD_STAGGERED_COUPLING_MANAGER_HH

#include <dumux/multidomain/couplingmanager.hh>
#include <dumux/assembly/numericepsilon.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/discretization/method.hh>

namespace Dumux {

/*!
 * \ingroup MultiDomain
 * \ingroup StaggeredDiscretization
 * \brief Base coupling manager for the staggered discretization.
 */
template<class MDTraits>
class StaggeredCouplingManager: virtual public CouplingManager<MDTraits>
{
    template<std::size_t id> using GridGeometry = typename MDTraits::template SubDomain<id>::GridGeometry;

public:

    static constexpr auto cellCenterIdx = GridGeometry<0>::cellCenterIdx();
    static constexpr auto faceIdx = GridGeometry<0>::faceIdx();

private:

    using ParentType = CouplingManager<MDTraits>;

    template<std::size_t id> using Problem =  typename MDTraits::template SubDomain<id>::Problem;
    template<std::size_t id> using GridView = typename GridGeometry<id>::GridView;

    template<std::size_t id> using FVElementGeometry = typename GridGeometry<id>::LocalView;
    template<std::size_t id> using SubControlVolume = typename FVElementGeometry<id>::SubControlVolume;
    template<std::size_t id> using SubControlVolumeFace = typename FVElementGeometry<id>::SubControlVolumeFace;
    template<std::size_t id> using Element = typename GridView<id>::template Codim<0>::Entity;

    template<std::size_t id> using GridVariables = typename MDTraits::template SubDomain<id>::GridVariables;
    template<std::size_t id> using GridVolumeVariables = typename MDTraits::template SubDomain<id>::GridVariables::GridVolumeVariables;
    template<std::size_t id> using ElementVolumeVariables = typename GridVolumeVariables<id>::LocalView;


    using GridIndexType = typename IndexTraits< GridView<0> >::GridIndex;
    using CouplingStencil = std::vector<GridIndexType>;

    using GridVariablesTuple = typename MDTraits::template TupleOfSharedPtr<GridVariables>;
    using ProblemTuple = typename MDTraits::template TupleOfSharedPtrConst<Problem>;
    using SolutionVector = typename MDTraits::SolutionVector;

    struct FaceCouplingContext
    {
        FVElementGeometry<cellCenterIdx> fvGeometry;
        ElementVolumeVariables<cellCenterIdx> elemVolVars;
    };


public:

    using ParentType::evalCouplingResidual;

    using Traits = MDTraits;

    //! Initialize the coupling manager
    void init(ProblemTuple&& problem,
              GridVariablesTuple&& gridVariables,
              const SolutionVector& curSol)
    {
        this->setSubProblems(problem);

        if (Dune::FloatCmp::ne(this->problem(cellCenterIdx).gravity(), this->problem(faceIdx).gravity()))
            DUNE_THROW(Dune::InvalidStateException, "Both models must use the same gravity vector");

        gridVariables_ = gridVariables;
        this->curSol() = curSol;
        // couplingData_ = std::make_shared<CouplingData>(*this);
        // TODO stencils
    }

    //! Initialize the coupling manager
    void init(ProblemTuple&& problem,
              GridVariablesTuple&& gridVariables,
              const SolutionVector& curSol,
              const SolutionVector& prevSol)
    {
        prevSol_ = &prevSol;
        init(std::forward<decltype(problem)>(problem), std::forward<decltype(gridVariables)>(gridVariables), curSol);
    }

    /*!
     * \copydoc Dumux::CouplingManager::updateCouplingContext()
     */
    template<std::size_t i, std::size_t j, class LocalAssemblerI, class PriVarsJ>
    void updateCouplingContext(Dune::index_constant<i> domainI,
                               const LocalAssemblerI& localAssemblerI,
                               Dune::index_constant<j> domainJ,
                               const std::size_t dofIdxGlobalJ,
                               const PriVarsJ& priVarsJ,
                               int pvIdxJ)
    {
        auto& curSol = this->curSol()[domainJ];

        // only proceed if the solution vector has actually been set in the base class
        // (which is not the case if the staggered model is not coupled to another domain)
        if (curSol.size() > 0)
            curSol[dofIdxGlobalJ][pvIdxJ] = priVarsJ[pvIdxJ];
    }

    /*!
     * \copydoc Dumux::CouplingManager::updateCouplingContext()
     */
    template<class LocalAssemblerI, class PriVarsJ>
    void updateCouplingContext(Dune::index_constant<faceIdx> domainI,
                               const LocalAssemblerI& localAssemblerI,
                               Dune::index_constant<cellCenterIdx> domainJ,
                               const std::size_t dofIdxGlobalJ,
                               const PriVarsJ& priVarsJ,
                               int pvIdxJ)
    {
        auto& curSol = this->curSol()[domainJ];

        // only proceed if the solution vector has actually been set in the base class
        // (which is not the case if the staggered model is not coupled to another domain)
        if (curSol.size() > 0)
            curSol[dofIdxGlobalJ][pvIdxJ] = priVarsJ[pvIdxJ];

        const auto& fvGeometry = faceCouplingContext_->fvGeometry;
        const auto& scv = fvGeometry.scv(dofIdxGlobalJ);
        auto& volVars = getVolVarAccess(cellCenterIdx, scv);
        using VolumeVariables = std::decay_t<decltype(volVars)>;
        using PrimaryVariables = typename VolumeVariables::PrimaryVariables;
        PrimaryVariables priVars = makePriVarsFromCellCenterPriVars<PrimaryVariables>(priVarsJ);
        auto elemSol = elementSolution<FVElementGeometry<cellCenterIdx>>(std::move(priVars));
        volVars.update(elemSol, this->problem(domainI), localAssemblerI.element(), scv);
    }

    /*!
     * \brief returns an iteratable container of all indices of degrees of freedom of domain j
     *        that couple with / influence the element residual of the given element of domain i
     *
     * \param domainI the domain index of domain i
     * \param elementI the coupled element of domain í
     * \param domainJ the domain index of domain j

     * \note  this is a specialization for getting the indices of the coupled staggered face dofs
     */
    const CouplingStencil& couplingStencil(Dune::index_constant<cellCenterIdx> domainI,
                                           const Element<cellCenterIdx>& elementI,
                                           Dune::index_constant<faceIdx> domainJ) const
    {
        const auto& connectivityMap = this->problem(domainI).gridGeometry().connectivityMap();
        const auto eIdx = this->problem(domainI).gridGeometry().elementMapper().index(elementI);
        return connectivityMap(domainI, domainJ, eIdx);
    }

    /*!
     * \brief returns an iteratable container of all indices of degrees of freedom of domain j
     *        that couple with / influence the residual of the given subcontrolvolume face of domain i
     *
     * \param domainI the domain index of domain i
     * \param scvfI the coupled subcontrolvolume face of domain í
     * \param domainJ the domain index of domain j

     * \note  This function has to be implemented by all coupling managers for all combinations of i and j
     */
    template<std::size_t i, std::size_t j>
    const CouplingStencil couplingStencil(Dune::index_constant<i> domainI,
                                          const SubControlVolumeFace<i>& scvfI,
                                          Dune::index_constant<j> domainJ) const
    {
        static_assert(i != j, "Domain i cannot be coupled to itself!");
        static_assert(AlwaysFalse<Dune::index_constant<i>>::value,
                      "The coupling manager does not implement the couplingStencil() function" );

        return CouplingStencil(); // supress compiler warning of function having no return statement
    }

    /*!
     * \brief returns an iteratable container of all indices of degrees of freedom of domain j
     *        that couple with / influence the residual of the given subcontrolvolume face of domain i
     *
     * \param domainI the domain index of domain i
     * \param scvfI the coupled subcontrolvolume face of domain í
     * \param domainJ the domain index of domain j

     * \note this is a specialization for getting the indices of the coupled cellcentered dofs
     */
    const CouplingStencil& couplingStencil(Dune::index_constant<faceIdx> domainI,
                                           const SubControlVolume<faceIdx>& scvI,
                                           Dune::index_constant<cellCenterIdx> domainJ) const
    {
        const auto& connectivityMap = this->problem(domainI).gridGeometry().connectivityMap();
        return connectivityMap(domainI, domainJ, scvI.correspondingCellCenterScvfIndex());
    }

    /*!
     * \copydoc Dumux::CouplingManager::evalCouplingResidual()
     *
     * \note this is a specialization for calculating the coupled residual for cellcentered dofs
     *       w.r.t. staggered face dof changes
     */
    template<class LocalAssemblerI, std::size_t j>
    decltype(auto) evalCouplingResidual(Dune::index_constant<cellCenterIdx> domainI,
                                        const LocalAssemblerI& localAssemblerI,
                                        Dune::index_constant<j> domainJ,
                                        std::size_t dofIdxGlobalJ) const
    {
        static_assert(domainI != domainJ, "Domain i cannot be coupled to itself!");
        return localAssemblerI.evalLocalResidualForCellCenter();
    }

     /*!
      * \ingroup MultiDomain
      * \brief evaluates the face residual of a coupled face of domain i which depends on the variables
      *        at the degree of freedom with index dofIdxGlobalJ of domain j
      *
      * \param domainI the domain index of domain i
      * \param scvfI the subcontrol volume face whose residual shall be evaluated of domain i
      * \param localAssemblerI the local assembler assembling the element residual of an element of domain i
      * \param domainJ the domain index of domain j
      * \param dofIdxGlobalJ the index of the degree of freedom of domain j which has an influence on the element residual of domain i
      *
      * \note  the default implementation evaluates the complete face residual
      *        if only certain terms of the residual of the residual are coupled
      *        to dof with index dofIdxGlobalJ the function can be overloaded in the coupling manager
      * \return the face residual
      */
    template<class LocalAssemblerI, std::size_t j>
    decltype(auto) evalCouplingResidual(Dune::index_constant<faceIdx> domainI,
                                        const SubControlVolumeFace<cellCenterIdx>& scvfI, // TODO
                                        const LocalAssemblerI& localAssemblerI,
                                        Dune::index_constant<j> domainJ,
                                        std::size_t dofIdxGlobalJ) const
    {
        static_assert(domainI != domainJ, "Domain i cannot be coupled to itself!");
        return localAssemblerI.evalLocalResidualForFace(scvfI);
    }

    /*!
     * \brief return the numeric epsilon used for deflecting primary variables of coupled domain i.
     * \note  specialization for non-staggered schemes
     */
    template<std::size_t i, typename std::enable_if_t<(GridGeometry<i>::discMethod != DiscretizationMethod::staggered), int> = 0>
    decltype(auto) numericEpsilon(Dune::index_constant<i> id,
                                  const std::string& paramGroup) const
    {
        return ParentType::numericEpsilon(id, paramGroup);
    }

    /*!
     * \brief return the numeric epsilon used for deflecting primary variables of coupled domain i.
     * \note  specialization for non-staggered schemes
     */
    template<std::size_t i, typename std::enable_if_t<(GridGeometry<i>::discMethod == DiscretizationMethod::staggered), int> = 0>
    decltype(auto) numericEpsilon(Dune::index_constant<i>,
                                  const std::string& paramGroup) const
    {
        constexpr std::size_t numEqCellCenter = Traits::template SubDomain<cellCenterIdx>::PrimaryVariables::dimension;
        constexpr std::size_t numEqFace = Traits::template SubDomain<faceIdx>::PrimaryVariables::dimension;
        constexpr bool isCellCenter = GridGeometry<i>::isCellCenter();
        constexpr std::size_t numEq = isCellCenter ? numEqCellCenter : numEqFace;
        constexpr auto prefix = isCellCenter ? "CellCenter" : "Face";

        try {
            if(paramGroup.empty())
                return NumericEpsilon<typename Traits::Scalar, numEq>(prefix);
            else
                return NumericEpsilon<typename Traits::Scalar, numEq>(paramGroup + "." + prefix);
        }
        catch (Dune::RangeError& e)
        {
            DUNE_THROW(Dumux::ParameterException, "For the staggered model, you can to specify \n\n"
                        "  CellCenter.Assembly.NumericDifference.PriVarMagnitude = mCC\n"
                        "  Face.Assembly.NumericDifference.PriVarMagnitude = mFace\n"
                        "  CellCenter.Assembly.NumericDifference.BaseEpsilon = eCC_0 ... eCC_numEqCellCenter-1\n"
                        "  Face.Assembly.NumericDifference.BaseEpsilon = eFace_0 ... eFace_numEqFace-1\n\n"
                        "Wrong numer of values set for " << prefix  << " (has " << numEq << " primary variable(s))\n\n" << e);
        }

    }

    /*!
     * \brief Access the coupling context needed for the Darcy domain
     */
    const auto& faceCouplingContext(const Element<faceIdx>& element, const SubControlVolume<faceIdx>& scv) const
    {
        const auto eIdx = this->problem(faceIdx).gridGeometry().elementMapper().index(element);
        if (!faceCouplingContext_ || boundFaceElemIdx_ != eIdx)
        {
            boundFaceElemIdx_ = eIdx;
            const auto& curSol = this->curSol();
            auto cellCenterFvGeometry = localView(this->problem(cellCenterIdx).gridGeometry());
            auto cellCenterElemVolVars = localView(gridVariables(cellCenterIdx).curGridVolVars());
            cellCenterFvGeometry.bind(element);
            cellCenterElemVolVars.bind(element, cellCenterFvGeometry, curSol);
            faceCouplingContext_.emplace(FaceCouplingContext{cellCenterFvGeometry, cellCenterElemVolVars});
        }
        return *faceCouplingContext_;
    }

protected:
    //! the grid variables of domain i
    template<std::size_t i>
    GridVariables<i>& gridVariables(Dune::index_constant<i> domainId) const
    { return *std::get<domainId>(gridVariables_); }

    template<std::size_t i>
    auto& getVolVarAccess(Dune::index_constant<i> domainId, const SubControlVolume<i>& scv)
    {
        GridVolumeVariables<i>& gridVolVars = gridVariables(domainId).curGridVolVars();
        ElementVolumeVariables<i>& elemVolVars = faceCouplingContext_->elemVolVars;

        if constexpr (GridVolumeVariables<i>::cachingEnabled)
            return gridVolVars.volVars(scv);
        else
            return elemVolVars[scv];
    }

private:


    /*!
     * \brief A tuple of std::shared_ptrs to the grid variables of the sub problems
     */
     mutable GridVariablesTuple gridVariables_;

     const SolutionVector* prevSol_ = nullptr;

     mutable std::optional<FaceCouplingContext> faceCouplingContext_;

     mutable std::size_t boundFaceElemIdx_;

};

} //end namespace Dumux

#endif
