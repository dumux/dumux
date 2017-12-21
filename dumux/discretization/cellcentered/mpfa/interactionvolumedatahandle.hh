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
 * \brief Data handle class for interaction volumes of mpfa methods.
 *        This class is passed to interaction volumes to store the necessary data in it.
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_INTERACTIONVOLUMEDATAHANDLE_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_INTERACTIONVOLUMEDATAHANDLE_HH

#include <dune/common/exceptions.hh>
#include <dumux/common/properties.hh>

namespace Dumux
{
    //! Empty data handle class
    template<class TypeTag>
    class EmptyDataHandle
    {
        //! we use the dynamic types here to be compatible on the boundary
        using InteractionVolume = typename GET_PROP_TYPE(TypeTag, PrimaryInteractionVolume);
        using DirichletDataContainer = typename InteractionVolume::DirichletDataContainer;
        using GlobalIndexContainer = typename InteractionVolume::Traits::DynamicGlobalIndexContainer;
        using Matrix = typename InteractionVolume::Traits::DynamicMatrix;
        using Vector = typename InteractionVolume::Traits::DynamicVector;

    public:
        //! diffusion caches need to set phase and component index
        void setDiffusionContext(unsigned int phaseIdx, unsigned int compIdx) {}

        //! functions to set the size of the matrices
        void resizeT(unsigned int n, unsigned int m) {}
        void resizeAB(unsigned int n, unsigned int m) {}
        void resizeOutsideTij(unsigned int n, unsigned int m) {}

        //! functions to set the pointers to the stencil
        void setVolVarsStencilPointer(const GlobalIndexContainer& stencil) {}

        //! return functions for the stored data
        const GlobalIndexContainer& volVarsStencil() const { return throw_<const GlobalIndexContainer&>(); }
        const DirichletDataContainer& dirichletData() const { return throw_<DirichletDataContainer&>(); }

        const Matrix& T() const { return throw_<const Matrix&>(); }
        Matrix& T() { return throw_<Matrix&>(); }

        const Matrix& AB() const { return throw_<const Matrix&>(); }
        Matrix& AB() { return throw_<Matrix&>(); }

        const Matrix& outsideTij() const { return throw_<const Matrix&>(); }
        Matrix& outsideTij() { return throw_<Matrix&>(); }

    private:
        template<class ReturnType>
        ReturnType throw_() const { DUNE_THROW(Dune::InvalidStateException, "Trying to access data for a deactivated physical process"); }
    };

    //! Data handle for quantities related to advection
    template<class TypeTag, bool EnableAdvection>
    class AdvectionDataHandle
    {
        //! we use the dynamic types here to be compatible on the boundary
        using InteractionVolume = typename GET_PROP_TYPE(TypeTag, PrimaryInteractionVolume);
        using DirichletDataContainer = typename InteractionVolume::DirichletDataContainer;
        using GlobalIndexContainer = typename InteractionVolume::Traits::DynamicGlobalIndexContainer;
        using Matrix = typename InteractionVolume::Traits::DynamicMatrix;
        using Vector = typename InteractionVolume::Traits::DynamicVector;

    public:
        //! set the sizes of the matrices
        void resizeT(unsigned int n, unsigned int m) { advectionT_.resize(n, m); }
        void resizeAB(unsigned int n, unsigned int m) { advectionAB_.resize(n, m); }
        void resizeOutsideTij(unsigned int n, unsigned int m) { advectionTout_.resize(n, m); }

        //! sets the pointer to the stencil
        void setVolVarsStencilPointer(const GlobalIndexContainer& stencil) { advectionVolVarsStencil_ = &stencil; }

        //! return functions for the stored data
        const GlobalIndexContainer& volVarsStencil() const { return *advectionVolVarsStencil_; }

        const Matrix& T() const { return advectionT_; }
        Matrix& T() { return advectionT_; }

        const Matrix& AB() const { return advectionAB_; }
        Matrix& AB() { return advectionAB_; }

        const Matrix& outsideTij() const { return advectionTout_; }
        Matrix& outsideTij() { return advectionTout_; }

    private:
        // advection-related variables
        const GlobalIndexContainer* advectionVolVarsStencil_;  //!< Pointer to the global volvar indices (stored in the interaction volume)
        Matrix advectionT_;                                    //!< The transmissibilities
        Matrix advectionAB_;                                   //!< Coefficients for gradient reconstruction
        Matrix advectionTout_;                                 //!< The transmissibilities associated with "outside" faces (only necessary on surface grids)
    };

    //! Data handle for quantities related to diffusion
    template<class TypeTag, bool EnableDiffusion>
    class DiffusionDataHandle
    {
        //! we use the dynamic types here to be compatible on the boundary
        using InteractionVolume = typename GET_PROP_TYPE(TypeTag, PrimaryInteractionVolume);
        using DirichletDataContainer = typename InteractionVolume::DirichletDataContainer;
        using GlobalIndexContainer = typename InteractionVolume::Traits::DynamicGlobalIndexContainer;
        using Matrix = typename InteractionVolume::Traits::DynamicMatrix;
        using Vector = typename InteractionVolume::Traits::DynamicVector;

        static constexpr int numPhases = GET_PROP_VALUE(TypeTag, NumPhases);
        static constexpr int numComponents = GET_PROP_VALUE(TypeTag, NumComponents);

    public:
        //! diffusion caches need to set phase and component index
        void setDiffusionContext(unsigned int phaseIdx, unsigned int compIdx)
        {
            contextPhaseIdx_ = phaseIdx;
            contextCompIdx_ = compIdx;
        }

        //! set the sizes of the matrices
        void resizeT(unsigned int n, unsigned int m)
        {
            for (auto& array : diffusionT_)
                for (auto& matrix : array)
                    matrix.resize(n, m);
        }

        void resizeAB(unsigned int n, unsigned int m)
        {
            for (auto& array : diffusionAB_)
                for (auto& matrix : array)
                    matrix.resize(n, m);
        }

        void resizeOutsideTij(unsigned int n, unsigned int m)
        {
            for (auto& array : diffusionTout_)
                for (auto& matrix : array)
                    matrix.resize(n, m);
        }

        //! sets the pointer to stencil
        void setVolVarsStencilPointer(const GlobalIndexContainer& stencil)
        {
            diffusionVolVarsStencil_[contextPhaseIdx_][contextCompIdx_] = &stencil;
        }

        //! return functions for the stored data
        const GlobalIndexContainer& volVarsStencil() const
        { return *diffusionVolVarsStencil_[contextPhaseIdx_][contextCompIdx_]; }

        const Matrix& T() const { return diffusionT_[contextPhaseIdx_][contextCompIdx_]; }
        Matrix& T() { return diffusionT_[contextPhaseIdx_][contextCompIdx_]; }

        const Matrix& AB() const { return diffusionAB_[contextPhaseIdx_][contextCompIdx_]; }
        Matrix& AB() { return diffusionAB_[contextPhaseIdx_][contextCompIdx_]; }

        const Matrix& outsideTij() const { return diffusionTout_[contextPhaseIdx_][contextCompIdx_]; }
        Matrix& outsideTij() { return diffusionTout_[contextPhaseIdx_][contextCompIdx_]; }

    private:
        // diffusion-related variables (see comments in AdvectionDataHandle)
        unsigned int contextPhaseIdx_;                         //!< The phase index set for the context
        unsigned int contextCompIdx_;                          //!< The component index set for the context
        std::array<std::array<const GlobalIndexContainer*, numComponents>, numPhases> diffusionVolVarsStencil_;
        std::array<std::array<Matrix, numComponents>, numPhases> diffusionT_;
        std::array<std::array<Matrix, numComponents>, numPhases> diffusionAB_;
        std::array<std::array<Matrix, numComponents>, numPhases> diffusionTout_;
    };

    //! Data handle for quantities related to advection
    template<class TypeTag, bool EnableHeatConduction>
    class HeatConductionDataHandle
    {
        //! we use the dynamic types here to be compatible on the boundary
        using InteractionVolume = typename GET_PROP_TYPE(TypeTag, PrimaryInteractionVolume);
        using DirichletDataContainer = typename InteractionVolume::DirichletDataContainer;
        using GlobalIndexContainer = typename InteractionVolume::Traits::DynamicGlobalIndexContainer;
        using Matrix = typename InteractionVolume::Traits::DynamicMatrix;
        using Vector = typename InteractionVolume::Traits::DynamicVector;

    public:
        //! set the sizes of the matrices
        void resizeT(unsigned int n, unsigned int m) { heatConductionT_.resize(n, m); }
        void resizeAB(unsigned int n, unsigned int m) { heatConductionAB_.resize(n, m); }
        void resizeOutsideTij(unsigned int n, unsigned int m) { heatConductionTout_.resize(n, m); }

        //! sets the pointer to the stencil
        void setVolVarsStencilPointer(const GlobalIndexContainer& stencil) { heatConductionVolVarsStencil_ = &stencil; }

        //! return functions for the stored data
        const GlobalIndexContainer& volVarsStencil() const { return *heatConductionVolVarsStencil_; }

        const Matrix& T() const { return heatConductionT_; }
        Matrix& T() { return heatConductionT_; }

        const Matrix& AB() const { return heatConductionAB_; }
        Matrix& AB() { return heatConductionAB_; }

        const Matrix& outsideTij() const { return heatConductionTout_; }
        Matrix& outsideTij() { return heatConductionTout_; }

    private:
        // heat conduction-related variables
        const GlobalIndexContainer* heatConductionVolVarsStencil_;  //!< Pointer to the global volvar indices (stored in the interaction volume)
        Matrix heatConductionT_;                                    //!< The transmissibilities
        Matrix heatConductionAB_;                                   //!< Coefficients for gradient reconstruction
        Matrix heatConductionTout_;                                 //!< The transmissibilities associated with "outside" faces (only necessary on surface grids)
    };

    //! Process-dependet data handle when related process is disabled
    template<class TypeTag> class AdvectionDataHandle<TypeTag, false> : public EmptyDataHandle<TypeTag> {};
    template<class TypeTag> class DiffusionDataHandle<TypeTag, false> : public EmptyDataHandle<TypeTag> {};
    template<class TypeTag> class HeatConductionDataHandle<TypeTag, false> : public EmptyDataHandle<TypeTag> {};

    //! Interaction volume data handle class
    template<class TypeTag>
    class InteractionVolumeDataHandle : public AdvectionDataHandle<TypeTag, GET_PROP_VALUE(TypeTag, EnableAdvection)>,
                                        public DiffusionDataHandle<TypeTag, GET_PROP_VALUE(TypeTag, EnableMolecularDiffusion)>,
                                        public HeatConductionDataHandle<TypeTag, GET_PROP_VALUE(TypeTag, EnableEnergyBalance)>
    {
        using AdvectionHandle = AdvectionDataHandle<TypeTag, GET_PROP_VALUE(TypeTag, EnableAdvection)>;
        using DiffusionHandle = DiffusionDataHandle<TypeTag, GET_PROP_VALUE(TypeTag, EnableMolecularDiffusion)>;
        using HeatConductionHandle = HeatConductionDataHandle<TypeTag, GET_PROP_VALUE(TypeTag, EnableEnergyBalance)>;

        //! we use the dynamic types here to be compatible on the boundary
        using InteractionVolume = typename GET_PROP_TYPE(TypeTag, PrimaryInteractionVolume);
        using DirichletDataContainer = typename InteractionVolume::DirichletDataContainer;
        using GlobalIndexContainer = typename InteractionVolume::Traits::DynamicGlobalIndexContainer;
        using Matrix = typename InteractionVolume::Traits::DynamicMatrix;
        using Vector = typename InteractionVolume::Traits::DynamicVector;

    public:
        enum class Contexts : unsigned int
        {
            undefined,
            advection,
            diffusion,
            heatConduction
        };

        //! The constructor
        InteractionVolumeDataHandle() : context_(Contexts::undefined) {}

        //! set the context of the cache
        void setAdvectionContext() { context_ = Contexts::advection; }
        void setHeatConductionContext() { context_ = Contexts::heatConduction; }
        void setDiffusionContext(unsigned int phaseIdx, unsigned int compIdx)
        {
            context_ = Contexts::diffusion;
            DiffusionHandle::setDiffusionContext(phaseIdx, compIdx);
        }

        //! returns the current context
        Contexts getContext() const { return context_; }

        //! set the sizes of the matrices
        void resizeT(unsigned int n, unsigned int m)
        {
            AdvectionHandle::resizeT(n, m);
            DiffusionHandle::resizeT(n, m);
            HeatConductionHandle::resizeT(n, m);
        }

        void resizeAB(unsigned int n, unsigned int m)
        {
            AdvectionHandle::resizeAB(n, m);
            DiffusionHandle::resizeAB(n, m);
            HeatConductionHandle::resizeAB(n, m);
        }

        void resizeOutsideTij(unsigned int n, unsigned int m)
        {
            AdvectionHandle::resizeOutsideTij(n, m);
            DiffusionHandle::resizeOutsideTij(n, m);
            HeatConductionHandle::resizeOutsideTij(n, m);
        }

        //! sets the pointer to the stencil
        void setVolVarsStencilPointer(const GlobalIndexContainer& stencil)
        {
            if (context_ == Contexts::advection)
                AdvectionHandle::setVolVarsStencilPointer(stencil);
            else if (context_ == Contexts::diffusion)
                DiffusionHandle::setVolVarsStencilPointer(stencil);
            else if (context_ == Contexts::heatConduction)
                HeatConductionHandle::setVolVarsStencilPointer(stencil);
            else
                DUNE_THROW(Dune::InvalidStateException, "No valid context set!");
        }

        //! sets the dirichlet data container
        void setDirichletData(DirichletDataContainer&& data)
        {
            dirichletData_ = std::move(data);
        }

        //! return functions for the stored data
        const DirichletDataContainer& dirichletData() const { return dirichletData_; }

        const GlobalIndexContainer& volVarsStencil() const
        {
            if (context_ == Contexts::advection)
                return AdvectionHandle::volVarsStencil();
            else if (context_ == Contexts::diffusion)
                return DiffusionHandle::volVarsStencil();
            else if (context_ == Contexts::heatConduction)
                return HeatConductionHandle::volVarsStencil();
            else
                DUNE_THROW(Dune::InvalidStateException, "No valid context set!");
        }

        const Matrix& T() const
        {
            if (context_ == Contexts::advection)
                return AdvectionHandle::T();
            else if (context_ == Contexts::diffusion)
                return DiffusionHandle::T();
            else if (context_ == Contexts::heatConduction)
                return HeatConductionHandle::T();
            else
                DUNE_THROW(Dune::InvalidStateException, "No valid context set!");
        }

        Matrix& T()
        {
            if (context_ == Contexts::advection)
                return AdvectionHandle::T();
            else if (context_ == Contexts::diffusion)
                return DiffusionHandle::T();
            else if (context_ == Contexts::heatConduction)
                return HeatConductionHandle::T();
            else
                DUNE_THROW(Dune::InvalidStateException, "No valid context set!");
        }

        const Matrix& AB() const
        {
            if (context_ == Contexts::advection)
                return AdvectionHandle::AB();
            else if (context_ == Contexts::diffusion)
                return DiffusionHandle::AB();
            else if (context_ == Contexts::heatConduction)
                return HeatConductionHandle::AB();
            else
                DUNE_THROW(Dune::InvalidStateException, "No valid context set!");
        }

        Matrix& AB()
        {
            if (context_ == Contexts::advection)
                return AdvectionHandle::AB();
            else if (context_ == Contexts::diffusion)
                return DiffusionHandle::AB();
            else if (context_ == Contexts::heatConduction)
                return HeatConductionHandle::AB();
            else
                DUNE_THROW(Dune::InvalidStateException, "No valid context set!");
        }

        const Matrix& outsideTij() const
        {
            if (context_ == Contexts::advection)
                return AdvectionHandle::outsideTij();
            else if (context_ == Contexts::diffusion)
                return DiffusionHandle::outsideTij();
            else if (context_ == Contexts::heatConduction)
                return HeatConductionHandle::outsideTij();
            else
                DUNE_THROW(Dune::InvalidStateException, "No valid context set!");
        }

        Matrix& outsideTij()
        {
            if (context_ == Contexts::advection)
                return AdvectionHandle::outsideTij();
            else if (context_ == Contexts::diffusion)
                return DiffusionHandle::outsideTij();
            else if (context_ == Contexts::heatConduction)
                return HeatConductionHandle::outsideTij();
            else
                DUNE_THROW(Dune::InvalidStateException, "No valid context set!");
        }

    private:
        Contexts context_;                     //!< The context variable
        DirichletDataContainer dirichletData_; //!< The dirichlet data container of this iv
    };

} // end namespace Dumux

#endif
