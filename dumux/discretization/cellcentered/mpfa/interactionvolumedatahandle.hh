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

namespace Dumux
{
    //! forward declaration of the implementation class
    //! Struct to store the transmissibilities and other data separate from the interaction volume
    template<class TypeTag, bool Advection, bool Diffusion, bool Energy>
    struct InteractionVolumeDataHandleImplementation;

    template<class TypeTag>
    using InteractionVolumeDataHandle = InteractionVolumeDataHandleImplementation<TypeTag,
                                                                                  GET_PROP_VALUE(TypeTag, EnableAdvection),
                                                                                  GET_PROP_VALUE(TypeTag, EnableMolecularDiffusion),
                                                                                  GET_PROP_VALUE(TypeTag, EnableEnergyBalance)>;

    //! Specialization for purely advective problems
    template<class TypeTag>
    struct InteractionVolumeDataHandleImplementation<TypeTag, true, false, false>
    {
    private:
        using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
        using GridView = typename GET_PROP_TYPE(TypeTag, GridView);

        //! we use the dynamic types here to be compatible on the boundary
        using InteractionVolume = typename GET_PROP_TYPE(TypeTag, PrimaryInteractionVolume);
        using DirichletDataContainer = typename InteractionVolume::DirichletDataContainer;
        using GlobalIndexContainer = typename InteractionVolume::Traits::DynamicGlobalIndexContainer;
        using Matrix = typename InteractionVolume::Traits::DynamicMatrix;
        using Vector = typename InteractionVolume::Traits::DynamicVector;

    public:
        enum class Contexts : unsigned int
        {
            advection
        };

        //! The constructor
        InteractionVolumeDataHandleImplementation() {}

        //! The context is always advection in this specialization
        void setAdvectionContext() {}

        //! returns the current context
        Contexts getContext() const { return Contexts::advection; }

        //! functions to set the size of the matrices
        void resizeT(unsigned int n, unsigned int m) { T_.resize(n, m); }
        void resizeAB(unsigned int n, unsigned int m) { AB_.resize(n, m); }
        void resizeOutsideTij(unsigned int n, unsigned int m) { tijOut_.resize(n, m); }

        //! functions to set the pointers to stencil and Dirichlet data
        void setVolVarsStencilPointer(const GlobalIndexContainer& stencil) { volVarsStencil_ = &stencil; }
        void setDirichletDataPointer(const DirichletDataContainer& data) { dirichletData_ = &data; }

        //! return functions for the stored data
        const GlobalIndexContainer& volVarsStencil() const { return *volVarsStencil_; }
        const DirichletDataContainer& dirichletData() const { return *dirichletData_; }

        const Matrix& T() const { return T_; }
        Matrix& T() { return T_; }

        const Matrix& AB() const { return AB_; }
        Matrix& AB() { return AB_; }

        const Matrix& outsideTij() const { return tijOut_; }
        Matrix& outsideTij() { return tijOut_; }

    private:
        const GlobalIndexContainer* volVarsStencil_;  //! Pointer to the global volvar indices (stored in the interaction volume)
        const DirichletDataContainer* dirichletData_; //! Container with dirichlet data of the iv

        Matrix T_;                                    //! The transmissibilities
        Matrix AB_;                                   //! Coefficients for gradient reconstruction
        Matrix tijOut_;                               //! The transmissibilities associated with "outside" faces (only necessary on surface grids)
    };

    //! Specialization for advective-diffusive problems
    template<class TypeTag>
    struct InteractionVolumeDataHandleImplementation<TypeTag, true, true, false>
    {
    private:
        using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
        using GridView = typename GET_PROP_TYPE(TypeTag, GridView);

        //! we use the dynamic types here to be compatible on the boundary
        using InteractionVolume = typename GET_PROP_TYPE(TypeTag, PrimaryInteractionVolume);
        using DirichletDataContainer = typename InteractionVolume::DirichletDataContainer;
        using GlobalIndexContainer = typename InteractionVolume::Traits::DynamicGlobalIndexContainer;
        using Matrix = typename InteractionVolume::Traits::DynamicMatrix;
        using Vector = typename InteractionVolume::Traits::DynamicVector;

        static constexpr int numPhases = GET_PROP_VALUE(TypeTag, NumPhases);
        static constexpr int numComponents = GET_PROP_VALUE(TypeTag, NumComponents);

    public:
        enum class Contexts : unsigned int
        {
            undefined,
            advection,
            diffusion
        };

        //! The constructor
        InteractionVolumeDataHandleImplementation() : context_(Contexts::undefined) {}

        //! sets the context of the cache to advection
        void setAdvectionContext()
        { context_ = Contexts::advection; }

        //! sets the context of the cache to diffusion
        void setDiffusionContext(unsigned int phaseIdx, unsigned int compIdx)
        {
            context_ = Contexts::diffusion;
            contextPhaseIdx_ = phaseIdx;
            contextCompIdx_ = compIdx;
        }

        //! returns the current context
        Contexts getContext() const
        { return context_; }

        //! functions to set the size of the matrices
        void resizeT(unsigned int n, unsigned int m)
        {
            advectionT_.resize(n, m);
            for (auto& array : diffusionT_)
                for (auto& matrix : array)
                    matrix.resize(n, m);

        }

        void resizeAB(unsigned int n, unsigned int m)
        {
            advectionAB_.resize(n, m);
            for (auto& array : diffusionAB_)
                for (auto& matrix : array)
                    matrix.resize(n, m);
        }

        void resizeOutsideTij(unsigned int n, unsigned int m)
        {
            advectionTout_.resize(n, m);
            for (auto& array : diffusionTout_)
                for (auto& matrix : array)
                    matrix.resize(n, m);
        }

        //! functions to set the pointers to stencil and Dirichlet data
        void setVolVarsStencilPointer(const GlobalIndexContainer& stencil)
        {
            assert(context_ != Contexts::undefined && "No valid context set!");
            if (context_ == Contexts::advection) advectionVolVarsStencil_ = &stencil;
            else diffusionVolVarsStencil_[contextPhaseIdx_][contextCompIdx_] = &stencil;
        }

        void setDirichletDataPointer(const DirichletDataContainer& data)
        {
            assert(context_ != Contexts::undefined && "No valid context set!");
            if (context_ == Contexts::advection) advectionDirichletData_ = &data;
            else diffusionDirichletData_[contextPhaseIdx_][contextCompIdx_] = &data;
        }

        //! return functions for the stored data
        const GlobalIndexContainer& volVarsStencil() const
        {
            assert(context_ != Contexts::undefined && "No valid context set!");
            return context_ == Contexts::advection ? *advectionVolVarsStencil_ : *(diffusionVolVarsStencil_[contextPhaseIdx_][contextCompIdx_]);
        }

        const DirichletDataContainer& dirichletData() const
        {
            assert(context_ != Contexts::undefined && "No valid context set!");
            return context_ == Contexts::advection ? *advectionDirichletData_ : *(diffusionDirichletData_[contextPhaseIdx_][contextCompIdx_]);
        }

        const Matrix& T() const
        {
            assert(context_ != Contexts::undefined && "No valid context set!");
            return context_ == Contexts::advection ? advectionT_ : diffusionT_[contextPhaseIdx_][contextCompIdx_];
        }

        Matrix& T()
        {
            assert(context_ != Contexts::undefined && "No valid context set!");
            return context_ == Contexts::advection ? advectionT_ : diffusionT_[contextPhaseIdx_][contextCompIdx_];
        }

        const Matrix& AB() const
        {
            assert(context_ != Contexts::undefined && "No valid context set!");
            return context_ == Contexts::advection ? advectionAB_ : diffusionAB_[contextPhaseIdx_][contextCompIdx_];
        }

        Matrix& AB()
        {
            assert(context_ != Contexts::undefined && "No valid context set!");
            return context_ == Contexts::advection ? advectionAB_ : diffusionAB_[contextPhaseIdx_][contextCompIdx_];
        }

        const Matrix& outsideTij() const
        {
            assert(context_ != Contexts::undefined && "No valid context set!");
            return context_ == Contexts::advection ? advectionTout_ : diffusionTout_[contextPhaseIdx_][contextCompIdx_];
        }

        Matrix& outsideTij()
        {
            assert(context_ != Contexts::undefined && "No valid context set!");
            return context_ == Contexts::advection ? advectionTout_ : diffusionTout_[contextPhaseIdx_][contextCompIdx_];
        }

    private:
        Contexts context_;                                     //! The context variable

        // advection-related variables
        const GlobalIndexContainer* advectionVolVarsStencil_;  //! Pointer to the global volvar indices (stored in the interaction volume)
        const DirichletDataContainer* advectionDirichletData_; //! Pointer to the container with dirichlet data of the iv
        Matrix advectionT_;                                    //! The transmissibilities
        Matrix advectionAB_;                                   //! Coefficients for gradient reconstruction
        Matrix advectionTout_;                                 //! The transmissibilities associated with "outside" faces (only necessary on surface grids)

        // diffusion-related variables (see comments above)
        unsigned int contextPhaseIdx_;                         //! The phase index set for the context
        unsigned int contextCompIdx_;                          //! The component index set for the context
        std::array<std::array<const GlobalIndexContainer*, numComponents>, numPhases> diffusionVolVarsStencil_;
        std::array<std::array<const DirichletDataContainer*, numComponents>, numPhases> diffusionDirichletData_;
        std::array<std::array<Matrix, numComponents>, numPhases> diffusionT_;
        std::array<std::array<Matrix, numComponents>, numPhases> diffusionAB_;
        std::array<std::array<Matrix, numComponents>, numPhases> diffusionTout_;
    };

    // //! Specialization for problems involving advection, diffusion and heat conduction
    // template<class TypeTag>
    // struct InteractionVolumeDataHandleImplementation<TypeTag, true, true, true>
    // {
    // private:
    //     using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    //     using GridView = typename GET_PROP_TYPE(TypeTag, GridView);

    //     //! we use the types from the boundary interaction volume traits to be compatible on the boundary
    //     using BoundaryInteractionVolume = typename GET_PROP_TYPE(TypeTag, BoundaryInteractionVolume);
    //     using GlobalIndexSet = typename BoundaryInteractionVolume::Traits::GlobalIndexSet;
    //     using PositionVector = typename BoundaryInteractionVolume::Traits::PositionVector;
    //     using Matrix = typename BoundaryInteractionVolume::Traits::Matrix;
    //     using Vector = typename BoundaryInteractionVolume::Traits::Vector;

    //     static constexpr int numPhases = GET_PROP_VALUE(TypeTag, NumPhases);
    //     static constexpr int numComponents = GET_PROP_VALUE(TypeTag, NumComponents);

    // public:
    //     enum class Contexts : unsigned int
    //     {
    //         undefined,
    //         advection,
    //         diffusion,
    //         heatConduction
    //     };

    //     //! The constructor
    //     InteractionVolumeDataHandleImplementation() : context_(Contexts::undefined) {}

    //     //! sets the context of the cache to advection
    //     void setAdvectionContext()
    //     {
    //         context_ = Contexts::advection;
    //     }

    //     //! sets the context of the cache to diffusion
    //     void setDiffusionContext(unsigned int phaseIdx, unsigned int compIdx)
    //     {
    //         context_ = Contexts::diffusion;
    //         contextPhaseIdx_ = phaseIdx;
    //         contextCompIdx_ = compIdx;
    //     }

    //     //! sets the context of the cache to advection
    //     void setHeatConductionContext()
    //     {
    //         context_ = Contexts::heatConduction;
    //     }

    //     //! returns the current context
    //     Contexts getContext() const
    //     { return context_; }

    //     //! functions to set the size of the matrices
    //     void resizeT(unsigned int n, unsigned int m)
    //     {
    //         advectionT_.resize(n, m);
    //         heatConductionT_.resize(n, m);
    //         for (auto& array : diffusionT_)
    //             for (auto& matrix : array)
    //                 matrix.resize(n, m);

    //     }

    //     void resizeCA(unsigned int n, unsigned int m)
    //     {
    //         advectionCA_.resize(n, m);
    //         heatConductionCA_.resize(n, m);
    //         for (auto& array : diffusionCA_)
    //             for (auto& matrix : array)
    //                 matrix.resize(n, m);
    //     }

    //     void resizeAB(unsigned int n, unsigned int m)
    //     {
    //         advectionAB_.resize(n, m);
    //         heatConductionAB_.resize(n, m);
    //         for (auto& array : diffusionAB_)
    //             for (auto& matrix : array)
    //                 matrix.resize(n, m);
    //     }

    //     void resizeOutsideTij(unsigned int n, unsigned int m)
    //     {
    //         advectionTout_.resize(n, m);
    //         heatConductionTout_.resize(n, m);
    //         for (auto& array : diffusionTout_)
    //             for (auto& matrix : array)
    //                 matrix.resize(n, m);
    //     }

    //     //! functions to set the pointers to stencil and positions
    //     void setVolVarsStencilPointer(const GlobalIndexSet& stencil)
    //     {
    //         assert(context_ != Contexts::undefined && "No valid context set!");
    //         if (context_ == Contexts::diffusion)
    //             diffusionVolVarsStencil_[contextPhaseIdx_][contextCompIdx_] = &stencil;
    //         else if (context_ == Contexts::advection)
    //             advectionVolVarsStencil_ = &stencil;
    //         else
    //             heatConductionVolVarsStencil_ = &stencil;
    //     }

    //     void setVolVarsPositionsPointer(const PositionVector& pos)
    //     {
    //         assert(context_ != Contexts::undefined && "No valid context set!");
    //         if (context_ == Contexts::diffusion)
    //             diffusionVolVarsPositions_[contextPhaseIdx_][contextCompIdx_] = &pos;
    //         else if (context_ == Contexts::advection)
    //             advectionVolVarsPositions_ = &pos;
    //         else
    //             heatConductionVolVarsPositions_ = &pos;
    //     }

    //     //! return functions for the stored data
    //     const GlobalIndexSet& volVarsStencil() const
    //     {
    //         assert(context_ != Contexts::undefined && "No valid context set!");
    //         if (context_ == Contexts::diffusion)
    //             return *(diffusionVolVarsStencil_[contextPhaseIdx_][contextCompIdx_]);
    //         else if (context_ == Contexts::advection)
    //             return *advectionVolVarsStencil_;
    //         else
    //             return *heatConductionVolVarsStencil_;
    //     }

    //     const PositionVector& volVarsPositions() const
    //     {
    //         assert(context_ != Contexts::undefined && "No valid context set!");
    //         if (context_ == Contexts::diffusion)
    //             return *(diffusionVolVarsPositions_[contextPhaseIdx_][contextCompIdx_]);
    //         else if (context_ == Contexts::advection)
    //             return *advectionVolVarsPositions_;
    //         else
    //             return *heatConductionVolVarsPositions_;
    //     }

    //     //! return functions for the stored matrices
    //     const Matrix& T() const
    //     {
    //         assert(context_ != Contexts::undefined && "No valid context set!");
    //         if (context_ == Contexts::diffusion)
    //             return diffusionT_[contextPhaseIdx_][contextCompIdx_];
    //         else if (context_ == Contexts::advection)
    //             return advectionT_;
    //         else
    //             return heatConductionT_;
    //     }

    //     Matrix& T()
    //     {
    //         assert(context_ != Contexts::undefined && "No valid context set!");
    //         if (context_ == Contexts::diffusion)
    //             return diffusionT_[contextPhaseIdx_][contextCompIdx_];
    //         else if (context_ == Contexts::advection)
    //             return advectionT_;
    //         else
    //             return heatConductionT_;
    //     }

    //     const Matrix& CA() const
    //     {
    //         assert(context_ != Contexts::undefined && "No valid context set!");
    //         if (context_ == Contexts::diffusion)
    //             return diffusionCA_[contextPhaseIdx_][contextCompIdx_];
    //         else if (context_ == Contexts::advection)
    //             return advectionCA_;
    //         else
    //             return heatConductionCA_;
    //     }

    //     Matrix& CA()
    //     {
    //         assert(context_ != Contexts::undefined && "No valid context set!");
    //         if (context_ == Contexts::diffusion)
    //             return diffusionCA_[contextPhaseIdx_][contextCompIdx_];
    //         else if (context_ == Contexts::advection)
    //             return advectionCA_;
    //         else
    //             return heatConductionCA_;
    //     }

    //     const Matrix& AB() const
    //     {
    //         assert(context_ != Contexts::undefined && "No valid context set!");
    //         if (context_ == Contexts::diffusion)
    //             return diffusionAB_[contextPhaseIdx_][contextCompIdx_];
    //         else if (context_ == Contexts::advection)
    //             return advectionAB_;
    //         else
    //             return heatConductionAB_;
    //     }

    //     Matrix& AB()
    //     {
    //         assert(context_ != Contexts::undefined && "No valid context set!");
    //         if (context_ == Contexts::diffusion)
    //             return diffusionAB_[contextPhaseIdx_][contextCompIdx_];
    //         else if (context_ == Contexts::advection)
    //             return advectionAB_;
    //         else
    //             return heatConductionAB_;
    //     }

    //     const Matrix& outsideTij() const
    //     {
    //         assert(context_ != Contexts::undefined && "No valid context set!");
    //         if (context_ == Contexts::diffusion)
    //             return diffusionTout_[contextPhaseIdx_][contextCompIdx_];
    //         else if (context_ == Contexts::advection)
    //             return advectionTout_;
    //         else
    //             return heatConductionTout_;
    //     }

    //     Matrix& outsideTij()
    //     {
    //         assert(context_ != Contexts::undefined && "No valid context set!");
    //         if (context_ == Contexts::diffusion)
    //             return diffusionTout_[contextPhaseIdx_][contextCompIdx_];
    //         else if (context_ == Contexts::advection)
    //             return advectionTout_;
    //         else
    //             return heatConductionTout_;
    //     }

    // private:
    //     Contexts context_;                                   //! The context variable

    //     // advection-related variables
    //     const GlobalIndexSet* advectionVolVarsStencil_;      //! Pointer to the global volvar indices (stored in the interaction volume)
    //     const PositionVector* advectionVolVarsPositions_;    //! Pointer to the positions of the vol vars (stored in the interaction volume)
    //     Matrix advectionT_;                                  //! The transmissibilities
    //     Matrix advectionCA_;                                 //! The coefficients for the neumann flux transformations
    //     Matrix advectionAB_;                                 //! Coefficients for gradient reconstruction
    //     Matrix advectionTout_;                  //! The transmissibilities associated with "outside" faces (only necessary on surface grids)

    //     // diffusion-related variables (see comments above)
    //     unsigned int contextPhaseIdx_;                       //! The phase index set for the context
    //     unsigned int contextCompIdx_;                        //! The component index set for the context
    //     std::array<std::array<const GlobalIndexSet*, numComponents>, numPhases> diffusionVolVarsStencil_;
    //     std::array<std::array<const PositionVector*, numComponents>, numPhases> diffusionVolVarsPositions_;
    //     std::array<std::array<Matrix, numComponents>, numPhases> diffusionT_;
    //     std::array<std::array<Matrix, numComponents>, numPhases> diffusionCA_;
    //     std::array<std::array<Matrix, numComponents>, numPhases> diffusionAB_;
    //     std::array<std::array<Matrix, numComponents>, numPhases> diffusionTout_;

    //     // heat conduction-related variables (see comments above)
    //     const GlobalIndexSet* heatConductionVolVarsStencil_;
    //     const PositionVector* heatConductionVolVarsPositions_;
    //     Matrix heatConductionT_;
    //     Matrix heatConductionCA_;
    //     Matrix heatConductionAB_;
    //     Matrix heatConductionTout_;
    // };

    // //! Specialization for problems involving advection, diffusion and heat conduction
    // template<class TypeTag>
    // struct InteractionVolumeDataHandleImplementation<TypeTag, true, false, true>
    // {
    // private:
    //     using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    //     using GridView = typename GET_PROP_TYPE(TypeTag, GridView);

    //     //! we use the types from the boundary interaction volume traits to be compatible on the boundary
    //     using BoundaryInteractionVolume = typename GET_PROP_TYPE(TypeTag, BoundaryInteractionVolume);
    //     using GlobalIndexSet = typename BoundaryInteractionVolume::Traits::GlobalIndexSet;
    //     using PositionVector = typename BoundaryInteractionVolume::Traits::PositionVector;
    //     using Matrix = typename BoundaryInteractionVolume::Traits::Matrix;
    //     using Vector = typename BoundaryInteractionVolume::Traits::Vector;

    // public:
    //     enum class Contexts : unsigned int
    //     {
    //         undefined,
    //         advection,
    //         heatConduction
    //     };

    //     //! The constructor
    //     InteractionVolumeDataHandleImplementation() : context_(Contexts::undefined) {}

    //     //! sets the context of the cache to advection
    //     void setAdvectionContext()
    //     {
    //         context_ = Contexts::advection;
    //     }

    //     //! sets the context of the cache to advection
    //     void setHeatConductionContext()
    //     {
    //         context_ = Contexts::heatConduction;
    //     }

    //     //! returns the current context
    //     Contexts getContext() const
    //     { return context_; }

    //     //! functions to set the size of the matrices
    //     void resizeT(unsigned int n, unsigned int m)
    //     {
    //         advectionT_.resize(n, m);
    //         heatConductionT_.resize(n, m);
    //     }

    //     void resizeCA(unsigned int n, unsigned int m)
    //     {
    //         advectionCA_.resize(n, m);
    //         heatConductionCA_.resize(n, m);
    //     }

    //     void resizeAB(unsigned int n, unsigned int m)
    //     {
    //         advectionAB_.resize(n, m);
    //         heatConductionAB_.resize(n, m);
    //     }

    //     void resizeOutsideTij(unsigned int n, unsigned int m)
    //     {
    //         advectionTout_.resize(n, m);
    //         heatConductionTout_.resize(n, m);
    //     }

    //     //! functions to set the pointers to stencil and positions
    //     void setVolVarsStencilPointer(const GlobalIndexSet& stencil)
    //     {
    //         assert(context_ != Contexts::undefined && "No valid context set!");
    //         if (context_ == Contexts::advection)
    //             advectionVolVarsStencil_ = &stencil;
    //         else
    //             heatConductionVolVarsStencil_ = &stencil;
    //     }

    //     void setVolVarsPositionsPointer(const PositionVector& pos)
    //     {
    //         assert(context_ != Contexts::undefined && "No valid context set!");
    //         if (context_ == Contexts::advection)
    //             advectionVolVarsPositions_ = &pos;
    //         else
    //             heatConductionVolVarsPositions_ = &pos;
    //     }

    //     //! return functions for the stored data
    //     const GlobalIndexSet& volVarsStencil() const
    //     {
    //         assert(context_ != Contexts::undefined && "No valid context set!");
    //         if (context_ == Contexts::advection)
    //             return *advectionVolVarsStencil_;
    //         else
    //             return *heatConductionVolVarsStencil_;
    //     }

    //     const PositionVector& volVarsPositions() const
    //     {
    //         assert(context_ != Contexts::undefined && "No valid context set!");
    //         if (context_ == Contexts::advection)
    //             return *advectionVolVarsPositions_;
    //         else
    //             return *heatConductionVolVarsPositions_;
    //     }

    //     //! return functions for the stored matrices
    //     const Matrix& T() const
    //     {
    //         assert(context_ != Contexts::undefined && "No valid context set!");
    //         if (context_ == Contexts::advection)
    //             return advectionT_;
    //         else
    //             return heatConductionT_;
    //     }

    //     Matrix& T()
    //     {
    //         assert(context_ != Contexts::undefined && "No valid context set!");
    //         if (context_ == Contexts::advection)
    //             return advectionT_;
    //         else
    //             return heatConductionT_;
    //     }

    //     const Matrix& CA() const
    //     {
    //         assert(context_ != Contexts::undefined && "No valid context set!");
    //         if (context_ == Contexts::advection)
    //             return advectionCA_;
    //         else
    //             return heatConductionCA_;
    //     }

    //     Matrix& CA()
    //     {
    //         assert(context_ != Contexts::undefined && "No valid context set!");
    //         if (context_ == Contexts::advection)
    //             return advectionCA_;
    //         else
    //             return heatConductionCA_;
    //     }

    //     const Matrix& AB() const
    //     {
    //         assert(context_ != Contexts::undefined && "No valid context set!");
    //         if (context_ == Contexts::advection)
    //             return advectionAB_;
    //         else
    //             return heatConductionAB_;
    //     }

    //     Matrix& AB()
    //     {
    //         assert(context_ != Contexts::undefined && "No valid context set!");
    //         if (context_ == Contexts::advection)
    //             return advectionAB_;
    //         else
    //             return heatConductionAB_;
    //     }

    //     const Matrix& outsideTij() const
    //     {
    //         assert(context_ != Contexts::undefined && "No valid context set!");
    //         if (context_ == Contexts::advection)
    //             return advectionTout_;
    //         else
    //             return heatConductionTout_;
    //     }

    //     Matrix& outsideTij()
    //     {
    //         assert(context_ != Contexts::undefined && "No valid context set!");
    //         if (context_ == Contexts::advection)
    //             return advectionTout_;
    //         else
    //             return heatConductionTout_;
    //     }

    // private:
    //     Contexts context_;                                   //! The context variable

    //     // advection-related variables
    //     const GlobalIndexSet* advectionVolVarsStencil_;      //! Pointer to the global volvar indices (stored in the interaction volume)
    //     const PositionVector* advectionVolVarsPositions_;    //! Pointer to the positions of the vol vars (stored in the interaction volume)
    //     Matrix advectionT_;                                  //! The transmissibilities
    //     Matrix advectionCA_;                                 //! The coefficients for the neumann flux transformations
    //     Matrix advectionAB_;                                 //! Coefficients for gradient reconstruction
    //     Matrix advectionTout_;                  //! The transmissibilities associated with "outside" faces (only necessary on surface grids)

    //     // heat conduction-related variables (see comments above)
    //     const GlobalIndexSet* heatConductionVolVarsStencil_;
    //     const PositionVector* heatConductionVolVarsPositions_;
    //     Matrix heatConductionT_;
    //     Matrix heatConductionCA_;
    //     Matrix heatConductionAB_;
    //     Matrix heatConductionTout_;
    // };
} // end namespace

#endif
