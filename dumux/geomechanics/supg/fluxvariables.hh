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
//STOKES
/*!
 * \file
 * \brief This file contains the data which is required to calculate
 *        the fluxes of the Stokes model over a face of a finite volume.
 *
 * This means pressure gradients, phase densities at the integration point, etc.
 */
#ifndef DUMUX_STOKES_FLUX_VARIABLES_HH
#define DUMUX_STOKES_FLUX_VARIABLES_HH

#include <dune/common/exceptions.hh>

#include <dumux/common/math.hh>
#include <dumux/common/valgrind.hh>

#include "properties.hh"

namespace Dumux
{

/*!
 * \ingroup BoxStokesModel
 * \ingroup ImplicitFluxVariables
 * \brief This template class contains the data which is required to
 *        calculate the mass and momentum fluxes over the face of a
 *        sub-control volume for the Stokes model.
 *
 * This means pressure gradients, phase densities, viscosities, etc.
 * at the integration point of the sub-control-volume face.
 */
template <class TypeTag>
class StokesFluxVariables
{
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) Implementation;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;

    enum { dim = GridView::dimension };

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef Dune::FieldVector<Scalar, dim> DimVector;
    typedef Dune::FieldMatrix<Scalar, dim, dim> DimMatrix;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename FVElementGeometry::SubControlVolumeFace SCVFace;

public:
    /*!
     * \brief Compute / update the flux variables
     *
     * \param problem The problem
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param fIdx The local index of the SCV (sub-control-volume) face
     * \param elemVolVars The volume variables of the current element
     * \param onBoundary A boolean variable to specify whether the flux variables
     * are calculated for interior SCV faces or boundary faces, default=false
     * \todo The fvGeometry should be better initialized, passed and stored as an std::shared_ptr
     */
    void update(const Problem &problem,
                const Element &element,
                const FVElementGeometry &fvGeometry,
                const int fIdx,
                const ElementVolumeVariables &elemVolVars,
                const bool onBoundary = false)
    {
        fvGeometryPtr_ = &fvGeometry;
        onBoundary_ = onBoundary;
        fIdx_ = fIdx;

        asImp_().calculateValues_(problem, element, elemVolVars);
       // asImp_().determineUpwindDirection_(elemVolVars);
    }


//ELASTIC
    /*!
     * \brief Returns the sub-control-volume face.
     */
     const SCVFace &face() const
     { return fvGeometry_().subContVolFace[faceIdx_]; }

//STOKES
protected:
    //! Returns the implementation of the flux variables (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

    void calculateValues_(const Problem &problem,
                          const Element &element,
                          const ElementVolumeVariables &elemVolVars)
    {
        // calculate gradients and secondary variables at IPs
        DimVector tmp(0.0);

        density_ = Scalar(0);
        dynamicViscosity_ = Scalar(0);
        pressure_ = Scalar(0);
        normalvelocity_ = Scalar(0);
        velocity_ = Scalar(0);
        pressureGrad_ = Scalar(0);
        velocityGrad_ = Scalar(0);

        for (int scvIdx = 0;
             scvIdx < fvGeometry_().numScv;
             scvIdx++) // loop over adjacent vertices
        {
            // phase density and viscosity at IP
            density_ += elemVolVars[scvIdx].density() *
                face().shapeValue[scvIdx];
            dynamicViscosity_ += elemVolVars[scvIdx].dynamicViscosity() *
                face().shapeValue[scvIdx];
            pressure_ += elemVolVars[scvIdx].pressure() *
                face().shapeValue[scvIdx];

            // velocity at the IP (fluxes)
            DimVector velocityTimesShapeValue = elemVolVars[scvIdx].velocity();
            velocityTimesShapeValue *= face().shapeValue[scvIdx];
            velocity_ += velocityTimesShapeValue;

            // the pressure gradient
            tmp = face().grad[scvIdx];
            tmp *= elemVolVars[scvIdx].pressure();
            pressureGrad_ += tmp;
            // take gravity into account
            tmp = problem.gravity();
            tmp *= density_;
            // pressure gradient including influence of gravity
            pressureGrad_ -= tmp;

            // the velocity gradients and divergence
            for (int dimIdx = 0; dimIdx < dim; ++dimIdx)
            {
                tmp = face().grad[scvIdx];
                tmp *= elemVolVars[scvIdx].velocity()[dimIdx];
                velocityGrad_[dimIdx] += tmp;
            }
        }

        normalvelocity_ = velocity_ * face().normal;

        Valgrind::CheckDefined(density_);
        Valgrind::CheckDefined(dynamicViscosity_);
        Valgrind::CheckDefined(normalvelocity_);
        Valgrind::CheckDefined(velocity_);
        Valgrind::CheckDefined(pressureGrad_);
        Valgrind::CheckDefined(velocityGrad_);
    }


public:
    /*!
     * \brief Return the pressure \f$\mathrm{[Pa]}\f$ at the integration
     *        point.
     */
    Scalar pressure() const
    { return pressure_; }

    /*!
     * \brief Return the mass density \f$ \mathrm{[kg/m^3]} \f$ at the integration
     *        point.
     */
    Scalar density() const
    { return density_; }

    /*!
     * \brief Return the dynamic viscosity \f$ \mathrm{[Pa\cdot s]} \f$ at the integration
     *        point.
     */
    Scalar dynamicViscosity() const
    { return dynamicViscosity_; }

    /*!
     * \brief Returns the kinematic viscosity \f$ \frac{m^2}{s} \f$ of the fluid in
     *        the sub-control volume.
     */
    Scalar kinematicViscosity() const
    { return dynamicViscosity_ / density_; }

    /*!
     * \brief Return the velocity \f$ \mathrm{[m/s]} \f$ at the integration
     *        point multiplied by the normal and the area.
     */
    Scalar normalVelocity() const
    { return normalvelocity_; }

    /*!
     * \brief Return the pressure gradient at the integration point.
     */
    const DimVector &pressureGrad() const
    { return pressureGrad_; }

    /*!
     * \brief Return the velocity vector at the integration point.
     */
    const DimVector &velocity() const
    { return velocity_; }

    /*!
     * \brief Return the velocity gradient at the integration
     *        point of a face.
     */
    const DimMatrix &velocityGrad() const
    { return velocityGrad_; }

    /*!
     * \brief Return the dynamic eddy viscosity
     *        \f$\mathrm{[Pa \cdot s]} = \mathrm{[N \cdot s/m^2]}\f$ (if implemented).
     */
    const Scalar dynamicEddyViscosity() const
    { return kinematicEddyViscosity() * density(); }

    /*!
     * \brief Return the kinematic eddy viscosity
     *        \f$\mathrm{[m^2/s]}\f$ (if implemented).
     */
    const Scalar kinematicEddyViscosity() const
    { return 0; }


protected:
    // return const reference to fvGeometry
    const FVElementGeometry& fvGeometry_() const
    { return *fvGeometryPtr_; }

//STOKES: fIdx_ instead of faceIdx_ as in ELASTIC
    // the index of the considered face
    int fIdx_;
    bool onBoundary_;

    // values at the integration point
    Scalar density_;
    Scalar dynamicViscosity_;
    Scalar pressure_;
    Scalar normalvelocity_;
    DimVector velocity_;

    // gradients at the IPs
    DimVector pressureGrad_;
    DimMatrix velocityGrad_;



private:
    const FVElementGeometry* fvGeometryPtr_; //!< Information about the geometry of discretization
};

} // end namespace

#endif








 /*!
 * \file
 *
 * \brief This file contains the data which is required to calculate the gradients
 *        over a face of a finite volume that are needed for the momentum balance
 *        of a linear-elastic solid.
 *
 * This means gradients of solid displacement vectors, strains and stresses at
 * the integration point
 *
 * This class is also used as a base class for the one-phase and two-phase
 * linear-elastic models.
 */
#ifndef DUMUX_ELASTIC_FLUX_VARIABLES_HH
#define DUMUX_ELASTIC_FLUX_VARIABLES_HH

#include "properties.hh"

namespace Dumux
{

/*!
 * \ingroup ElasticBoxModel
 * \ingroup ImplicitFluxVariables
 * \brief This template class contains the data which is required to
 *        calculate the gradients over a face of a finite volume for
 *        the linear elasticity model.
 *
 * This means gradients of solid displacement vectors, strains and stresses at
 * the integration point
 */
template<class TypeTag>
class ElasticFluxVariablesBase
{
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables) ElementVolumeVariables;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity Element;
    enum {
        dim = GridView::dimension
    };

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::FieldVector<Scalar, dim> DimVector;
    typedef Dune::FieldMatrix<Scalar, dim, dim> DimMatrix;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename FVElementGeometry::SubControlVolumeFace SCVFace;

public:
    /*!
     * \brief Compute / update the flux variables
     *
     * \param problem The problem
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry
     * \param fIdx The local index of the SCV (sub-control-volume) face
     * \param elemVolVars The volume variables of the current element
     * \param onBoundary A boolean variable to specify whether the flux variables
     * are calculated for interior SCV faces or boundary faces, default=false
     * \todo The fvGeometry should be better initialized, passed and stored as an std::shared_ptr
     */
    void update(const Problem &problem,
                const Element &element,
                const FVElementGeometry &fvGeometry,
                const int fIdx,
                const ElementVolumeVariables &elemVolVars,
                const bool onBoundary = false)
    {
        fvGeometryPtr_ = &fvGeometry;
        onBoundary_ = onBoundary;
        faceIdx_ = fIdx;

        gradU_ = 0.0;
        gradUTransposed_ = 0.0;
        epsilon_ = 0.0;
        sigma_ = 0.0;

        lambda_ = 0.0;
        mu_ = 0.0;
        divU_ = 0.0;

        calculateGradients_(problem, element, elemVolVars);
        calculateStrain_(problem, element, elemVolVars);
        calculateStress_(problem, element, elemVolVars);
    }

        /*!
         * \brief Return a stress tensor component [Pa] at the integration point.
         */
        Scalar sigma(int row, int col) const
        { return sigma_[row][col]; }

        /*!
         * \brief Return the stress tensor [Pa] at the integration point.
         */
        DimMatrix sigma() const
        { return sigma_; }

        /*!
         * \brief Return the volumetric strain i.e. the divergence of the solid displacement
         *  vector at the integration point.
         */
        Scalar divU() const
        { return divU_; }

        /*!
          * \brief Returns the Lame parameter mu at integration point.
          */
        Scalar mu() const
        { return mu_; }

        /*!
         * \brief Returns the pressure value at integration point.
         */
        Scalar pressure() const
        { return pressure_; }

        /*!
          * \brief Returns the sub-control-volume face.
          */
        const SCVFace &face() const
        { return fvGeometry_().subContVolFace[faceIdx_]; }

protected:
    /*!
     * \brief Calculation of the solid displacement gradients.
     *
     *        \param problem The considered problem file
     *        \param element The considered element of the grid
     *        \param elemVolVars The parameters stored in the considered element
     */
    void calculateGradients_(const Problem &problem,
            const Element &element,
            const ElementVolumeVariables &elemVolVars)
    {
        const VolumeVariables &volVarsI = elemVolVars[face().i];
        const VolumeVariables &volVarsJ = elemVolVars[face().j];

        // calculate gradients
        DimVector tmp(0.0);
        for (int idx = 0;
                idx < fvGeometry_().numScv;
                idx++) // loop over adjacent vertices
        {
            // FE gradient at vertex idx
            const DimVector &feGrad = face().grad[idx];

            // the displacement vector gradient
            for (int coordIdx = 0; coordIdx < dim; ++coordIdx) {
                tmp = feGrad;
                tmp *= elemVolVars[idx].displacement(coordIdx);
                gradU_[coordIdx] += tmp;
            }
        }

        // average the Lame parameters at integration point
        // note: it still needs to be checked which mean (arithmetic, harmonic.. is appropriate
      //  lambda_ = (volVarsI.lambda() + volVarsJ.lambda()) / 2.;
      //  mu_ = (volVarsI.mu() + volVarsJ.mu()) / 2.;

        for(int col=0; col < dim; col++)
        {
         //?   divU_ += gradU_[col][col];  //needed for divergence of stress tensor?

            for(int row=0; row<dim; row++)
                gradUTransposed_[row][col] = gradU_[col][row];
        }
    }

    /*!
     * \brief Calculation of the strain tensor.
     *
     *        \param problem The considered problem file
     *        \param element The considered element of the grid
     *        \param elemVolVars The parameters stored in the considered element
     */
    void calculateStrain_(const Problem &problem,
            const Element &element,
            const ElementVolumeVariables &elemVolVars)
    {
        // calculate the strain tensor
        epsilon_ += gradU_;
        epsilon_ += gradUTransposed_;
        epsilon_ *= 0.5;
    }

    /*!
     * \brief Calculation of the stress tensor.
     *
     *        \param problem The considered problem file
     *        \param element The considered element of the grid
     *        \param elemVolVars The parameters stored in the considered element
     */
    void calculateStress_(const Problem &problem,
            const Element &element,
            const ElementVolumeVariables &elemVolVars)
    {
        DimMatrix firstTerm(0.0), secondTerm(0.0);

        epsilonTimesIdentity_ = divU_;

        firstTerm += epsilon_;
        firstTerm *= 2;
        firstTerm *= mu_;

        for (int i = 0; i < dim; ++i)
        secondTerm[i][i] += 1.0;
        secondTerm *= -pressure_;
       /& secondTerm *= epsilonTimesIdentity_;

        // calculate the stress tensor
        sigma_ += firstTerm;
        sigma_ += secondTerm;
    }

    // return const reference to fvGeometry
    const FVElementGeometry& fvGeometry_() const
    { return *fvGeometryPtr_; }

    int faceIdx_;
    bool onBoundary_;

    // Lame parameter mu at the integration point
    Scalar mu_;
    // Lame parameter lambda at the integration point
    //Scalar lambda_;
    // divergence of the solid displacement vector at the integration point
    Scalar divU_;
    // volumetric strain at the integration point
    Scalar epsilonTimesIdentity_;
    // gradient and transposed gradient of the solid displacement vector
    // at the integration point
    DimMatrix gradU_, gradUTransposed_;
    // pressure value at the integration point
    DimMatrix pressure_;
    // strain tensor at the integration point
    DimMatrix epsilon_;
    // stress tensor at the integration point
    DimMatrix sigma_;

private:
    const FVElementGeometry* fvGeometryPtr_; //!< Information about the geometry of discretization

};

} // end namespace

#endif
