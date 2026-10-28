// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_FOEPPL_VON_KARMAN_WRINKLING_TEST_PROBLEM_HH
#define DUMUX_FOEPPL_VON_KARMAN_WRINKLING_TEST_PROBLEM_HH

#include <cmath>
#include <memory>
#include <random>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <dune/common/exceptions.hh>

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/fvproblem.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/math.hh>

namespace Dumux {

//! In-plane loading of the rectangular sheet
enum class SheetLoading { shear, stretch };

/*!
 * \brief Geometry, material and imperfection of a wrinkling rectangular sheet
 *
 * The sheet occupies \f$ [0,L]\times[0,W] \f$. Under `shear` loading the two long
 * edges are held by rigid blocks a fixed distance apart and one of them is moved
 * along \f$ x \f$ by the shear displacement. Under `stretch` loading the two short
 * edges are gripped and pulled apart. In both cases the edges held by the rig are
 * clamped and the remaining two edges are free, in plane and out of plane.
 *
 * In units of \f$ (W, t, t^2/W) \f$ for \f$ (\mathbf{x}, w, \mathbf{u}) \f$ the
 * Föppl-von Kármán system has no parameter but the Poisson ratio, the aspect ratio and the
 * load, so the input is expected in those units: width and thickness one,
 * \f$ E = 12(1-\nu^2) \f$ so that \f$ D = 1 \f$ and \f$ C = 12 \f$, and a shear
 * displacement or strain scaled by \f$ W/t^2 \f$ and \f$ (W/t)^2 \f$. This keeps every
 * field of order one, which the numerically differentiated Jacobian needs, and makes the
 * load read directly as the distance from the buckling threshold.
 *
 * The flat state is an exact solution of the Föppl-von Kármán system, so a small
 * transverse load is applied to select a branch. It is a fixed superposition of the
 * modes admitted by the clamped edges with pseudo-random amplitudes, so that no
 * wavelength is preferred a priori and the wrinkle spacing is selected by the
 * instability itself.
 */
template<class Scalar, class GlobalPosition>
class WrinklingSheet
{
public:
    WrinklingSheet()
    : length_(getParam<Scalar>("Problem.Length"))
    , width_(getParam<Scalar>("Problem.Width"))
    , thickness_(getParam<Scalar>("Problem.Thickness"))
    , poissonRatio_(getParam<Scalar>("Problem.PoissonRatio"))
    , imperfection_(getParam<Scalar>("Problem.Imperfection"))
    , seed_(getParam<Scalar>("Problem.SeedAmplitude", 0.0))
    , seedWavelength_(getParam<Scalar>("Problem.SeedWavelength", 0.0))
    , seedTensionField_(getParam<bool>("Problem.SeedTensionField", false))
    , seedCompressiveStrain_(getParam<Scalar>("Problem.SeedCompressiveStrain", 0.0))
    , damping_(getParam<Scalar>("Problem.Damping", 0.0))
    , gravity_(getParam<Scalar>("Problem.Gravity", 0.0))
    , preStrain_(getParam<Scalar>("Problem.PreStrain", 0.0))
    , numModes_(getParam<int>("Problem.ImperfectionModes", 24))
    , load_(0.0)
    {
        const auto loading = getParam<std::string>("Problem.Loading", "Shear");
        if (loading == "Shear")
            loading_ = SheetLoading::shear;
        else if (loading == "Stretch")
            loading_ = SheetLoading::stretch;
        else
            DUNE_THROW(Dune::IOError, "Problem.Loading must be Shear or Stretch, not " << loading);

        youngsModulus_ = getParam<Scalar>("Problem.E");
        const auto E = youngsModulus_;
        const auto t = thickness_;
        const auto nu = poissonRatio_;
        bendingStiffness_ = E*t*t*t/(12*(1-nu*nu));
        membraneStiffness_ = E*t/(1-nu*nu);
        shearModulus_ = E/(2*(1 + nu));

        std::mt19937 generator(getParam<unsigned>("Problem.ImperfectionSeed", 42));
        std::uniform_real_distribution<Scalar> distribution(-1.0, 1.0);
        modeAmplitudes_.resize(2*numModes_);
        for (auto& c : modeAmplitudes_)
            c = distribution(generator);

        // normalize so that the imperfection is independent of the number of modes
        const auto norm = std::sqrt(0.5*modeAmplitudes_.size());
        for (auto& c : modeAmplitudes_)
            c /= norm;
    }

    SheetLoading loading() const { return loading_; }
    Scalar length() const { return length_; }
    Scalar width() const { return width_; }
    Scalar thickness() const { return thickness_; }

    //! The edges held by the loading rig, which are clamped; the other two are free
    bool isClamped(const GlobalPosition& globalPos) const
    {
        const auto eps = 1e-7*width_;
        return loading_ == SheetLoading::shear
            ? globalPos[1] < eps || globalPos[1] > width_ - eps
            : globalPos[0] < eps || globalPos[0] > length_ - eps;
    }

    //! The second clamp, on which the potential is one unknown constant
    bool isSecondClamp(const GlobalPosition& globalPos) const
    {
        const auto eps = 1e-7*width_;
        return loading_ == SheetLoading::shear ? globalPos[1] > width_ - eps : globalPos[0] > length_ - eps;
    }

    /*!
     * \brief Whether the potential on the second clamp floats
     *
     * The two clamps are separate pieces of the boundary joined by free edges. Prescribing
     * \f$ \varphi = 0 \f$ on both drops the net flux of \f$ \nabla w - \boldsymbol{\theta} \f$
     * through one of them and lets both clamps release the same uniform slope, so on the
     * second clamp \f$ \varphi \f$ is one unknown constant (`Problem.ClampGauge = floating`).
     */
    bool floatingClampGauge() const
    { return getParam<std::string>("Problem.ClampGauge", "floating") == "floating"; }

    //! Shear resultant of the flat sheet in simple shear, \f$ G\gamma t \f$
    Scalar flatShearResultant() const
    { return shearModulus_*(load_/width_)*thickness_; }

    /*!
     * \brief Difference of the principal strains of the flat sheet
     *
     * With \f$ \varepsilon_{xx}=0 \f$, \f$ \varepsilon_{yy}=\varepsilon_0 \f$ and
     * \f$ \varepsilon_{xy}=\gamma/2 \f$ the principal strains are
     * \f$ (\varepsilon_0\pm R)/2 \f$ with \f$ R=\sqrt{\varepsilon_0^2+\gamma^2} \f$.
     */
    Scalar principalStrainDifference() const
    {
        using std::hypot;
        return hypot(preStrain_, load_/width_);
    }

    /*!
     * \brief Angle of the crests to the \f$ x \f$ axis in degrees
     *
     * The wrinkles run along the tensile principal direction of the flat sheet,
     * \f$ \tan2\theta = \gamma/(-\varepsilon_0) \f$, which is \f$ 45^\circ \f$ in simple
     * shear and turns towards \f$ y \f$ as the strain locked in across the sheet grows.
     */
    Scalar tensionFieldCrestAngle() const
    {
        using std::atan2;
        const auto gamma = load_/width_;
        if (gamma == 0.0 && preStrain_ == 0.0)
            return 45.0;
        return 0.5*atan2(gamma, -preStrain_)*180.0/M_PI;
    }

    //! Angle of the wrinkle normal to the \f$ x \f$ axis in degrees
    Scalar tensionFieldNormalAngle() const
    { return tensionFieldCrestAngle() - 90.0; }

    /*!
     * \brief Compressive strain the wrinkles of the developed state absorb
     *
     * The uniaxial tension along the crests contracts the sheet across them by
     * \f$ \nu\varepsilon_1 \f$, so the wrinkles take up the rest of the transverse
     * principal strain, \f$ -\varepsilon_2-\nu\varepsilon_1
     * = [(1-\nu)R-(1+\nu)\varepsilon_0]/2 \f$.
     */
    Scalar tensionFieldCompression() const
    {
        return 0.5*((1.0 - poissonRatio_)*principalStrainDifference() - (1.0 + poissonRatio_)*preStrain_);
    }

    /*!
     * \brief Shear resultant of the wrinkled sheet
     *
     * Once the wrinkles have relieved the compressive principal stress the sheet carries a
     * uniaxial tension along the crests, \f$ N_1 = Et(\varepsilon_0+R)/2 \f$, whose shear
     * component is \f$ (N_1/2)\sin2\theta \f$ with \f$ \sin2\theta = \gamma/R \f$: the
     * flat value times \f$ (1+\nu)/2 \f$ and \f$ 1 + \varepsilon_0/R \f$.
     */
    Scalar tensionFieldShearResultant() const
    {
        const auto R = principalStrainDifference();
        return 0.5*(1.0 + poissonRatio_)*flatShearResultant()*(R > 0.0 ? 1.0 + preStrain_/R : 1.0);
    }

    /*!
     * \brief Damping coefficient of the pseudo-dynamic relaxation
     *
     * Relaxing \f$ c\,\dot{w} + \mathcal{R}(w) = 0 \f$ in pseudo-time descends the energy,
     * so the sheet settles into a stable equilibrium and cannot come to rest on the flat
     * state, which at these loads is a saddle of very high index. A Newton solve has no
     * such preference, which is why it wanders off.
     */
    Scalar damping() const { return damping_; }

    /*!
     * \brief Uniform transverse load per unit area, the weight of a horizontal sheet
     *
     * In the scaled units this is \f$ qW^4/(Dt) \f$; a horizontal Kapton sheet of
     * \f$ 25\,\mu \f$m weighs \f$ 0.35 \f$ Pa, which is \f$ 7.4\times10^5 \f$.
     */
    Scalar gravity() const { return gravity_; }

    /*!
     * \brief Strain across the sheet locked in before the shear is applied
     *
     * The rig tensions the sheet in \f$ y \f$ before the fixed edge is clamped, so the
     * distance between the clamps holds a strain \f$ \varepsilon_0 \f$ under the shear. In
     * the scaled units a physical strain is multiplied by \f$ (W/t)^2 \f$.
     */
    Scalar preStrain() const { return preStrain_; }
    Scalar youngsModulus() const { return youngsModulus_; }
    Scalar bendingStiffness() const { return bendingStiffness_; }
    Scalar membraneStiffness() const { return membraneStiffness_; }
    Scalar shearModulus() const { return shearModulus_; }
    Scalar poissonRatio() const { return poissonRatio_; }

    //! Shear displacement of the moving edge, or applied longitudinal strain
    void setLoad(Scalar load) { load_ = load; }

    //! In-plane displacement of the homogeneous state the rig imposes at the given load
    GlobalPosition affineDisplacement(const GlobalPosition& globalPos, Scalar load) const
    {
        GlobalPosition u(0.0);
        u[0] = loading_ == SheetLoading::shear ? load*globalPos[1]/width_ : load*globalPos[0];
        return u;
    }

    //! In-plane displacement that holds the strain locked in across the sheet
    GlobalPosition preStrainDisplacement(const GlobalPosition& globalPos) const
    {
        GlobalPosition u(0.0);
        if (loading_ == SheetLoading::shear)
            u[1] = preStrain_*globalPos[1];
        return u;
    }
    Scalar load() const { return load_; }

    //! Magnitude of the transverse load that breaks the symmetry of the flat state
    void setImperfection(Scalar imperfection) { imperfection_ = imperfection; }
    Scalar imperfection() const { return imperfection_; }

    int numModes() const { return numModes_; }

    /*!
     * \brief Shape of the imperfection
     *
     * A product of two fixed pseudo-random mode sums, one per direction, which excites
     * every admissible mode pair without preferring a wavelength. The separable form
     * keeps it cheap: this is evaluated in the innermost assembly loop.
     */
    Scalar imperfectionShape(const GlobalPosition& globalPos) const
    {
        using std::sin;
        Scalar sx(0.0), sy(0.0);
        for (int m = 1; m <= numModes_; ++m)
        {
            sx += modeAmplitudes_[m-1]*sin(m*M_PI*globalPos[0]/length_);
            sy += modeAmplitudes_[numModes_ + m-1]*sin(m*M_PI*globalPos[1]/width_);
        }
        return sx*sy;
    }

    //! Transverse imperfection load per unit area
    Scalar imperfectionLoad(const GlobalPosition& globalPos) const
    { return imperfection_*imperfectionShape(globalPos); }

    /*!
     * \brief The wrinkle wave that is seeded: its unit normal, wavenumber and phase
     *
     * Under shear the crests run at 45 degrees, so the wave varies along
     * \f$ (x-y)/\sqrt{2} \f$. Under stretch they run along the sheet and the wave varies
     * across it, where it also has to fit a whole number of half-waves between the
     * clamped long edges.
     */
    std::tuple<GlobalPosition, Scalar, Scalar> seedWave(const GlobalPosition& globalPos) const
    {
        using std::sqrt;
        using std::round;
        using std::max;

        const auto lambda = seedWavelengthOrPredicted_();

        if (loading_ == SheetLoading::shear)
        {
            using std::sin;
            using std::cos;
            const auto crest = tensionFieldCrestAngle()*M_PI/180.0;
            GlobalPosition normal;
            normal[0] = sin(crest);
            normal[1] = -cos(crest);
            const auto k = 2*M_PI/lambda;
            return { normal, k, k*(normal[0]*globalPos[0] + normal[1]*globalPos[1]) };
        }

        GlobalPosition normal;
        normal[0] = 0.0;
        normal[1] = 1.0;
        const auto halfWaves = max(1.0, round(2.0*width_/lambda));
        const auto k = halfWaves*M_PI/width_;
        return { normal, k, k*globalPos[1] };
    }

    //! Envelope that takes the seed smoothly into the clamped edges
    Scalar seedEnvelope(const GlobalPosition& globalPos) const
    {
        using std::sin;
        using std::tanh;
        if (loading_ == SheetLoading::shear)
        {
            // a boundary layer of about half a wavelength with zero slope at the clamped
            // edge, so that the seed carries no strain there that the mesh has to resolve
            const auto layer = seedWavelengthOrPredicted_()/M_PI;
            const auto lower = tanh(globalPos[1]/layer);
            const auto upper = tanh((width_ - globalPos[1])/layer);
            return lower*lower*upper*upper;
        }

        // under stretch only the middle of the sheet is in compression
        return sin(M_PI*globalPos[0]/length_);
    }

    /*!
     * \brief Amplitude of the seeded wrinkle field
     *
     * In the developed state the wrinkles have taken up the compressive strain, which is
     * what relaxes the compression and stops the fine modes being unstable. Averaging
     * \f$ \tfrac{1}{2}\nabla w\otimes\nabla w \f$ over a wave of amplitude \f$ A \f$ and
     * wavenumber \f$ k \f$ gives \f$ A^2k^2/4 \f$ along the direction the wave runs across,
     * so absorbing a compressive strain \f$ \Delta \f$ needs \f$ A = 2\sqrt{\Delta}/k \f$.
     * Under shear \f$ \Delta \f$ is the strain the tension field leaves to the wrinkles;
     * under stretch it has to be supplied, as minus the transverse resultant the flat sheet
     * would carry divided by \f$ Et \f$.
     */
    Scalar seedAmplitude(Scalar wavenumber) const
    {
        if (!seedTensionField_)
            return seed_;

        using std::sqrt;
        using std::max;
        const auto compression = loading_ == SheetLoading::shear
            ? max(0.0, tensionFieldCompression())
            : seedCompressiveStrain_;
        return 2.0*sqrt(compression)/wavenumber;
    }

    /*!
     * \brief Deflection the continuation is started from
     *
     * Deeply post-buckled wrinkling cannot be reached by loading up from the flat
     * state, because every mode between the finest unstable one and the one that is
     * finally selected goes unstable on the way. The continuation is therefore started
     * from a wrinkle field of the expected orientation and wavelength, whose amplitude
     * follows from the excess length the wrinkles have to take up. The envelope makes
     * the seed satisfy the clamped edge conditions.
     */
    Scalar seedDeflection(const GlobalPosition& globalPos) const
    {
        if (seed_ == 0.0 && !seedTensionField_)
            return 0.0;

        using std::sin;
        const auto [normal, k, phase] = seedWave(globalPos);
        return seedAmplitude(k)*seedEnvelope(globalPos)*sin(phase);
    }

    //! Slope of the seeded deflection, which the rotations start from
    GlobalPosition seedRotation(const GlobalPosition& globalPos) const
    {
        GlobalPosition gradient(0.0);
        if (seed_ == 0.0 && !seedTensionField_)
            return gradient;

        const auto h = 1e-6*width_;
        for (int dir = 0; dir < GlobalPosition::dimension; ++dir)
        {
            auto forward = globalPos, backward = globalPos;
            forward[dir] += h;
            backward[dir] -= h;
            gradient[dir] = (seedDeflection(forward) - seedDeflection(backward))/(2.0*h);
        }
        return gradient;
    }

    /*!
     * \brief Shear gradient potential that goes with the seeded wrinkle field
     *
     * The rotation balance ties the potential to the deflection through
     * \f$ \nabla\varphi = -D\nabla\Delta w \f$, so a wave of wavenumber \f$ k \f$
     * carries \f$ \varphi = Dk^2 w \f$. Starting from it leaves the transverse balance
     * with the bending force of the seed to absorb instead of its membrane traction,
     * which is larger by \f$ N/(Dk^2) \f$.
     */
    Scalar seedShearGradPotential(const GlobalPosition& globalPos) const
    {
        if (seed_ == 0.0 && !seedTensionField_)
            return 0.0;

        const auto [normal, k, phase] = seedWave(globalPos);
        return bendingStiffness_*k*k*seedDeflection(globalPos);
    }

    /*!
     * \brief In-plane displacement that goes with the seeded wrinkle field
     *
     * Choosing the amplitude so that the wrinkles absorb the compressive strain only
     * balances the *mean*. The wrinkle strain \f$ \tfrac{1}{2}(\partial_n w)^2 \f$ also
     * oscillates at twice the wavenumber, and the in-plane field has to carry a matching
     * component or the mid-surface is left with an oscillating membrane stress of order
     * \f$ CA^2k^2/4 \f$, which can exceed the applied resultant itself. Cancelling it needs
     * \f$ u_n = -(A^2k/8)\sin(2\mathbf{k}\cdot\mathbf{x}) \f$.
     */
    GlobalPosition seedInPlaneDisplacement(const GlobalPosition& globalPos) const
    {
        GlobalPosition u(0.0);
        if (!seedTensionField_)
            return u;

        using std::sin;
        const auto [normal, k, phase] = seedWave(globalPos);
        const auto amplitude = seedAmplitude(k)*seedEnvelope(globalPos);

        u = normal;
        u *= -amplitude*amplitude*k/8.0*sin(2.0*phase);
        return u;
    }

    /*!
     * \brief Wrinkle wavelength expected from the bending-tension balance
     *
     * \f$ \lambda = 2\pi (B/K)^{1/4} \f$ with the bending stiffness \f$ B \f$ and the
     * stiffness \f$ K = N\pi^2/\ell^2 \f$ of the tension \f$ N \f$ acting along the
     * wrinkles over their length \f$ \ell \f$ \cite CerdaMahadevan2003. Under shear the
     * wrinkles run at 45 degrees, so \f$ \ell = \sqrt{2}W \f$ and the tension along them
     * equals the compression across them, \f$ N = G\gamma t \f$ with the shear strain
     * \f$ \gamma = \delta/W \f$. Under stretch, \f$ \ell = L \f$ and \f$ N = Et\varepsilon \f$.
     */
    Scalar predictedWavelength(Scalar load) const
    {
        using std::pow;
        using std::sqrt;
        const auto tension = loading_ == SheetLoading::shear
            ? shearModulus_*(load/width_)*thickness_
            : youngsModulus_*thickness_*load;
        const auto wrinkleLength = loading_ == SheetLoading::shear
            ? sqrt(2.0)*width_
            : length_;
        return 2*M_PI*pow(bendingStiffness_*wrinkleLength*wrinkleLength/(tension*M_PI*M_PI), 0.25);
    }

private:
    Scalar seedWavelengthOrPredicted_() const
    { return seedWavelength_ > 0.0 ? seedWavelength_ : predictedWavelength(load_); }

    Scalar length_, width_, thickness_, poissonRatio_, imperfection_, seed_, seedWavelength_;
    bool seedTensionField_;
    Scalar seedCompressiveStrain_, damping_, gravity_, preStrain_;
    int numModes_;
    Scalar load_;
    SheetLoading loading_;
    Scalar youngsModulus_, bendingStiffness_, membraneStiffness_, shearModulus_;
    std::vector<Scalar> modeAmplitudes_;
};

/*!
 * \brief Rotation sub-problem of the wrinkling sheet
 */
template<class TypeTag>
class WrinklingTestProblemRotation : public FVProblem<TypeTag>
{
    using ParentType = FVProblem<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    using Sheet = WrinklingSheet<Scalar, GlobalPosition>;

public:
    WrinklingTestProblemRotation(std::shared_ptr<const GridGeometry> gridGeometry,
                                 std::shared_ptr<CouplingManager> couplingManager,
                                 std::shared_ptr<const Sheet> sheet)
    : ParentType(gridGeometry)
    , couplingManager_(couplingManager)
    , sheet_(sheet)
    {}

    Scalar D(const GlobalPosition& globalPos) const { return sheet_->bendingStiffness(); }
    Scalar poissonRatio(const GlobalPosition& globalPos) const { return sheet_->poissonRatio(); }

    BoundaryTypes boundaryTypes(const Element& element, const SubControlVolume& scv) const
    {
        BoundaryTypes values;
        if (sheet_->isClamped(scv.dofPosition()))
            values.setAllDirichlet();
        else
            values.setAllNeumann();
        return values;
    }

    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

    //! The traction \f$ \mathbf{T}\mathbf{n} = -\varphi\mathbf{n} \f$ of a free edge
    template<class ElementVolumeVariables, class ElementFluxVariablesCache>
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVariablesCache& elemFluxVarsCache,
                        const SubControlVolumeFace& scvf) const
    {
        NumEqVector values(0.0);
        const auto vars = this->couplingManager().deformationAndPotentials(fvGeometry, scvf);
        values.axpy(vars[this->couplingManager().shearGradPotentialIdx()], scvf.unitOuterNormal());
        return values;
    }

    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

private:
    std::shared_ptr<CouplingManager> couplingManager_;
    std::shared_ptr<const Sheet> sheet_;
};

/*!
 * \brief Deformation-and-potentials sub-problem of the wrinkling sheet
 *
 * On the edges held by the rig the deflection and the shear gradient potential vanish, the
 * latter fixing the gauge of the Helmholtz decomposition. The other two edges are free, which
 * in the mixed formulation is a natural condition on all three equations; the shear curl
 * potential is then determined by the free-edge traction of the rotations and needs no
 * interior constraint.
 */
template<class TypeTag>
class WrinklingTestProblemDeformation : public FVProblem<TypeTag>
{
    using ParentType = FVProblem<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    using Sheet = WrinklingSheet<Scalar, GlobalPosition>;

public:
    WrinklingTestProblemDeformation(std::shared_ptr<const GridGeometry> gridGeometry,
                                    std::shared_ptr<CouplingManager> couplingManager,
                                    std::shared_ptr<const Sheet> sheet)
    : ParentType(gridGeometry)
    , couplingManager_(couplingManager)
    , sheet_(sheet)
    , floatingClampGauge_(sheet->floatingClampGauge())
    {}

    Scalar damping(const GlobalPosition& globalPos) const { return sheet_->damping(); }

    Scalar membraneStiffness(const GlobalPosition& globalPos) const { return sheet_->membraneStiffness(); }
    Scalar poissonRatio(const GlobalPosition& globalPos) const { return sheet_->poissonRatio(); }

    BoundaryTypes boundaryTypes(const Element& element, const SubControlVolume& scv) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        const auto& pos = scv.dofPosition();
        if (sheet_->isClamped(pos))
        {
            values.setDirichlet(Indices::verticalDeformationIdx, Indices::shearGradPotentialEqIdx);
            if (!sheet_->isSecondClamp(pos) || !floatingClampGauge_)
                values.setDirichlet(Indices::shearGradPotentialIdx, Indices::deformationEqIdx);
        }
        return values;
    }

    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

    template<class ElementVolumeVariables, class ElementFluxVariablesCache>
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVariablesCache& elemFluxVarsCache,
                        const SubControlVolumeFace& scvf) const
    {
        NumEqVector values(0.0);
        const auto rotation = this->couplingManager().rotation(fvGeometry, scvf);
        const auto tangent = [&](){
            auto tangent = scvf.unitOuterNormal();
            std::swap(tangent[0], tangent[1]);
            tangent[1] = -tangent[1];
            return tangent;
        }();
        values[Indices::shearCurlPotentialEqIdx] = -vtmv(tangent, 1.0, rotation);
        return values;
    }

    NumEqVector sourceAtPos(const GlobalPosition& globalPos) const
    {
        NumEqVector source(0.0);
        source[Indices::shearGradPotentialEqIdx] = sheet_->imperfectionLoad(globalPos) + sheet_->gravity();
        return source;
    }

    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values(0.0);
        values[Indices::verticalDeformationIdx] = sheet_->seedDeflection(globalPos);
        values[Indices::shearGradPotentialIdx] = sheet_->seedShearGradPotential(globalPos);
        return values;
    }

    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

private:
    std::shared_ptr<CouplingManager> couplingManager_;
    std::shared_ptr<const Sheet> sheet_;
    bool floatingClampGauge_;
};

/*!
 * \brief In-plane sub-problem of the wrinkling sheet
 *
 * Under shear the two long edges are held by the rigid blocks of the shear rig, one
 * of which is moved along \f$ x \f$; the two short edges are traction free. Under
 * stretch the two short edges are gripped, which prevents the lateral Poisson
 * contraction in their vicinity and is what puts the centre of the sheet into
 * transverse compression, and the two long edges are traction free.
 */
template<class TypeTag>
class WrinklingTestProblemInPlane : public FVProblem<TypeTag>
{
    using ParentType = FVProblem<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    using Sheet = WrinklingSheet<Scalar, GlobalPosition>;

public:
    WrinklingTestProblemInPlane(std::shared_ptr<const GridGeometry> gridGeometry,
                                std::shared_ptr<CouplingManager> couplingManager,
                                std::shared_ptr<const Sheet> sheet)
    : ParentType(gridGeometry)
    , couplingManager_(couplingManager)
    , sheet_(sheet)
    , eps_(1e-7*sheet->width())
    {}

    Scalar damping(const GlobalPosition& globalPos) const { return sheet_->damping(); }

    Scalar membraneStiffness(const GlobalPosition& globalPos) const { return sheet_->membraneStiffness(); }
    Scalar poissonRatio(const GlobalPosition& globalPos) const { return sheet_->poissonRatio(); }

    BoundaryTypes boundaryTypes(const Element& element, const SubControlVolume& scv) const
    {
        BoundaryTypes values;
        if (sheet_->isClamped(scv.dofPosition()))
            values.setAllDirichlet();
        else
            values.setAllNeumann();
        return values;
    }

    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values(0.0);
        if (sheet_->loading() == SheetLoading::shear)
        {
            if (globalPos[1] > sheet_->width() - eps_)
            {
                values[0] = sheet_->load();
                values[1] = sheet_->preStrain()*sheet_->width();
            }
        }
        else if (globalPos[0] > sheet_->length() - eps_)
            values[0] = sheet_->load()*sheet_->length();
        return values;
    }

    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values(0.0);
        const auto affine = sheet_->affineDisplacement(globalPos, sheet_->load());
        const auto locked = sheet_->preStrainDisplacement(globalPos);
        const auto seed = sheet_->seedInPlaneDisplacement(globalPos);
        for (int dir = 0; dir < GlobalPosition::dimension; ++dir)
            values[dir] = affine[dir] + locked[dir] + seed[dir];

        return values;
    }

    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

private:
    std::shared_ptr<CouplingManager> couplingManager_;
    std::shared_ptr<const Sheet> sheet_;
    Scalar eps_;
};

} // end namespace Dumux

#endif
