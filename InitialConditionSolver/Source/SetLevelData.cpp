/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "SetLevelData.H"
#include "AMRIO.H"
#include "BCFunc.H"
#include "BRMeshRefine.H"
#include "BiCGStabSolver.H"
#include "BoxIterator.H"
#include "CONSTANTS.H"
#include "CoarseAverage.H"
#include "LoadBalance.H"
#include "MyPhiFunction.H"
#include "PoissonParameters.H"
#include "SetBinaryBH.H"
#include "VariableCoeffPoissonOperatorFactory.H"
#include "computeNorm.H"
#include "parstream.H"
#include <cmath>

// Set various LevelData functions across the grid

// This takes an IntVect and writes the physical coordinates to a RealVect
inline void get_loc(RealVect &a_out_loc, const IntVect &a_iv,
                    const RealVect &a_dx, const PoissonParameters &a_params)
{
    a_out_loc = a_iv + 0.5 * RealVect::Unit;
    a_out_loc *= a_dx;
    a_out_loc -= a_params.domainLength / 2.0;
}

// set initial guess value for the conformal factor psi
// defined by \gamma_ij = \psi^4 \tilde \gamma_ij, scalar field phi
// and \bar Aij = psi^2 A_ij.
// For now the default setup is 2 Bowen York BHs plus a scalar field
// with some initial user specified configuration
void set_initial_conditions(LevelData<FArrayBox> &a_multigrid_vars,
                            LevelData<FArrayBox> &a_dpsi, const RealVect &a_dx,
                            const PoissonParameters &a_params)
{

    CH_assert(a_multigrid_vars.nComp() == NUM_MULTIGRID_VARS);

    DataIterator dit = a_multigrid_vars.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
        FArrayBox &multigrid_vars_box = a_multigrid_vars[dit()];
        FArrayBox &dpsi_box = a_dpsi[dit()];
        Box b = multigrid_vars_box.box();
        BoxIterator bit(b);
        for (bit.begin(); bit.ok(); ++bit)
        {

            // work out location on the grid
            IntVect iv = bit();

            // set psi to 1.0 and zero dpsi
            // note that we don't include the singular part of psi
            // for the BHs - this is added at the output data stage
            // and when we calculate psi_0 in the rhs etc
            // as it already satisfies Laplacian(psi) = 0
            multigrid_vars_box(iv, c_psi_reg) = 1.0;
            dpsi_box(iv, 0) = 0.0;

            // set the phi value - need the distance from centre
            RealVect loc;
            get_loc(loc, iv, a_dx, a_params);

            // set phi according to user defined function
            multigrid_vars_box(iv, c_phi_0) =
                my_phi_function(loc, a_params.phi_amplitude,
                                a_params.phi_wavelength, a_params.domainLength);

            // set Aij for spin and momentum according to BH params
            set_binary_bh_Aij(multigrid_vars_box, iv, loc, a_params);
        }
    }
} // end set_initial_conditions

// computes the Laplacian of psi at a point in a box
inline Real get_laplacian_psi(const IntVect &a_iv, const FArrayBox &a_psi_fab,
                              const RealVect &a_dx)
{
    Real laplacian_of_psi = 0.0;
    for (int idir = 0; idir < SpaceDim; ++idir)
    {
        IntVect iv_offset1 = a_iv;
        IntVect iv_offset2 = a_iv;
        iv_offset1[idir] -= 1;
        iv_offset2[idir] += 1;

        // 2nd order stencil for now
        Real d2psi_dxdx = 1.0 / (a_dx[idir] * a_dx[idir]) *
                          (+1.0 * a_psi_fab(iv_offset2) -
                           2.0 * a_psi_fab(a_iv) + 1.0 * a_psi_fab(iv_offset1));
        laplacian_of_psi += d2psi_dxdx;
    }
    return laplacian_of_psi;
} // end get_laplacian_psi

// computes the gradient of the scalar field squared at a point in a box
// i.e. \delta^{ij} d_i phi d_j phi
inline Real get_grad_phi_sq(const IntVect &a_iv, const FArrayBox &a_phi_fab,
                            const RealVect &a_dx)
{
    Real grad_phi_sq = 0.0;
    for (int idir = 0; idir < SpaceDim; ++idir)
    {
        IntVect iv_offset1 = a_iv;
        IntVect iv_offset2 = a_iv;
        iv_offset1[idir] -= 1;
        iv_offset2[idir] += 1;

        // 2nd order stencils for now
        Real dphi_dx =
            0.5 * (a_phi_fab(iv_offset2) - a_phi_fab(iv_offset1)) / a_dx[idir];

        grad_phi_sq += dphi_dx * dphi_dx;
    }
    return grad_phi_sq;
} // end get_grad_phi_sq

// set the rhs source for the poisson eqn
void set_rhs(LevelData<FArrayBox> &a_rhs,
             LevelData<FArrayBox> &a_multigrid_vars, const RealVect &a_dx,
             const PoissonParameters &a_params, const Real constant_K)
{

    CH_assert(a_multigrid_vars.nComp() == NUM_MULTIGRID_VARS);

    DataIterator dit = a_rhs.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
        FArrayBox &multigrid_vars_box = a_multigrid_vars[dit()];
        FArrayBox &rhs_box = a_rhs[dit()];
        rhs_box.setVal(0.0, 0);
        Box this_box = rhs_box.box(); // no ghost cells

        FArrayBox phi_fab(Interval(c_phi_0, c_phi_0), multigrid_vars_box);
        FArrayBox psi_fab(Interval(c_psi_reg, c_psi_reg), multigrid_vars_box);

        BoxIterator bit(this_box);
        for (bit.begin(); bit.ok(); ++bit)
        {
            IntVect iv = bit();
            RealVect loc;
            get_loc(loc, iv, a_dx, a_params);

            // rhs = m/8 psi_0^5 - 2 pi rho_grad psi_0  - laplacian(psi_0)
            Real m =
                get_m(multigrid_vars_box(iv, c_phi_0), a_params, constant_K);
            Real grad_phi_sq = get_grad_phi_sq(iv, phi_fab, a_dx);
            Real laplacian_of_psi = get_laplacian_psi(iv, psi_fab, a_dx);

            // Also \bar  A_ij \bar A^ij
            Real A2 = 0.0;
            A2 = pow(multigrid_vars_box(iv, c_A11_0), 2.0) +
                 pow(multigrid_vars_box(iv, c_A22_0), 2.0) +
                 pow(multigrid_vars_box(iv, c_A33_0), 2.0) +
                 2 * pow(multigrid_vars_box(iv, c_A12_0), 2.0) +
                 2 * pow(multigrid_vars_box(iv, c_A13_0), 2.0) +
                 2 * pow(multigrid_vars_box(iv, c_A23_0), 2.0);

            Real psi_bl = get_psi_brill_lindquist(loc, a_params);
            Real psi_0 = multigrid_vars_box(iv, c_psi_reg) + psi_bl;

            rhs_box(iv, 0) = 0.125 * m * pow(psi_0, 5.0) -
                             0.125 * A2 * pow(psi_0, -7.0) -
                             M_PI * a_params.G_Newton * grad_phi_sq * psi_0 -
                             laplacian_of_psi;
        }
    }
} // end set_rhs

// Set the integrand for the integrability condition for constant K
// when periodic BCs set
void set_constant_K_integrand(LevelData<FArrayBox> &a_integrand,
                              LevelData<FArrayBox> &a_multigrid_vars,
                              const RealVect &a_dx,
                              const PoissonParameters &a_params)
{

    CH_assert(a_multigrid_vars.nComp() == NUM_MULTIGRID_VARS);

    DataIterator dit = a_integrand.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
        FArrayBox &multigrid_vars_box = a_multigrid_vars[dit()];
        FArrayBox &integrand_box = a_integrand[dit()];
        integrand_box.setVal(0.0, 0);
        Box this_box = integrand_box.box(); // no ghost cells

        FArrayBox phi_fab(Interval(c_phi_0, c_phi_0), multigrid_vars_box);
        FArrayBox psi_fab(Interval(c_psi_reg, c_psi_reg), multigrid_vars_box);

        BoxIterator bit(this_box);
        for (bit.begin(); bit.ok(); ++bit)
        {
            IntVect iv = bit();
            RealVect loc;
            get_loc(loc, iv, a_dx, a_params);

            // integrand = -1.5*m + 1.5 * \bar A_ij \bar A^ij psi_0^-12 +
            // 24 pi rho_grad psi_0^-4  + 12*laplacian(psi_0)*psi^-5
            Real m = get_m(multigrid_vars_box(iv, c_phi_0), a_params, 0.0);
            Real grad_phi_sq = get_grad_phi_sq(iv, phi_fab, a_dx);
            Real laplacian_of_psi = get_laplacian_psi(iv, psi_fab, a_dx);

            // Also \bar  A_ij \bar A^ij
            Real A2 = 0.0;
            A2 = pow(multigrid_vars_box(iv, c_A11_0), 2.0) +
                 pow(multigrid_vars_box(iv, c_A22_0), 2.0) +
                 pow(multigrid_vars_box(iv, c_A33_0), 2.0) +
                 2 * pow(multigrid_vars_box(iv, c_A12_0), 2.0) +
                 2 * pow(multigrid_vars_box(iv, c_A13_0), 2.0) +
                 2 * pow(multigrid_vars_box(iv, c_A23_0), 2.0);

            Real psi_bl = get_psi_brill_lindquist(loc, a_params);
            Real psi_0 = multigrid_vars_box(iv, c_psi_reg) + psi_bl;

            integrand_box(iv, 0) = -1.5 * m + 1.5 * A2 * pow(psi_0, -12.0) +
                                   12.0 * M_PI * a_params.G_Newton *
                                       grad_phi_sq * pow(psi_0, -4.0) +
                                   12.0 * laplacian_of_psi * pow(psi_0, -5.0);
        }
    }
} // end set_constant_K_integrand

// set the regrid condition - abs value of this drives AMR regrid
void set_regrid_condition(LevelData<FArrayBox> &a_condition,
                          LevelData<FArrayBox> &a_multigrid_vars,
                          const RealVect &a_dx,
                          const PoissonParameters &a_params)
{

    CH_assert(a_multigrid_vars.nComp() == NUM_MULTIGRID_VARS);

    DataIterator dit = a_condition.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
        FArrayBox &multigrid_vars_box = a_multigrid_vars[dit()];
        FArrayBox &condition_box = a_condition[dit()];
        condition_box.setVal(0.0, 0);
        Box this_box = condition_box.box(); // no ghost cells

        FArrayBox phi_fab(Interval(c_phi_0, c_phi_0), multigrid_vars_box);

        BoxIterator bit(this_box);
        for (bit.begin(); bit.ok(); ++bit)
        {
            IntVect iv = bit();
            RealVect loc;
            get_loc(loc, iv, a_dx, a_params);

            // calculate contributions
            Real m = get_m(multigrid_vars_box(iv, c_phi_0), a_params, 0.0);
            Real grad_phi_sq = get_grad_phi_sq(iv, phi_fab, a_dx);

            // Also \bar  A_ij \bar A^ij
            Real A2 = 0.0;
            A2 = pow(multigrid_vars_box(iv, c_A11_0), 2.0) +
                 pow(multigrid_vars_box(iv, c_A22_0), 2.0) +
                 pow(multigrid_vars_box(iv, c_A33_0), 2.0) +
                 2 * pow(multigrid_vars_box(iv, c_A12_0), 2.0) +
                 2 * pow(multigrid_vars_box(iv, c_A13_0), 2.0) +
                 2 * pow(multigrid_vars_box(iv, c_A23_0), 2.0);

            // the condition is similar to the rhs but we take abs
            // value of the contributions and add in BH criteria
            Real psi_bl = get_psi_brill_lindquist(loc, a_params);
            Real psi_0 = multigrid_vars_box(iv, c_psi_reg) + psi_bl;
            condition_box(iv, 0) = 1.5 * abs(m) + 1.5 * A2 * pow(psi_0, -7.0) +
                                   12.0 * M_PI * a_params.G_Newton *
                                       abs(grad_phi_sq) * pow(psi_0, 1.0) +
                                   log(psi_0);
        }
    }
} // end set_regrid_condition

// Add the correction to psi0 after the solver operates
void set_update_psi0(LevelData<FArrayBox> &a_multigrid_vars,
                     LevelData<FArrayBox> &a_dpsi,
                     const Copier &a_exchange_copier)
{

    // first exchange ghost cells for dpsi so they are filled with the correct
    // values
    a_dpsi.exchange(a_dpsi.interval(), a_exchange_copier);

    DataIterator dit = a_multigrid_vars.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
        FArrayBox &multigrid_vars_box = a_multigrid_vars[dit()];
        FArrayBox &dpsi_box = a_dpsi[dit()];

        Box this_box = multigrid_vars_box.box();
        BoxIterator bit(this_box);
        for (bit.begin(); bit.ok(); ++bit)
        {
            IntVect iv = bit();
            multigrid_vars_box(iv, c_psi_reg) += dpsi_box(iv, 0);
        }
    }
}

// m(K, rho) = 2/3K^2 - 16piG rho
inline Real get_m(const Real &phi_here, const PoissonParameters &a_params,
                  const Real constant_K)
{

    // KC TODO:
    // For now rho is just the gradient term which is kept separate
    // ... may want to add V(phi) and phidot/Pi here later though
    Real Pi_field = 0.0;
    Real V_of_phi = 0.0;
    Real rho = 0.5 * Pi_field * Pi_field + V_of_phi;

    return (2.0 / 3.0) * (constant_K * constant_K) -
           16.0 * M_PI * a_params.G_Newton * rho;
}

// The coefficient of the I operator on dpsi
void set_a_coef(LevelData<FArrayBox> &a_aCoef,
                LevelData<FArrayBox> &a_multigrid_vars,
                const PoissonParameters &a_params, const RealVect &a_dx,
                const Real constant_K)
{

    CH_assert(a_multigrid_vars.nComp() == NUM_MULTIGRID_VARS);

    DataIterator dit = a_aCoef.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
        FArrayBox &aCoef_box = a_aCoef[dit()];
        FArrayBox &multigrid_vars_box = a_multigrid_vars[dit()];
        Box this_box = aCoef_box.box();

        FArrayBox phi_fab(Interval(c_phi_0, c_phi_0), multigrid_vars_box);

        BoxIterator bit(this_box);
        for (bit.begin(); bit.ok(); ++bit)
        {
            IntVect iv = bit();
            RealVect loc;
            get_loc(loc, iv, a_dx, a_params);
            // m(K, phi) = 2/3 K^2 - 16 pi G rho
            Real m =
                get_m(multigrid_vars_box(iv, c_phi_0), a_params, constant_K);
            Real grad_phi_sq = get_grad_phi_sq(iv, phi_fab, a_dx);

            // Also \bar  A_ij \bar A^ij
            Real A2 = 0.0;
            A2 = pow(multigrid_vars_box(iv, c_A11_0), 2.0) +
                 pow(multigrid_vars_box(iv, c_A22_0), 2.0) +
                 pow(multigrid_vars_box(iv, c_A33_0), 2.0) +
                 2 * pow(multigrid_vars_box(iv, c_A12_0), 2.0) +
                 2 * pow(multigrid_vars_box(iv, c_A13_0), 2.0) +
                 2 * pow(multigrid_vars_box(iv, c_A23_0), 2.0);

            Real psi_bl = get_psi_brill_lindquist(loc, a_params);
            Real psi_0 = multigrid_vars_box(iv, c_psi_reg) + psi_bl;
            aCoef_box(iv, 0) = -0.625 * m * pow(psi_0, 4.0) -
                               0.875 * A2 * pow(psi_0, -8.0) +
                               M_PI * a_params.G_Newton * grad_phi_sq;
        }
    }
}

// The coefficient of the Laplacian operator, for now set to constant 1
// Note that beta = -1 so this sets the sign
// the rhs source of the Poisson eqn
void set_b_coef(LevelData<FArrayBox> &a_bCoef,
                const PoissonParameters &a_params, const RealVect &a_dx)
{

    CH_assert(a_bCoef.nComp() == 1);
    int comp_number = 0;

    for (DataIterator dit = a_bCoef.dataIterator(); dit.ok(); ++dit)
    {
        FArrayBox &bCoef_box = a_bCoef[dit()];
        bCoef_box.setVal(1.0, comp_number);
    }
}

// used to set output data for all ADM Vars for GRChombo restart
void set_output_data(LevelData<FArrayBox> &a_grchombo_vars,
                     LevelData<FArrayBox> &a_multigrid_vars,
                     const PoissonParameters &a_params, const RealVect &a_dx,
                     const Real &constant_K)
{

    CH_assert(a_grchombo_vars.nComp() == NUM_GRCHOMBO_VARS);
    CH_assert(a_multigrid_vars.nComp() == NUM_MULTIGRID_VARS);

    DataIterator dit = a_grchombo_vars.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
        FArrayBox &grchombo_vars_box = a_grchombo_vars[dit()];
        FArrayBox &multigrid_vars_box = a_multigrid_vars[dit()];

        // first set everything to zero
        for (int comp = 0; comp < NUM_GRCHOMBO_VARS; comp++)
        {
            grchombo_vars_box.setVal(0.0, comp);
        }

        // now set non zero terms - const across whole box
        // Conformally flat, and lapse = 1
        grchombo_vars_box.setVal(1.0, c_h11);
        grchombo_vars_box.setVal(1.0, c_h22);
        grchombo_vars_box.setVal(1.0, c_h33);
        grchombo_vars_box.setVal(1.0, c_lapse);

        // constant K
        grchombo_vars_box.setVal(constant_K, c_K);

        // now non constant terms by location
        Box this_box = grchombo_vars_box.box();
        BoxIterator bit(this_box);
        for (bit.begin(); bit.ok(); ++bit)
        {
            IntVect iv = bit();
            RealVect loc;
            get_loc(loc, iv, a_dx, a_params);

            // GRChombo conformal factor chi = psi^-4
            Real psi_bl = get_psi_brill_lindquist(loc, a_params);
            Real chi = pow(multigrid_vars_box(iv, c_psi_reg) + psi_bl, -4.0);
            grchombo_vars_box(iv, c_chi) = chi;
            Real factor = pow(chi, 1.5);

            // Copy phi and Aij across - note this is now \tilde Aij not \bar
            // Aij
            grchombo_vars_box(iv, c_phi) = multigrid_vars_box(iv, c_phi_0);
            grchombo_vars_box(iv, c_A11) =
                multigrid_vars_box(iv, c_A11_0) * factor;
            grchombo_vars_box(iv, c_A12) =
                multigrid_vars_box(iv, c_A12_0) * factor;
            grchombo_vars_box(iv, c_A13) =
                multigrid_vars_box(iv, c_A13_0) * factor;
            grchombo_vars_box(iv, c_A22) =
                multigrid_vars_box(iv, c_A22_0) * factor;
            grchombo_vars_box(iv, c_A23) =
                multigrid_vars_box(iv, c_A23_0) * factor;
            grchombo_vars_box(iv, c_A33) =
                multigrid_vars_box(iv, c_A33_0) * factor;
        }
    }
}
