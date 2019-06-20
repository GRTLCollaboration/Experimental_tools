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
#include "LoHiSide.H"
#include "PoissonParameters.H"
#include "SetBinaryBH.H"
#include "SetLevelDataF_F.H"
#include "VariableCoeffPoissonOperatorFactory.H"
#include "computeNorm.H"
#include "parstream.H"
#include <cmath>

// Set various LevelData functions across the grid

// This takes an IntVect and writes the physical coordinates to a RealVect
void get_loc(RealVect &a_out_loc, const IntVect &a_iv,
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
                            LevelData<FArrayBox> &a_dpsi,
                            GRChomboBCs &a_grchombo_boundaries,
                            const RealVect &a_dx,
                            const PoissonParameters &a_params)
{

    CH_assert(a_multigrid_vars.nComp() == NUM_MULTIGRID_VARS);

    DataIterator dit = a_multigrid_vars.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
        FArrayBox &multigrid_vars_box = a_multigrid_vars[dit()];
        FArrayBox &dpsi_box = a_dpsi[dit()];
        Box this_box = multigrid_vars_box.box();
        BoxIterator bit(this_box);
        for (bit.begin(); bit.ok(); ++bit)
        {
            set_initial_multigrid_cell(multigrid_vars_box, dpsi_box,
                bit(), a_dx, a_params);
        }

        // now fill boundary ghost cells if using nonperiodic boundaries in
        // GRChombo. Note that these cells are unused in the
        IntVect offset_lo, offset_hi;
        a_grchombo_boundaries.get_box_offsets(offset_lo, offset_hi, this_box);

        // reduce box to the intersection of the box and the
        // problem domain ie remove all outer ghost cells
        a_grchombo_boundaries.remove_outer_ghost_cells(this_box);

        for (int idir = 0; idir < SpaceDim; ++idir)
        {
            if (!a_params.periodic[idir])
            {
                for (SideIterator sit; sit.ok(); ++sit)
                {
                    Box boundary_box = a_grchombo_boundaries.get_boundary_box(
                        sit(), idir, offset_lo, offset_hi, this_box);

                    // now we have the appropriate box, fill it!
                    BoxIterator bbit(boundary_box);
                    for (bbit.begin(); bbit.ok(); ++bbit)
                    {
                        set_initial_multigrid_cell(multigrid_vars_box, dpsi_box,
                            bbit(), a_dx, a_params);
                    } // end loop through boundary box
                } // end loop over sides
            } // end if (periodic[idir])
        } // end loop over directions
    }
} // end set_initial_conditions

void set_initial_multigrid_cell(FArrayBox &a_multigrid_vars_box,
                                FArrayBox &a_dpsi_box,
                                const IntVect &a_iv,
                                const RealVect &a_dx,
                                const PoissonParameters &a_params)
{
    RealVect loc;
    get_loc(loc, a_iv, a_dx, a_params);

    // set psi to 1.0 and zero dpsi
    // note that we don't include the singular part of psi
    // for the BHs - this is added at the output data stage
    // and when we calculate psi_0 in the rhs etc
    // as it already satisfies Laplacian(psi) = 0
    a_multigrid_vars_box(a_iv, c_psi) = 1.0;
    a_dpsi_box(a_iv, 0) = 0.0;

    // set phi according to user defined function
    a_multigrid_vars_box(a_iv, c_phi_0) = my_phi_function(
        loc, a_params.phi_amplitude,
        a_params.phi_wavelength, a_params.domainLength);

    // set Aij for spin and momentum according to BH params
    set_binary_bh_Aij(a_multigrid_vars_box, a_iv, loc,
                      a_params);
} //end set set_initial_multigrid_cell

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

        // calculate the laplacian of psi across the box
        FArrayBox laplacian_of_psi(this_box, 1);
        FORT_GETLAPLACIANPSIF(CHF_FRA1(laplacian_of_psi, 0),
                              CHF_CONST_FRA1(multigrid_vars_box, c_psi),
                              CHF_CONST_REAL(a_dx[0]), CHF_BOX(this_box));

        // calculate the rho contribution from gradients of phi
        FArrayBox rho_gradient(this_box, 1);
        FORT_GETRHOGRADPHIF(CHF_FRA1(rho_gradient, 0),
                            CHF_CONST_FRA1(multigrid_vars_box, c_phi_0),
                            CHF_CONST_REAL(a_dx[0]), CHF_BOX(this_box));

        BoxIterator bit(this_box);
        for (bit.begin(); bit.ok(); ++bit)
        {
            IntVect iv = bit();
            RealVect loc;
            get_loc(loc, iv, a_dx, a_params);

            // rhs = m/8 psi_0^5 - 2 pi rho_grad psi_0  - laplacian(psi_0)
            Real m = 0;
            set_m_value(m, multigrid_vars_box(iv, c_phi_0), a_params,
                        constant_K);

            // Also \bar  A_ij \bar A^ij
            Real A2 = 0.0;
            A2 = pow(multigrid_vars_box(iv, c_A11_0), 2.0) +
                 pow(multigrid_vars_box(iv, c_A22_0), 2.0) +
                 pow(multigrid_vars_box(iv, c_A33_0), 2.0) +
                 2 * pow(multigrid_vars_box(iv, c_A12_0), 2.0) +
                 2 * pow(multigrid_vars_box(iv, c_A13_0), 2.0) +
                 2 * pow(multigrid_vars_box(iv, c_A23_0), 2.0);

            Real psi_bh = set_binary_bh_psi(loc, a_params);
            Real psi_0 = multigrid_vars_box(iv, c_psi) + psi_bh;

            rhs_box(iv, 0) =
                0.125 * m * pow(psi_0, 5.0) - 0.125 * A2 * pow(psi_0, -7.0) -
                2.0 * M_PI * a_params.G_Newton * rho_gradient(iv, 0) * psi_0 -
                laplacian_of_psi(iv, 0);
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

        // calculate the laplacian of psi across the box
        FArrayBox laplacian_of_psi(this_box, 1);
        FORT_GETLAPLACIANPSIF(CHF_FRA1(laplacian_of_psi, 0),
                              CHF_CONST_FRA1(multigrid_vars_box, c_psi),
                              CHF_CONST_REAL(a_dx[0]), CHF_BOX(this_box));

        // calculate the rho contribution from gradients of phi
        FArrayBox rho_gradient(this_box, 1);
        FORT_GETRHOGRADPHIF(CHF_FRA1(rho_gradient, 0),
                            CHF_CONST_FRA1(multigrid_vars_box, c_phi_0),
                            CHF_CONST_REAL(a_dx[0]), CHF_BOX(this_box));

        BoxIterator bit(this_box);
        for (bit.begin(); bit.ok(); ++bit)
        {
            IntVect iv = bit();
            RealVect loc;
            get_loc(loc, iv, a_dx, a_params);

            // integrand = -1.5*m + 1.5 * \bar A_ij \bar A^ij psi_0^-12 +
            // 24 pi rho_grad psi_0^-4  + 12*laplacian(psi_0)*psi^-5
            Real m = 0;
            set_m_value(m, multigrid_vars_box(iv, c_phi_0), a_params, 0.0);

            // Also \bar  A_ij \bar A^ij
            Real A2 = 0.0;
            A2 = pow(multigrid_vars_box(iv, c_A11_0), 2.0) +
                 pow(multigrid_vars_box(iv, c_A22_0), 2.0) +
                 pow(multigrid_vars_box(iv, c_A33_0), 2.0) +
                 2 * pow(multigrid_vars_box(iv, c_A12_0), 2.0) +
                 2 * pow(multigrid_vars_box(iv, c_A13_0), 2.0) +
                 2 * pow(multigrid_vars_box(iv, c_A23_0), 2.0);

            Real psi_bh = set_binary_bh_psi(loc, a_params);
            Real psi_0 = multigrid_vars_box(iv, c_psi) + psi_bh;

            integrand_box(iv, 0) =
                -1.5 * m + 1.5 * A2 * pow(psi_0, -12.0) +
                24.0 * M_PI * a_params.G_Newton * rho_gradient(iv, 0) *
                    pow(psi_0, -4.0) +
                12.0 * laplacian_of_psi(iv, 0) * pow(psi_0, -5.0);
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

        // calculate the rho contribution from gradients of phi
        FArrayBox rho_gradient(this_box, 1);
        FORT_GETRHOGRADPHIF(CHF_FRA1(rho_gradient, 0),
                            CHF_CONST_FRA1(multigrid_vars_box, c_phi_0),
                            CHF_CONST_REAL(a_dx[0]), CHF_BOX(this_box));

        BoxIterator bit(this_box);
        for (bit.begin(); bit.ok(); ++bit)
        {
            IntVect iv = bit();
            RealVect loc;
            get_loc(loc, iv, a_dx, a_params);

            // calculate contributions
            Real m = 0;
            set_m_value(m, multigrid_vars_box(iv, c_phi_0), a_params, 0.0);

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
            Real psi_bh = set_binary_bh_psi(loc, a_params);
            Real psi_0 = multigrid_vars_box(iv, c_psi) + psi_bh;
            condition_box(iv, 0) = 1.5 * abs(m) + 1.5 * A2 * pow(psi_0, -7.0) +
                                   24.0 * M_PI * a_params.G_Newton *
                                       abs(rho_gradient(iv, 0)) *
                                       pow(psi_0, 1.0) +
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
            multigrid_vars_box(iv, c_psi) += dpsi_box(iv, 0);
        }
    }
}

// m(K, rho) = 2/3K^2 - 16piG rho
void set_m_value(Real &m, const Real &phi_here,
                 const PoissonParameters &a_params, const Real constant_K)
{

    // KC TODO:
    // For now rho is just the gradient term which is kept separate
    // ... may want to add V(phi) and phidot/Pi here later though
    Real Pi_field = 0.0;
    Real V_of_phi = 0.0;
    Real rho = 0.5 * Pi_field * Pi_field + V_of_phi;

    m = (2.0 / 3.0) * (constant_K * constant_K) -
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

        // calculate the rho contribution from gradients of phi
        FArrayBox rho_gradient(this_box, 1);
        FORT_GETRHOGRADPHIF(CHF_FRA1(rho_gradient, 0),
                            CHF_CONST_FRA1(multigrid_vars_box, c_phi_0),
                            CHF_CONST_REAL(a_dx[0]), CHF_BOX(this_box));

        BoxIterator bit(this_box);
        for (bit.begin(); bit.ok(); ++bit)
        {
            IntVect iv = bit();
            RealVect loc;
            get_loc(loc, iv, a_dx, a_params);
            // m(K, phi) = 2/3 K^2 - 16 pi G rho
            Real m;
            set_m_value(m, multigrid_vars_box(iv, c_phi_0), a_params,
                        constant_K);

            // Also \bar  A_ij \bar A^ij
            Real A2 = 0.0;
            A2 = pow(multigrid_vars_box(iv, c_A11_0), 2.0) +
                 pow(multigrid_vars_box(iv, c_A22_0), 2.0) +
                 pow(multigrid_vars_box(iv, c_A33_0), 2.0) +
                 2 * pow(multigrid_vars_box(iv, c_A12_0), 2.0) +
                 2 * pow(multigrid_vars_box(iv, c_A13_0), 2.0) +
                 2 * pow(multigrid_vars_box(iv, c_A23_0), 2.0);

            Real psi_bh = set_binary_bh_psi(loc, a_params);
            Real psi_0 = multigrid_vars_box(iv, c_psi) + psi_bh;
            aCoef_box(iv, 0) =
                -0.625 * m * pow(psi_0, 4.0) - A2 * pow(psi_0, -8.0) +
                2.0 * M_PI * a_params.G_Newton * rho_gradient(iv, 0);
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
                     GRChomboBCs &a_grchombo_boundaries,
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
            set_non_const_output_cell(multigrid_vars_box,
                grchombo_vars_box, bit(), a_dx, a_params);
        }

        // finally non-constant boundary ghosts
        IntVect offset_lo, offset_hi;
        a_grchombo_boundaries.get_box_offsets(offset_lo, offset_hi, this_box);

        // reduce box to the intersection of the box and the
        // problem domain ie remove all outer ghost cells
        a_grchombo_boundaries.remove_outer_ghost_cells(this_box);

        // get the boundary box (may be Empty)
        for (int idir = 0; idir < SpaceDim; ++idir)
        {
            if (!a_params.periodic[idir])
            {
                for (SideIterator sit; sit.ok(); ++sit)
                {
                    Box boundary_box = a_grchombo_boundaries.get_boundary_box(
                        sit(), idir, offset_lo, offset_hi, this_box);

                    // now we have the appropriate box, fill it!
                    BoxIterator bbit(boundary_box);
                    for (bbit.begin(); bbit.ok(); ++bbit)
                    {
                        set_non_const_output_cell(multigrid_vars_box,
                            grchombo_vars_box, bbit(), a_dx, a_params);
                    } // end loop through boundary box
                } // end loop over sides
            } // end if (periodic[idir])
        } // end loop over directions
    } // end loop over boxes
} // end set_output_data

void set_non_const_output_cell(const FArrayBox &a_multigrid_vars_box,
                               FArrayBox &a_grchombo_vars_box,
                               const IntVect &a_iv,
                               const RealVect &a_dx,
                               const PoissonParameters &a_params)
{
    RealVect loc;
    get_loc(loc, a_iv, a_dx, a_params);

    // GRChombo conformal factor chi = psi^-4
    Real psi_bh = set_binary_bh_psi(loc, a_params);
    Real chi = pow(a_multigrid_vars_box(a_iv, c_psi) + psi_bh, -4.0);
    a_grchombo_vars_box(a_iv, c_chi) = chi;
    Real factor = pow(chi, 1.5);

    // Copy phi and Aij across - note this is now \tilde Aij not \bar
    // Aij
    a_grchombo_vars_box(a_iv, c_phi) = a_multigrid_vars_box(a_iv, c_phi_0);
    a_grchombo_vars_box(a_iv, c_A11) =
        a_multigrid_vars_box(a_iv, c_A11_0) * factor;
    a_grchombo_vars_box(a_iv, c_A12) =
        a_multigrid_vars_box(a_iv, c_A12_0) * factor;
    a_grchombo_vars_box(a_iv, c_A13) =
        a_multigrid_vars_box(a_iv, c_A13_0) * factor;
    a_grchombo_vars_box(a_iv, c_A22) =
        a_multigrid_vars_box(a_iv, c_A22_0) * factor;
    a_grchombo_vars_box(a_iv, c_A23) =
        a_multigrid_vars_box(a_iv, c_A23_0) * factor;
    a_grchombo_vars_box(a_iv, c_A33) =
        a_multigrid_vars_box(a_iv, c_A33_0) * factor;
}
