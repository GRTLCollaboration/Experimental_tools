#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "PoissonParameters.H"
#include "AMRIO.H"
#include "BCFunc.H"
#include "BRMeshRefine.H"
#include "BiCGStabSolver.H"
#include "BoxIterator.H"
#include "CONSTANTS.H"
#include "CoarseAverage.H"
#include "LoadBalance.H"
#include "computeNorm.H"
#include "parstream.H"
#include <cmath>

// function to read in the key params for solver
void getPoissonParameters(PoissonParameters &a_params)
{
    ParmParse pp;

    // problem specific params
    pp.get("alpha", a_params.alpha);
    pp.get("beta", a_params.beta);
    // print out the overall coeffs just to be sure we have selected them
    // correctly
    pout() << "alpha, beta = " << a_params.alpha << ", " << a_params.beta
           << endl;

    // Initial conditions for the scalar field
    pp.get("G_Newton", a_params.G_Newton);
    pp.get("phi_amplitude", a_params.phi_amplitude);
    pp.get("phi_wavelength", a_params.phi_wavelength);

    if (abs(a_params.phi_amplitude) > 0.0)
    {
        pout() << "Spacetime contains scalar field of amplitude "
               << a_params.phi_amplitude << endl;
    }

    // Initial conditions for the black holes
    pp.get("bh1_bare_mass", a_params.bh1_bare_mass);
    pp.get("bh2_bare_mass", a_params.bh2_bare_mass);
    pp.get("bh1_spin", a_params.bh1_spin);
    pp.get("bh2_spin", a_params.bh2_spin);
    pp.get("bh1_offset", a_params.bh1_offset);
    pp.get("bh2_offset", a_params.bh2_offset);
    pp.get("bh1_momentum", a_params.bh1_momentum);
    pp.get("bh2_momentum", a_params.bh2_momentum);

    if (abs(a_params.bh1_bare_mass) > 0.0 || abs(a_params.bh2_bare_mass) > 0.0)
    {
        pout() << "Spacetime contains black holes with bare masses "
               << a_params.bh1_bare_mass << " and " << a_params.bh2_bare_mass
               << endl;
    }

    // Set verbosity
    a_params.verbosity = 3;
    pp.query("verbosity", a_params.verbosity);

    // Chombo grid params
    pp.get("max_level", a_params.maxLevel);
    a_params.numLevels = a_params.maxLevel + 1;
    std::vector<int> nCellsArray(SpaceDim);
    pp.getarr("N", nCellsArray, 0, SpaceDim);
    for (int idir = 0; idir < SpaceDim; idir++)
    {
        a_params.nCells[idir] = nCellsArray[idir];
    }

    // Enforce that dx is same in every directions
    // and that ref_ratio = 2 always as these conditions
    // are required in several places in our code
    a_params.refRatio.resize(a_params.numLevels);
    a_params.refRatio.assign(2);
    Real domain_length;
    pp.get("L", domain_length);
    a_params.coarsestDx = domain_length / a_params.nCells[0];
    for (int idir = 0; idir < SpaceDim; idir++)
    {
        a_params.domainLength[idir] =
            a_params.coarsestDx * a_params.nCells[idir];
    }

    // Chombo refinement and load balancing criteria
    pp.get("refine_threshold", a_params.refineThresh);
    pp.get("block_factor", a_params.blockFactor);
    pp.get("max_grid_size", a_params.maxGridSize);
    pp.get("fill_ratio", a_params.fillRatio);
    pp.get("buffer_size", a_params.bufferSize);

    // set average type -
    // set to a bogus default value, so we only break from solver
    // default if it's set to something real
    a_params.coefficient_average_type = -1;
    if (pp.contains("coefficient_average_type"))
    {
        std::string tempString;
        pp.get("coefficient_average_type", tempString);
        if (tempString == "arithmetic")
        {
            a_params.coefficient_average_type = CoarseAverage::arithmetic;
        }
        else if (tempString == "harmonic")
        {
            a_params.coefficient_average_type = CoarseAverage::harmonic;
        }
        else
        {
            MayDay::Error("bad coefficient_average_type in input");
        }
    } // end if an average_type is present in inputs

    // set up coarse domain box
    IntVect lo = IntVect::Zero;
    IntVect hi = a_params.nCells;
    hi -= IntVect::Unit;
    Box crseDomBox(lo, hi);
    a_params.probLo = RealVect::Zero;
    a_params.probHi = RealVect::Zero;
    a_params.probHi += a_params.domainLength;

    // Hardcode num_ghosts to 3 as this is what GRChombo needs
    a_params.num_ghosts = 3;

    // Periodicity - for the moment enforce same in all directions
    ProblemDomain crseDom(crseDomBox);
    int is_periodic;
    pp.get("is_periodic", is_periodic);
    a_params.periodic.resize(SpaceDim);
    a_params.periodic.assign(is_periodic);
    for (int dir = 0; dir < SpaceDim; dir++)
    {
        crseDom.setPeriodic(dir, is_periodic);
    }
    a_params.coarsestDomain = crseDom;

    // Load GRChombo boundary params
    Vector<int> grchombo_hi_boundary(SpaceDim, GRChomboBCs::STATIC_BC);
    Vector<int> grchombo_lo_boundary(SpaceDim, GRChomboBCs::STATIC_BC);
    pp.queryarr("hi_boundary", grchombo_hi_boundary, 0, SpaceDim);
    pp.queryarr("lo_boundary", grchombo_lo_boundary, 0, SpaceDim);

    // set defaults and override below
    Vector<int> vars_boundary_parity(NUM_MULTIGRID_VARS, GRChomboBCs::EVEN);
    pout() << "periodicity = " << is_periodic << endl;
    for(int idir = 0; idir < SpaceDim; ++idir)
    {
        a_params.grchombo_boundary_params.hi_boundary[idir] =
            grchombo_hi_boundary[idir];
        a_params.grchombo_boundary_params.lo_boundary[idir] =
            grchombo_lo_boundary[idir];
        a_params.grchombo_boundary_params.is_periodic[idir] =
            a_params.periodic[idir];

    }
    a_params.nonperiodic_boundaries_exist = false;
    a_params.symmetric_boundaries_exist = false;

    for (int idir = 0; idir < SpaceDim; ++idir)
    {
        if (!a_params.periodic[idir])
        {
            a_params.nonperiodic_boundaries_exist = true;
            if ((grchombo_hi_boundary[idir] ==
                 GRChomboBCs::REFLECTIVE_BC) ||
                (grchombo_lo_boundary[idir] ==
                 GRChomboBCs::REFLECTIVE_BC))
            {
                a_params.symmetric_boundaries_exist = true;
                pp.getarr("vars_parity", vars_boundary_parity, 0,
                          NUM_MULTIGRID_VARS);

            }
        }
    }

    for (int ivar = 0; ivar < NUM_MULTIGRID_VARS; ++ivar)
    {
        a_params.grchombo_boundary_params.vars_parity[ivar]
            = vars_boundary_parity[ivar];
    }

    if (a_params.nonperiodic_boundaries_exist)
    {
        // write out boundary conditions where non periodic - useful for
        // debug
        GRChomboBCs::write_boundary_conditions(
            a_params.grchombo_boundary_params);
    }


}
