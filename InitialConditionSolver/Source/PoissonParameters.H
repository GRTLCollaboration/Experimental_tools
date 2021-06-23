#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _POISSONPARAMETERS_H_
#define _POISSONPARAMETERS_H_

#include <iostream>
using std::cerr;

#include "AMRMultiGrid.H"
#include "DisjointBoxLayout.H"
#include "FArrayBox.H"
#include "ParmParse.H"
// #include "GRParmParse.hpp"
#include "VariableCoeffPoissonOperator.H"
#include "GRChomboBCs.hpp"

/// Class for solver input parameters
/**
 */
class PoissonParameters
{
  public:
    IntVect nCells;
    int maxGridSize;
    int blockFactor;
    int bufferSize;
    Real fillRatio;
    bool forceGRChomboRegrid;
    Real refineThresh;
    int coefficient_average_type;
    int verbosity;
    Vector<int> periodic;
    int domBcType;
    int maxLevel;
    int numLevels;
    Vector<int> refRatio;
    ProblemDomain coarsestDomain;
    Real coarsestDx;
    RealVect domainLength;
    RealVect probLo;
    RealVect probHi;
    int num_ghosts;
    GRChomboBCs::params_t grchombo_boundary_params;
    bool symmetric_boundaries_exist;
    bool nonperiodic_boundaries_exist;
    Real alpha;
    Real beta;
    Real G_Newton;
    Real phi_amplitude;
    Real phi_mu;
    Real bh1_bare_mass;
    Real bh2_bare_mass;
    Real bh1_spin;
    Real bh2_spin;
    RealVect bh1_momentum;
    RealVect bh2_momentum;
    RealVect bh1_offset;
    RealVect bh2_offset;
    string filename;
};

/// get stuff from input file
extern void getPoissonParameters(PoissonParameters &a_params);

#endif
