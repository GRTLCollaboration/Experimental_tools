/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _READHDF5_H_
#define _READHDF5_H_

#include <iostream>

#include "BoxIterator.H"
#include "CH_HDF5.H"
#include "DataIterator.H"
#include "DisjointBoxLayout.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "ParmParse.H"
#include "ProblemDomain.H"
#include <sys/time.h>

#include "PoissonParameters.H"
#include "SetLevelData.H"

void read_level_data(LevelData<FArrayBox> &lev_dat_1,
                     LevelData<FArrayBox> &lev_dat_2)
{
    CH_assert(lev_dat_1.disjointBoxLayout() == lev_dat_2.disjointBoxLayout());

    bool disagree = 0;
    int abort_counter = 0;

    DataIterator dit = lev_dat_1.disjointBoxLayout().dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
        const FArrayBox &fab_1 = lev_dat_1[dit];
        FArrayBox &fab_2 = lev_dat_2[dit];
        // BoxIterator bit(fab_1.box());
        BoxIterator bit(lev_dat_1.disjointBoxLayout()[dit]);
        for (bit.begin(); bit.ok(); ++bit)
        {
            IntVect iv = bit();
            for (int icomp = 0; icomp < lev_dat_1.nComp(); ++icomp)
            {
//                Real difference = fab_1(iv, icomp) - fab_2(iv, icomp);
                fab_2(iv, c_phi_0) = fab_1(iv, c_phi);
            }
        }
    }
}

void readHDF5(LevelData<FArrayBox> &a_multigrid_vars,
              const PoissonParameters &a_params, int ilev)
{
#ifdef CH_USE_HDF5

    HDF5Handle handle_1(a_params.read_from_file, HDF5Handle::OPEN_RDONLY);
    HDF5HeaderData header;
    header.readFromFile(handle_1);

    // read max level
    if (header.m_int.find("max_level") == header.m_int.end())
        MayDay::Error("File 1 does not contain max_level");
    int max_level_1 = header.m_int["max_level"];

    bool disagree = false;

    int max_comparison_level = max_level_1;
//    pp.query("max_comparison_level", max_comparison_level);


//    for (int ilev = 0; ilev <= max_comparison_level; ++ilev)
//    {

        // Read file 1
    LevelData<FArrayBox> level_data_1;
    Real dx_1, dt_1, time_1;
    Box box_1;
    int ref_ratio_1;
    readLevel(handle_1, ilev, level_data_1, dx_1, dt_1, time_1, box_1,
              ref_ratio_1, level_data_1.interval());

        // level_data_1 and level_data_2 might have different box layouts so
        // need to do a copy
    LevelData<FArrayBox> level_data_2_copy(level_data_1.disjointBoxLayout(),
                                           level_data_1.nComp(),
                                           level_data_1.ghostVect());
    a_multigrid_vars.copyTo(level_data_2_copy);
    read_level_data(level_data_1, level_data_2_copy);
//    }

    handle_1.close();

#endif
}
#endif