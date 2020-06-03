
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

void defineExchangeCopier(const DisjointBoxLayout &a_level_grids,
                          DisjointBoxLayout grown_grids, Copier exchange_copier)
{
    IntVect iv_ghosts = 3 * IntVect::Unit;
    exchange_copier.exchangeDefine(grown_grids, iv_ghosts);
}

void read_level_data(LevelData<FArrayBox> &lev_dat_1,
                     LevelData<FArrayBox> &lev_dat_2)
{
    CH_assert(lev_dat_1.disjointBoxLayout() == lev_dat_2.disjointBoxLayout());

    DataIterator dit = lev_dat_1.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
        FArrayBox &fab_1 = lev_dat_1[dit];
        FArrayBox &fab_2 = lev_dat_2[dit];
        Box b = fab_2.box();
        BoxIterator bit(b);
        for (bit.begin(); bit.ok(); ++bit)
        {
            IntVect iv = bit();
            //            for (int icomp = 0; icomp < lev_dat_1.nComp();
            //            ++icomp)
            //            {
            //                Real difference = fab_1(iv, icomp) - fab_2(iv,
            //                icomp); cout << "hdf5: " << fab_1(iv, c_phi) <<
            //                endl;
            fab_2(iv, c_phi_0) = fab_1(iv, c_phi);
            //                cout << "hdf5: " << fab_1(iv, c_phi) << endl;

            //            }
        }
    }
}

void readHDF5(LevelData<FArrayBox> &a_multigrid_vars,
              const PoissonParameters &a_params, int ilev)
{
#ifdef CH_USE_HDF5

    HDF5Handle handle(a_params.read_from_file, HDF5Handle::OPEN_RDONLY);

    pout() << "Reading data from " << a_params.read_from_file << endl;

    HDF5HeaderData header;
    header.readFromFile(handle);

    // Read file
    const DisjointBoxLayout &level_domain = a_multigrid_vars.disjointBoxLayout();
    Real dx, dt, time;
    Box box;
    int ref_ratio;

    Copier exchange_copier;
    LevelData<FArrayBox> level_data_noghost(level_domain, NUM_GRCHOMBO_VARS, 0 * IntVect::Unit);
    LevelData<FArrayBox> level_data_ghost(level_domain, NUM_GRCHOMBO_VARS, 3 * IntVect::Unit);

    readLevel(handle, ilev, level_data_noghost, dx, dt, time, box,
                  ref_ratio, level_data_noghost.interval(), false);

    level_data_noghost.copyTo(level_data_ghost);
//    defineExchangeCopier(level_domain, level_data_noghost.disjointBoxLayout(), exchange_copier);
//    level_data_ghost.exchange(exchange_copier);

    LevelData<FArrayBox> level_data_readin(a_multigrid_vars.disjointBoxLayout(),
                                           level_data_ghost.nComp(),
                                           level_data_ghost.ghostVect());
    level_data_ghost.copyTo(level_data_readin);

    read_level_data(level_data_readin, a_multigrid_vars);

    handle.close();

#endif
}
#endif
