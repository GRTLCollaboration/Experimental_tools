/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// This program can be used to compare two HDF5 files written with GRChombo.
// It is the strongest possible check whether any changes to the program affect
// the output many timesteps later.
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

bool check_level_data(LevelData<FArrayBox> &lev_dat_1,
                      LevelData<FArrayBox> &lev_dat_2)
{
    CH_assert(lev_dat_1.disjointBoxLayout() == lev_dat_2.disjointBoxLayout());

    bool disagree = 0;
    int abort_counter = 0;

    DataIterator dit = lev_dat_1.disjointBoxLayout().dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
        const FArrayBox &fab_1 = lev_dat_1[dit];
        const FArrayBox &fab_2 = lev_dat_2[dit];
        // BoxIterator bit(fab_1.box());
        BoxIterator bit(lev_dat_1.disjointBoxLayout()[dit]);
        for (bit.begin(); bit.ok(); ++bit)
        {
            IntVect iv = bit();
            for (int icomp = 0; icomp < lev_dat_1.nComp(); ++icomp)
            {
                Real difference = fab_1(iv, icomp) - fab_2(iv, icomp);
                if (difference != 0)
                {
                    abort_counter++;
                    pout() << "Values at position " << iv << " in component "
                           << icomp << " disagree by " << difference << endl;
                    pout() << "val 1: " << fab_1(iv, icomp);
                    pout() << " val 2: " << fab_2(iv, icomp) << endl;
                    disagree = true;
                    if (abort_counter > 50)
                    {
                        pout() << "Too many problems on this level. Giving up "
                                  "and going to the next... "
                               << endl;
                        pout() << lev_dat_1.ghostVect() << endl;
                        return disagree;
                    }
                }
            }
        }
    }
    return disagree;
}

int main(int argc, char *argv[])
{
#ifdef CH_MPI
    // Start MPI
    MPI_Init(&argc, &argv);
    int num_procs;
    MPI_Comm_size(Chombo_MPI::comm, &num_procs);
    if (num_procs > 1)
        MayDay::Error("This little tool is not suitable for more than 1 rank "
                      "(although it wouldn't be hard to change this...)");
#ifdef CH_AIX
    H5dont_atexit();
#endif
#endif
    char *in_file = argv[1];
    ParmParse pp(argc - 2, argv + 2, NULL, in_file);

    std::string file_1;
    pp.get("file_1", file_1);
    HDF5Handle handle_1(file_1, HDF5Handle::OPEN_RDONLY);

    std::string file_2;
    pp.get("file_2", file_2);
    HDF5Handle handle_2(file_2, HDF5Handle::OPEN_RDONLY);

    HDF5HeaderData header_1;
    header_1.readFromFile(handle_1);

    HDF5HeaderData header_2;
    header_2.readFromFile(handle_2);

    // read max level
    if (header_1.m_int.find("max_level") == header_1.m_int.end())
        MayDay::Error("File 1 does not contain max_level");
    int max_level_1 = header_1.m_int["max_level"];

    if (header_2.m_int.find("max_level") == header_2.m_int.end())
        MayDay::Error("File 1 does not contain max_level");
    // int max_level_2 = header_2.m_int["max_level"];
    // CH_assert(max_level_1 == max_level_2);

    bool disagree = false;

    int max_comparison_level = max_level_1;
    pp.query("max_comparison_level", max_comparison_level);

    for (int ilev = 0; ilev <= max_comparison_level; ++ilev)
    {
        // Read file 1
        LevelData<FArrayBox> level_data_1;
        Real dx_1, dt_1, time_1;
        Box box_1;
        int ref_ratio_1;
        readLevel(handle_1, ilev, level_data_1, dx_1, dt_1, time_1, box_1,
                  ref_ratio_1, level_data_1.interval());

        // Read file_2
        LevelData<FArrayBox> level_data_2;
        Real dx_2, dt_2, time_2;
        Box box_2;
        int ref_ratio_2;
        readLevel(handle_2, ilev, level_data_2, dx_2, dt_2, time_2, box_2,
                  ref_ratio_2, level_data_2.interval());

        CH_assert(dx_1 == dx_2);
        CH_assert(dt_1 == dt_2);
        CH_assert(time_1 == time_2);
        CH_assert(box_1 == box_2);
        CH_assert(ref_ratio_1 == ref_ratio_2);
        CH_assert(level_data_1.interval() == level_data_2.interval());
        CH_assert(level_data_1.ghostVect() == level_data_2.ghostVect());

        // level_data_1 and level_data_2 might have different box layouts so
        // need to do a copy
        LevelData<FArrayBox> level_data_2_copy(level_data_1.disjointBoxLayout(),
                                               level_data_1.nComp(),
                                               level_data_1.ghostVect());
        level_data_2.copyTo(level_data_2_copy);
        bool level_disagree = check_level_data(level_data_1, level_data_2_copy);
        cout << "Does level " << ilev << " agree? " << !level_disagree << endl;
        disagree |= level_disagree;
    }

    handle_2.close();
    handle_1.close();

#ifdef CH_MPI
    // Exit MPI
    // dumpmemoryatexit();
    MPI_Finalize();
#endif

    return disagree;
}
