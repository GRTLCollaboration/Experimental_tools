/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef GRCHOMBOBCS_HPP_
#define GRCHOMBOBCS_HPP_

#include "BoxIterator.H"
#include "Copier.H"
#include "FourthOrderInterpStencil.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "RealVect.H"
#include "GRChomboUserVariables.hpp"
#include "MultigridUserVariables.hpp"

/// Class that enables filling of the boundary ghost cells which GRChombo
/// requires for non-periodic boundary conditions
class GRChomboBCs
{
  public:
    /// enum for possible boundary states
    /// Note Static and Sommerfeld do the same thing but just keeping it like
    /// this for constitency with GRChombo
    enum
    {
        STATIC_BC,
        SOMMERFELD_BC,
        REFLECTIVE_BC
    };

    /// enum for possible parity states for Reflective case
    enum
    {
        EVEN,
        ODD_X,
        ODD_Y,
        ODD_Z,
        ODD_XY,
        ODD_YZ,
        ODD_XZ
    };

    /// Structure containing the boundary condition params
    struct params_t
    {
        std::array<int, CH_SPACEDIM> hi_boundary;
        std::array<int, CH_SPACEDIM> lo_boundary;
        std::array<bool, CH_SPACEDIM> is_periodic;
        std::array<int, NUM_GRCHOMBO_VARS> vars_parity;
        //std::array<double, NUM_VARS> vars_asymptotic_values;
    };

  protected:
    // Member values
    double m_dx;            // The grid spacing
    int m_num_ghosts;       // the number of ghosts (usually 3)
    params_t m_params;      // the boundary params
    ProblemDomain m_domain; // the problem domain (excludes boundary cells)
    Box m_domain_box;       // The box representing the domain
    bool is_defined; // whether the GRChomboBCs class members are defined

  public:
    /// Default constructor - need to call define afterwards
    GRChomboBCs() { is_defined = false; }

    /// define function sets members and is_defined set to true
    void define(double a_dx, params_t a_params, ProblemDomain a_domain,
                int a_num_ghosts);


    /// write out boundary params (used during setup for debugging)
    static void write_boundary_conditions(params_t a_params);

    /// The function which returns the parity of each of the vars in
    /// GRChomboUserVariables.hpp The parity should be defined in the params
    /// file, and will be output to the pout files for checking at start/restart
    /// of simulation (It is only required for reflective boundary conditions.)
    int get_vars_parity(int a_comp, int a_dir) const;

    /// static version used for initial output of boundary values
    static int get_vars_parity(int a_comp, int a_dir, params_t a_params);

    /*
    /// Fill the boundary ghosts appropriately based on the params set
    void fill_boundary_ghosts(const Side::LoHiSide a_side,
                              LevelData<FArrayBox> &a_soln);

    /// Fill the boundary values appropriately based on the params set
    /// in the direction dir
    void fill_boundary_ghosts_dir(const Side::LoHiSide a_side,
                                  LevelData<FArrayBox> &a_soln,
                                  const int dir);

    /// Copy the boundary values from src to dest
    /// NB assumes same box layout of input and output data
    void copy_boundary_cells(const Side::LoHiSide a_side,
                             const LevelData<FArrayBox> &a_src,
                             LevelData<FArrayBox> &a_dest);
    */

    /// enforce symmetric boundary conditions, e.g. after interpolation
    void enforce_symmetric_boundaries(const Side::LoHiSide a_side,
                                      LevelData<FArrayBox> &a_state);

    /// Fill the fine boundary values in a_state
    /// Required for interpolating onto finer levels at boundaries
    void interp_boundaries(LevelData<FArrayBox> &a_fine_state,
                           LevelData<FArrayBox> &a_coarse_state,
                           const Side::LoHiSide a_side);

    /// Get the boundary condition for a_dir and a_side
    int get_boundary_condition(const Side::LoHiSide a_side, const int a_dir);

    /// get the boundary box to fill if we are at a boundary
    Box get_boundary_box(const Side::LoHiSide a_side, const int a_dir,
                         const IntVect &offset_lo, const IntVect &offset_hi,
                         Box &this_ghostless_box, int shrink_for_coarse = 0);

    void fill_reflective_cell(FArrayBox &soln_box, const IntVect iv,
                              const Side::LoHiSide a_side, const int dir) const;

    inline void get_box_offsets(IntVect &a_out_offset_lo,
                                IntVect &a_out_offset_hi,
                                const Box &a_box)
    {
        a_out_offset_lo = -a_box.smallEnd() + m_domain_box.smallEnd();
        a_out_offset_hi = +a_box.bigEnd() - m_domain_box.bigEnd();
    }

    inline void remove_outer_ghost_cells(Box &a_box)
    {
        a_box &= m_domain_box;
    }

  private:
    /// write out reflective conditions
    static void write_reflective_conditions(int idir, params_t a_params);

    /*
    void fill_sommerfeld_cell(FArrayBox &soln_box,
                              const IntVect iv) const;
    */
};

#endif /* GRCHOMBOBCS_HPP_ */
