/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "GRChomboBCs.hpp"
#include "ProblemDomain.H"
#include <array>
#include <map>
#include <string>

/// define function sets members and is_defined set to true
void GRChomboBCs::define(double a_dx, params_t a_params, ProblemDomain a_domain,
                         int a_num_ghosts)
{
    m_dx = a_dx;
    m_params = a_params;
    m_domain = a_domain;
    m_domain_box = a_domain.domainBox();
    m_num_ghosts = a_num_ghosts;
    is_defined = true;
}

void GRChomboBCs::write_reflective_conditions(int idir,
                                              params_t a_params)
{
    pout() << "The variables that are parity odd in this direction are : "
           << endl;
    for (int icomp = 0; icomp < NUM_MULTIGRID_VARS; icomp++)
    {
        int parity = get_vars_parity(icomp, idir, a_params);
        if (parity == -1)
        {
            pout() << GRChomboUserVariables::variable_names[icomp] << "    ";
        }
    }
}

/// write out boundary params (used during setup for debugging)
void GRChomboBCs::write_boundary_conditions(params_t a_params)
{
    pout() << "You are using non periodic boundary conditions." << endl;
    pout() << "The GRChombo boundary params chosen are:  " << endl;
    pout() << "---------------------------------" << endl;

    std::map<int, std::string> bc_names = {{STATIC_BC, "Static"},
                                           {SOMMERFELD_BC, "Sommerfeld"},
                                           {REFLECTIVE_BC, "Reflective"}};
    for(int idir = 0; idir < SpaceDim; ++idir)
    {
        if (!a_params.is_periodic[idir])
        {
            pout() << "- " << bc_names[a_params.hi_boundary[idir]]
                   << " boundaries in direction high " << idir << endl;
            // high directions
            if (a_params.hi_boundary[idir] == REFLECTIVE_BC)
            {
                write_reflective_conditions(idir, a_params);
            }
            pout() << endl;

            // low directions
            pout() << "- " << bc_names[a_params.lo_boundary[idir]]
                   << " boundaries in direction low " << idir << endl;
            if (a_params.lo_boundary[idir] == REFLECTIVE_BC)
            {
                write_reflective_conditions(idir, a_params);
            }
            pout() << endl;
        }
    }
    pout() << "---------------------------------" << endl;
}

/// The function which returns the parity of each of the vars in
/// UserVariables.hpp The parity should be defined in the params file, and
/// will be output to the pout files for checking at start/restart of
/// simulation (It is only required for reflective boundary conditions.)
int GRChomboBCs::get_vars_parity(int a_comp, int a_dir) const
{
    int vars_parity = get_vars_parity(a_comp, a_dir, m_params);

    return vars_parity;
}

/// static version used for initial output of boundary values
int GRChomboBCs::get_vars_parity(int a_comp, int a_dir,
                                 params_t a_params)
{
    int vars_parity = 1;
    if ((a_dir == 0) && (a_params.vars_parity[a_comp] == ODD_X ||
                         a_params.vars_parity[a_comp] == ODD_XY ||
                         a_params.vars_parity[a_comp] == ODD_XZ))
    {
        vars_parity = -1;
    }
    else if ((a_dir == 1) && (a_params.vars_parity[a_comp] == ODD_Y ||
                              a_params.vars_parity[a_comp] == ODD_XY ||
                              a_params.vars_parity[a_comp] == ODD_YZ))
    {
        vars_parity = -1;
    }
    else if ((a_dir == 2) && (a_params.vars_parity[a_comp] == ODD_Z ||
                              a_params.vars_parity[a_comp] == ODD_XZ ||
                              a_params.vars_parity[a_comp] == ODD_YZ))
    {
        vars_parity = -1;
    }
    return vars_parity;
}

/*
/// Fill the boundary values appropriately based on the params set
void GRChomboBCs::fill_boundary_ghosts(const Side::LoHiSide a_side,
                                       LevelData<FArrayBox> &a_soln);
{
    CH_assert(is_defined);
    CH_TIME("GRChomboBCs::fill_boundary_ghosts");

    // cycle through the directions, filling the cells
    for(int idir; idir < SpaceDim; ++idir)
    {
        // only do something if this direction is not periodic
        if (!m_params.is_periodic[idir])
        {
            fill_boundary_ghosts_dir(a_side, a_soln, idir);
        }
    }
}

void GRChomboBCs::fill_sommerfeld_cell(FArrayBox &soln_box,
                                       const IntVect iv) const
{
    // assumes an asymptotic value + radial waves and permits them
    // to exit grid with minimal reflections
    // get real position on the grid
    RealVect loc(iv + 0.5 * RealVect::Unit);
    loc *= m_dx;
    loc -= m_center;
    double radius_squared = 0.0;
    for(int idir; idir < SpaceDim; ++idir)
    {
        radius_squared += loc[idir] * loc[idir];
    }
    double radius = sqrt(radius_squared);
    IntVect lo_local_offset = iv - soln_box.smallEnd();
    IntVect hi_local_offset = soln_box.bigEnd() - iv;

    // Apply Sommerfeld BCs to each variable
    for (int icomp = 0; icomp < NUM_VARS; icomp++)
    {
        soln_box(iv, icomp) = 0.0;
        for(int idir2; idir2 < SpaceDim; ++idir2)
        {
            IntVect iv_offset1 = iv;
            IntVect iv_offset2 = iv;
            double d1;
            // bit of work to get the right stencils for near
            // the edges of the domain, only using second order
            // stencils for now
            if (lo_local_offset[idir2] < 1)
            {
                // near lo end
                iv_offset1[idir2] += +1;
                iv_offset2[idir2] += +2;
                d1 = 1.0 / m_dx *
                     (-1.5 * soln_box(iv, icomp) +
                      2.0 * soln_box(iv_offset1, icomp) -
                      0.5 * soln_box(iv_offset2, icomp));
            }
            else if (hi_local_offset[idir2] < 1)
            {
                // near hi end
                iv_offset1[idir2] += -1;
                iv_offset2[idir2] += -2;
                d1 = 1.0 / m_dx *
                     (+1.5 * soln_box(iv, icomp) -
                      2.0 * soln_box(iv_offset1, icomp) +
                      0.5 * soln_box(iv_offset2, icomp));
            }
            else
            {
                // normal case
                iv_offset1[idir2] += +1;
                iv_offset2[idir2] += -1;
                d1 =
                    0.5 / m_dx *
                    (soln_box(iv_offset1, icomp) - soln_box(iv_offset2, icomp));
            }

            // for each direction add dphidx * x^i / r
            soln_box(iv, icomp) += -d1 * loc[idir2] / radius;
        }

        // asymptotic values - these need to have been set in
        // the params file
        soln_box(iv, icomp) +=
            (m_params.vars_asymptotic_values[icomp] - soln_box(iv, icomp)) /
            radius;
    }
}
*/

void GRChomboBCs::fill_reflective_cell(FArrayBox &soln_box,
                                       const IntVect iv,
                                       const Side::LoHiSide a_side,
                                       const int dir) const
{
    // assume boundary is a reflection of values within the grid
    // care must be taken with variable parity to maintain correct
    // values on reflection, e.g. x components of vectors are odd
    // parity in the x direction
    IntVect iv_copy = iv;
    /// where to copy the data from - mirror image in domain
    if (a_side == Side::Lo)
    {
        iv_copy[dir] = -iv[dir] - 1;
    }
    else
    {
        iv_copy[dir] = 2 * m_domain_box.bigEnd(dir) - iv[dir] + 1;
    }

    // replace value at iv with value at iv_copy
    for (int icomp = 0; icomp < NUM_MULTIGRID_VARS; icomp++)
    {
        int parity = get_vars_parity(icomp, dir);
        soln_box(iv, icomp) = parity * soln_box(iv_copy, icomp);
    }
}

/*
/// Fill the boundary values appropriately based on the params set
/// in the direction dir
void GRChomboBCs::fill_boundary_ghosts_dir(const Side::LoHiSide a_side,
                                           LevelData<FArrayBox> &a_soln,
                                           const int dir)
{
    // iterate through the boxes, shared amongst threads
    DataIterator dit = a_soln.dataIterator();
    int nbox = dit.size();
#pragma omp parallel for default(shared)
    for (int ibox = 0; ibox < nbox; ++ibox)
    {
        DataIndex dind = dit[ibox];
        FArrayBox &soln_box = a_soln[dind];
        Box this_box = soln_box.box();
        IntVect offset_lo = -this_box.smallEnd() + m_domain_box.smallEnd();
        IntVect offset_hi = +this_box.bigEnd() - m_domain_box.bigEnd();

        // reduce box to the intersection of the box and the
        // problem domain ie remove all outer ghost cells
        this_box &= m_domain_box;
        // get the boundary box (may be Empty) and the condition on it
        int boundary_condition = get_boundary_condition(a_side, dir);
        Box boundary_box =
            get_boundary_box(a_side, dir, offset_lo, offset_hi, this_box);

        // now we have the appropriate box, fill it!
        BoxIterator bit(boundary_box);
        for (bit.begin(); bit.ok(); ++bit)
        {
            IntVect iv = bit();
            switch (boundary_condition)
            {
            // simplest case - boundary values are fixed to the initial ones
            case STATIC_BC:
            {
                for (int icomp = 0; icomp < NUM_VARS; icomp++)
                {
                    soln_box(iv, icomp) = 0.0;
                }
                break;
            }
            case SOMMERFELD_BC:
            {
                fill_sommerfeld_cell(soln_box, iv);
                break;
            }
            case REFLECTIVE_BC:
            {
                fill_reflective_cell(soln_box, iv, a_side, dir);
                break;
            }
            default:
                MayDay::Error(
                    "BoundaryCondition::Supplied boundary not supported.");
            } // end switch
        }     // end iterate over box
    }         // end iterate over boxes
}

/// Copy the boundary values from src to dest
/// NB assumes same box layout of input and output data
void GRChomboBCs::copy_boundary_cells(const Side::LoHiSide a_side,
                                      const LevelData<FArrayBox> &a_src,
                                      LevelData<FArrayBox> &a_dest)
{
    CH_TIME("GRChomboBCs::copy_boundary_cells");

    CH_assert(is_defined);
    if (a_src.boxLayout() == a_dest.boxLayout())
    {
        // cycle through the directions
        for(int idir; idir < SpaceDim; ++idir)
        {
            // only do something if this direction is not periodic
            if (!m_params.is_periodic[idir])
            {
                // iterate through the boxes, shared amongst threads
                DataIterator dit = a_dest.dataIterator();
                int nbox = dit.size();
#pragma omp parallel for default(shared)
                for (int ibox = 0; ibox < nbox; ++ibox)
                {
                    DataIndex dind = dit[ibox];
                    FArrayBox &m_dest_box = a_dest[dind];
                    Box this_box = m_dest_box.box();
                    IntVect offset_lo =
                        -this_box.smallEnd() + m_domain_box.smallEnd();
                    IntVect offset_hi =
                        +this_box.bigEnd() - m_domain_box.bigEnd();

                    // reduce box to the intersection of the box and the
                    // problem domain ie remove all outer ghost cells
                    this_box &= m_domain_box;

                    // get the boundary box (may be Empty)
                    Box boundary_box = get_boundary_box(a_side, idir, offset_lo,
                                                        offset_hi, this_box);

                    BoxIterator bit(boundary_box);
                    for (bit.begin(); bit.ok(); ++bit)
                    {
                        IntVect iv = bit();
                        for (int icomp = 0; icomp < NUM_VARS; icomp++)
                        {
                            m_dest_box(iv, icomp) = a_src[dind](iv, icomp);
                        }
                    } // end iterate over box
                }     // end iterate over boxes
            }         // end if(not periodic)
        }             // end iterate over spacedims
    }                 // end test for same box layout
}
*/

/// enforce symmetric boundary conditions, e.g. after interpolation
void GRChomboBCs::enforce_symmetric_boundaries(
    const Side::LoHiSide a_side, LevelData<FArrayBox> &a_state)
{
    CH_assert(is_defined);
    CH_TIME("GRChomboBCs::enforce_symmetric_boundaries");

    // cycle through the directions
    for(int idir = 0; idir < SpaceDim; ++idir)
    {
        // only do something if this direction is not periodic and symmetric
        if (!m_params.is_periodic[idir])
        {
            int boundary_condition = get_boundary_condition(a_side, idir);

            if (boundary_condition == REFLECTIVE_BC)
            {
                // iterate through the boxes, shared amongst threads
                DataIterator dit = a_state.dataIterator();
                int nbox = dit.size();
                for (int ibox = 0; ibox < nbox; ++ibox)
                {
                    DataIndex dind = dit[ibox];
                    FArrayBox &soln_box = a_state[dind];
                    Box this_box = soln_box.box();
                    IntVect offset_lo = -this_box.smallEnd()
                                        + m_domain_box.smallEnd();
                    IntVect offset_hi = +this_box.bigEnd()
                                        - m_domain_box.bigEnd();

                    // reduce box to the intersection of the box and the
                    // problem domain ie remove all outer ghost cells
                    this_box &= m_domain_box;
                    // get the boundary box (may be Empty)
                    Box boundary_box = get_boundary_box(a_side, idir, offset_lo,
                                                        offset_hi, this_box);

                    // now we have the appropriate box, fill it!
                    BoxIterator bit(boundary_box);
                    for (bit.begin(); bit.ok(); ++bit)
                    {
                        IntVect iv = bit();
                        fill_reflective_cell(soln_box, iv, a_side, idir);
                    }     // end iterate over box
                }         // end iterate over boxes
            }
        }
    }
}

/// Fill the fine boundary values in a_state
/// Required for interpolating onto finer levels at boundaries
void GRChomboBCs::interp_boundaries(LevelData<FArrayBox> &a_fine_state,
                                    LevelData<FArrayBox> &a_coarse_state,
                                    const Side::LoHiSide a_side)
{
    CH_assert(is_defined);
    CH_TIME("GRChomboBCs::interp_boundaries");

    // cycle through the directions
    for(int idir = 0; idir < SpaceDim; ++idir)
    {
        // only do something if this direction is not periodic
        if (!m_params.is_periodic[idir])
        {
            // Ref ratio is always two
            int ref_ratio = 2;

            // create a coarsened fine layout and copy the coarse data onto
            // it
            DisjointBoxLayout coarsened_layout;
            coarsen(coarsened_layout, a_fine_state.disjointBoxLayout(),
                    ref_ratio * IntVect::Unit);
            LevelData<FArrayBox> coarsened_fine;
            coarsened_fine.define(coarsened_layout, NUM_MULTIGRID_VARS,
                                  m_num_ghosts * IntVect::Unit);
            Box coarse_domain_box = coarsen(m_domain_box, ref_ratio);

            // trick the copyTo into thinking the boundary cells are within
            // the domain by growing the domain
            Box grown_domain_box = coarse_domain_box;
            grown_domain_box.grow(m_num_ghosts * IntVect::Unit);
            Copier boundary_copier;
            boundary_copier.ghostDefine(
                a_coarse_state.disjointBoxLayout(),
                coarsened_fine.disjointBoxLayout(), grown_domain_box,
                m_num_ghosts * IntVect::Unit, m_num_ghosts * IntVect::Unit);
            a_coarse_state.copyTo(a_coarse_state.interval(), coarsened_fine,
                                  coarsened_fine.interval(), boundary_copier);

            // iterate through the coarse boxes, shared amongst threads
            DataIterator dit = coarsened_layout.dataIterator();
            int nbox = dit.size();
#pragma omp parallel for default(shared)
            for (int ibox = 0; ibox < nbox; ++ibox)
            {
                DataIndex dind = dit[ibox];
                FArrayBox &m_fine_box = a_fine_state[dind];
                FArrayBox &m_coarse_box = coarsened_fine[dind];
                Box this_box = m_coarse_box.box();
                Box fine_box = m_fine_box.box();
                IntVect offset_lo =
                    -this_box.smallEnd() + coarse_domain_box.smallEnd();
                IntVect offset_hi =
                    +this_box.bigEnd() - coarse_domain_box.bigEnd();

                // reduce box to the intersection of the box and the
                // problem domain ie remove all outer ghost cells
                this_box &= coarse_domain_box;

                // get the boundary box - remove one cell as we only want 2
                // coarse cells filled in each direction, to fill the 3 fine
                // cells on the level above
                Box boundary_box = get_boundary_box(a_side, idir, offset_lo,
                                                    offset_hi, this_box, 1);

                // define standard stencil for interp where not near
                // boundaries in other dirs
                IntVect default_offset =
                    IntVect::Zero + sign(a_side) * 2 * BASISV(idir);
                FourthOrderInterpStencil default_stencil(default_offset,
                                                         ref_ratio);

                // now interp the box from coarse to fine
                BoxIterator bit(boundary_box);
                for (bit.begin(); bit.ok(); ++bit)
                {
                    IntVect iv = bit();
                    IntVect lo_local_offset = iv - m_coarse_box.smallEnd();
                    IntVect hi_local_offset = m_coarse_box.bigEnd() - iv;

                    // bit of work to get the right stencils for near the
                    // edges of the box
                    bool near_boundary = false;
                    IntVect local_boundary_offset = IntVect::Zero;
                    for(int idir2; idir2 < SpaceDim; ++idir2)
                    {
                        if (idir2 == idir)
                        {
                            local_boundary_offset[idir2] =
                                default_offset[idir2];
                        }
                        else if ((idir2 != idir) &&
                                 (lo_local_offset[idir2] > 1) &&
                                 (hi_local_offset[idir2] > 1))
                        {
                            local_boundary_offset[idir2] = 0;
                        }
                        else if ((idir2 != idir) &&
                                 (lo_local_offset[idir2] == 1))
                        {
                            local_boundary_offset[idir2] = -2;
                            near_boundary = true;
                        }
                        else if ((idir2 != idir) &&
                                 (hi_local_offset[idir2] == 1))
                        {
                            local_boundary_offset[idir2] = +2;
                            near_boundary = true;
                        }
                        else
                        {
                            MayDay::Error(
                                "GRChomboBCs::define bad boxes");
                        }
                    }

                    // if not near the boundary use the default stencil,
                    // otherwise use the one calculated locally
                    if (!near_boundary)
                    {
                        default_stencil.fillFine(m_fine_box, m_coarse_box, iv);
                    }
                    else
                    {
                        FourthOrderInterpStencil local_stencil(
                            local_boundary_offset, ref_ratio);
                        local_stencil.fillFine(m_fine_box, m_coarse_box, iv);
                    }
                } // end loop box
            }     // end loop boxes
        }         // end if is_periodic
    }             // end loop idir
}

/// Get the boundary condition for a_dir and a_side
int GRChomboBCs::get_boundary_condition(const Side::LoHiSide a_side,
                                               const int a_dir)
{
    int boundary_condition = 0;
    if (a_side == Side::Lo)
    {
        boundary_condition = m_params.lo_boundary[a_dir];
    }
    else
    {
        boundary_condition = m_params.hi_boundary[a_dir];
    }
    return boundary_condition;
}

/// get the boundary box to fill if we are at a boundary
Box GRChomboBCs::get_boundary_box(
    const Side::LoHiSide a_side, const int a_dir, const IntVect &offset_lo,
    const IntVect &offset_hi, Box &this_ghostless_box, int shrink_for_coarse)
{
    // default constructor gives empty box
    Box boundary_box;

    // check if we are over the edges of the domain - are we a boundary box?
    // if so create the box of the cells we want to fill
    if (((a_side == Side::Hi) && (offset_hi[a_dir] > 0)) ||
        ((a_side == Side::Lo) && (offset_lo[a_dir] > 0)))
    {
        // Get just the boundary box to iterate over, m_num_ghosts ghost
        // cells unless we are filling the coarse cells in the interp case
        // where we want to fill only two coarse ghost cells (to cover 3
        // fine ones)
        if (a_side == Side::Lo)
        {
            boundary_box = adjCellLo(this_ghostless_box, a_dir,
                                     m_num_ghosts - shrink_for_coarse);
        }
        else
        {
            boundary_box = adjCellHi(this_ghostless_box, a_dir,
                                     m_num_ghosts - shrink_for_coarse);
        }

        // adjust for any offsets - catches the corners etc
        // but only want to fill them once, so y fills x, z fills y and x
        // etc. Required in periodic direction corners in cases where there
        // are mixed boundaries, (otherwise these corners are full of nans)
        for(int idir = 0; idir < SpaceDim; ++idir)
        {
            if (offset_lo[idir] > 0) // this direction is a low end boundary
            {
                if ((idir < a_dir) || (m_params.is_periodic[idir]))
                {
                    // grow it to fill the corners
                    boundary_box.growLo(idir, m_num_ghosts - shrink_for_coarse);
                }
            }
            else // cut off end ghost cell
            {
                if (idir != a_dir)
                {
                    boundary_box.growLo(idir, -shrink_for_coarse);
                }
            }

            if (offset_hi[idir] > 0) // this direction is a high end
                                     // boundary
            {
                if ((idir < a_dir) || (m_params.is_periodic[idir]))
                {
                    // grow it to fill the corners
                    boundary_box.growHi(idir, m_num_ghosts - shrink_for_coarse);
                }
            }
            else // cut off end ghost cell
            {
                if (idir != a_dir)
                {
                    boundary_box.growHi(idir, -shrink_for_coarse);
                }
            }
        }
    }
    return boundary_box;
}
