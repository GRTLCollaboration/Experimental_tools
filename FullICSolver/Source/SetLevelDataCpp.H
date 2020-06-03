/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// computes the gradient of the variable you want in a box, i.e. d_i phi
inline RealVect get_grad(const Box &a_b, const FArrayBox &multigrid_vars_box,
                         const Interval &a_ncomp, const RealVect &a_dx,
                         FArrayBox &a_grad_of_variable,
                         const PoissonParameters &a_params)
{
    BoxIterator bit(a_b);
    for (bit.begin(); bit.ok(); ++bit)
    {
        IntVect a_iv = bit();
        // set the phi value - need the distance from centre
        RealVect loc(a_iv + 0.5 * RealVect::Unit);
        loc *= a_dx;
        loc -= a_params.domainLength / 2.0;

        for (int icomp = a_ncomp.begin(); icomp <= a_ncomp.end(); ++icomp)
        {
            RealVect grad = {0.0, 0.0, 0.0};
            for (int idir = 0; idir < SpaceDim; ++idir)
            {
                IntVect iv_offset1 = a_iv;
                IntVect iv_offset2 = a_iv;
                iv_offset1[idir] -= 1;
                iv_offset2[idir] += 1;

                // 2nd order stencils for now
                grad[idir] = 0.5 *
                             (multigrid_vars_box(iv_offset2, icomp) -
                              multigrid_vars_box(iv_offset1, icomp)) /
                             a_dx[idir];

                a_grad_of_variable(a_iv, 3 * icomp + idir) = grad[idir];
                //				pout() <<
                // a_grad_of_variable(a_iv, 3*icomp + idir) << endl;
            }
        }
    }
} // end get_grad

// computes the second derivative of the variable you want in a box, i.e. d_i
// d_i phi
inline RealVect get_grad2(const Box &a_b, const FArrayBox &multigrid_vars_box,
                          const Interval &a_ncomp, const RealVect &a_dx,
                          FArrayBox &a_grad2_of_variable,
                          const PoissonParameters &a_params)
{
    BoxIterator bit(a_b);
    for (bit.begin(); bit.ok(); ++bit)
    {
        IntVect a_iv = bit();
        // set the phi value - need the distance from centre
        RealVect loc(a_iv + 0.5 * RealVect::Unit);
        loc *= a_dx;
        loc -= a_params.domainLength / 2.0;

        for (int icomp = a_ncomp.begin(); icomp <= a_ncomp.end(); ++icomp)
        {
            RealVect grad_2 = {0.0, 0.0, 0.0};
            for (int idir = 0; idir < SpaceDim; ++idir)
            {
                IntVect iv_offset0 = a_iv;
                IntVect iv_offset1 = a_iv;
                IntVect iv_offset2 = a_iv;
                iv_offset1[idir] -= 1;
                iv_offset2[idir] += 1;

                // 2nd order stencils for now
                grad_2[idir] = (multigrid_vars_box(iv_offset2, icomp) -
                                2. * multigrid_vars_box(iv_offset0, icomp) +
                                multigrid_vars_box(iv_offset1, icomp)) /
                               (a_dx[idir] * a_dx[idir]);

                a_grad2_of_variable(a_iv, 3 * icomp + idir) = grad_2[idir];
            }
        }
    }
} // end get_grad2

// computes the mixed derivative of the variable you want in a box, i.e. d_i d_j
// phi
inline RealVect get_mixed_grad2(const Box &a_b,
                                const FArrayBox &multigrid_vars_box,
                                const Interval &a_ncomp, const RealVect &a_dx,
                                FArrayBox &a_mixed_grad2_of_variable,
                                const PoissonParameters &a_params)
{
    BoxIterator bit(a_b);
    for (bit.begin(); bit.ok(); ++bit)
    {
        IntVect a_iv = bit();
        // set the phi value - need the distance from centre
        RealVect loc(a_iv + 0.5 * RealVect::Unit);
        loc *= a_dx;
        loc -= a_params.domainLength / 2.0;

        for (int icomp = a_ncomp.begin(); icomp <= a_ncomp.end(); ++icomp)
        {
            RealVect mixed_grad2 = {0.0, 0.0, 0.0};
            for (int idir1 = 0; idir1 < SpaceDim; ++idir1)
            {
                for (int idir2 = idir1 + 1; idir2 < SpaceDim; ++idir2)
                {
                    IntVect iv_offset1 = a_iv;
                    IntVect iv_offset2 = a_iv;
                    IntVect iv_offset3 = a_iv;
                    IntVect iv_offset4 = a_iv;
                    iv_offset1[idir1] -= 1;
                    iv_offset1[idir2] -= 1;
                    iv_offset2[idir1] += 1;
                    iv_offset2[idir2] += 1;
                    iv_offset3[idir1] += 1;
                    iv_offset3[idir2] -= 1;
                    iv_offset4[idir1] -= 1;
                    iv_offset4[idir2] += 1;

                    // 2nd order stencils for now
                    mixed_grad2[idir1 + idir2] =
                        (multigrid_vars_box(iv_offset1, icomp) +
                         multigrid_vars_box(iv_offset2, icomp) -
                         multigrid_vars_box(iv_offset3, icomp) -
                         multigrid_vars_box(iv_offset4, icomp)) /
                        (4 * a_dx[idir1] * a_dx[idir2]);

                    a_mixed_grad2_of_variable(a_iv, 3 * icomp + idir1 +idir2) =
                        mixed_grad2[idir1 + idir2];
                }
            }
        }
    }
} // end get_mixed_grad2

// computes the laplacian of the variable you want in a box, i.e. d_i d^i phi
inline RealVect get_laplacian(const Box &a_b, // Should this be RealVect or Real
                              const FArrayBox &multigrid_vars_box,
                              const Interval &a_ncomp, const RealVect &a_dx,
                              FArrayBox &a_lapl_of_variable,
                              const PoissonParameters &a_params)
{
    BoxIterator bit(a_b);
    for (bit.begin(); bit.ok(); ++bit)
    {
        IntVect a_iv = bit();
        // set the phi value - need the distance from centre
        RealVect loc(a_iv + 0.5 * RealVect::Unit);
        loc *= a_dx;
        loc -= a_params.domainLength / 2.0;

        for (int icomp = a_ncomp.begin(); icomp <= a_ncomp.end(); ++icomp)
        {
            Real laplacian_of_psi = 0.0;
            for (int idir = 0; idir < SpaceDim; ++idir)
            {
                IntVect iv_offset1 = a_iv;
                IntVect iv_offset2 = a_iv;
                iv_offset1[idir] -= 1;
                iv_offset2[idir] += 1;

                // 2nd order stencil for now
                Real d2psi_dxdx =
                    1.0 / (a_dx[idir] * a_dx[idir]) *
                    (+1.0 * multigrid_vars_box(iv_offset2, icomp) -
                     2.0 * multigrid_vars_box(a_iv, icomp) +
                     1.0 * multigrid_vars_box(iv_offset1, icomp));
                laplacian_of_psi += d2psi_dxdx;
            }
            a_lapl_of_variable(a_iv, icomp) = laplacian_of_psi;
        }
    }
} // end get_grad
