/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

// computes the gradient of the variable you want in a box, i.e. d_i phi
inline RealVect get_grad(const Box &a_b, const FArrayBox &multigrid_vars_box, const int &c_comp,
                            const RealVect &a_dx, FArrayBox &a_grad_of_variable, const PoissonParameters &a_params)
{
	BoxIterator bit(a_b);
	for (bit.begin(); bit.ok(); ++bit)
    {
        IntVect a_iv = bit();
        // set the phi value - need the distance from centre
        RealVect loc(a_iv + 0.5 * RealVect::Unit);
        loc *= a_dx;
    	loc -= a_params.domainLength / 2.0;

		RealVect grad = {0.0, 0.0, 0.0};
	    for (int idir = 0; idir < SpaceDim; ++idir)
  		{
        	IntVect iv_offset1 = a_iv;
        	IntVect iv_offset2 = a_iv;
        	iv_offset1[idir] -= 1;
        	iv_offset2[idir] += 1;

        	// 2nd order stencils for now
        	grad[idir] =
            	0.5 * (multigrid_vars_box(iv_offset2, c_comp) - multigrid_vars_box(iv_offset1, c_comp)) / a_dx[idir];
			
			a_grad_of_variable(a_iv, 3*c_comp + idir) = grad[idir]; 
		}

	}
} // end get_grad
