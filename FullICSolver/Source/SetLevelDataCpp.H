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

// computes the second derivative of the variable you want in a box, i.e. d_i d_i phi
inline RealVect get_grad_2(const Box &a_b, const FArrayBox &multigrid_vars_box, const int &c_comp,
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

		RealVect grad_2 = {0.0, 0.0, 0.0};
	    for (int idir = 0; idir < SpaceDim; ++idir)
  		{
			IntVect iv_offset0 = a_iv;
        	IntVect iv_offset1 = a_iv;
        	IntVect iv_offset2 = a_iv;
        	iv_offset1[idir] -= 1;
        	iv_offset2[idir] += 1;

        	// 2nd order stencils for now
        	grad_2[idir] =
            	 (multigrid_vars_box(iv_offset2, c_comp) - 2.*multigrid_vars_box(iv_offset0, c_comp) + multigrid_vars_box(iv_offset1, c_comp)) / (a_dx[idir] * a_dx[idir]);
			
			a_grad_of_variable(a_iv, 3*c_comp + idir) = grad_2[idir]; 
		}

	}
} // end get_grad2

// computes the mixed second derivatives of the variable you want in a box, i.e. d_i d_j phi
inline RealVect get_grad_mixed(const Box &a_b, const FArrayBox &multigrid_vars_box, const int &c_comp,
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

		RealVect grad_mixed = {0.0, 0.0, 0.0};
	    for (int idir = 0; idir < SpaceDim; ++idir)
  		{
			for (int idir2 = 0; idir2 < idir + 1; ++idir2)
			{
				if (idir != idir2)
				{
					IntVect iv_offset1 = a_iv;
		        	IntVect iv_offset2 = a_iv;
		        	IntVect iv_offset3 = a_iv;
					IntVect iv_offset4 = a_iv;
		        	iv_offset1[idir]  -= 1;
					iv_offset1[idir2] -= 1;
		         	iv_offset2[idir]  -= 1;
					iv_offset2[idir2] += 1;
		        	iv_offset3[idir]  += 1;
			     	iv_offset3[idir2] -= 1;
					iv_offset4[idir]  += 1;
		    		iv_offset4[idir2] += 1;

		        	// 2nd order stencils for now
		        	grad_mixed[idir+idir2-1] =
		            	 (multigrid_vars_box(iv_offset4, c_comp) - multigrid_vars_box(iv_offset3, c_comp) - multigrid_vars_box(iv_offset2, c_comp) 
						  										 + multigrid_vars_box(iv_offset1, c_comp)) / (4.*a_dx[idir] * a_dx[idir2]);
				}
			}
		}

		for (int idir = 0; idir < SpaceDim; ++idir)
		{		
			a_grad_of_variable(a_iv, 3*c_comp + idir) = grad_mixed[idir]; 
		}

	}
} // end get_grad_mixed


// computes the laplacian of the variable you want in a box, i.e. d_i d^i phi
inline RealVect get_laplacian(const Box &a_b, const FArrayBox &multigrid_vars_box, const int &c_comp,
                            const RealVect &a_dx, FArrayBox &a_lapl_of_variable, const PoissonParameters &a_params)
{
	BoxIterator bit(a_b);
	for (bit.begin(); bit.ok(); ++bit)
    {
        IntVect a_iv = bit();
        // set the phi value - need the distance from centre
        RealVect loc(a_iv + 0.5 * RealVect::Unit);
        loc *= a_dx;
    	loc -= a_params.domainLength / 2.0;

		Real laplacian_of_psi = 0.0;
    	for (int idir = 0; idir < SpaceDim; ++idir)
    	{   
        	IntVect iv_offset1 = a_iv;
        	IntVect iv_offset2 = a_iv;
        	iv_offset1[idir] -= 1;
        	iv_offset2[idir] += 1;
    
        	// 2nd order stencil for now
        	Real d2psi_dxdx = 1.0 / (a_dx[idir] * a_dx[idir]) *
            	              (+1.0 * multigrid_vars_box(iv_offset2, c_comp) -
                	           2.0 * multigrid_vars_box(a_iv, c_comp) + 1.0 * multigrid_vars_box(iv_offset1, c_comp));
        	laplacian_of_psi += d2psi_dxdx;
		}
		a_lapl_of_variable(a_iv, c_comp) = laplacian_of_psi;
	}
} // end get_grad
