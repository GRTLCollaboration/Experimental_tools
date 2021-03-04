/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef MYPOTENTIALFUNCTION_HPP
#define MYPOTENTIALFUNCTION_HPP

Real my_potential_function(Real phi_here) {
  Real V_of_phi = 1.0;//1e-4*pow(tanh(phi_here), 2.0);

  return V_of_phi;
}

#endif /* MYPOTENTIALFUNCTION_HPP */
