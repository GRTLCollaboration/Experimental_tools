/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef SETBINARYBH_HPP
#define SETBINARYBH_HPP

#include "PoissonParameters.H"

// Get the bh radius - input loc_bh should be location on the grid relative
// to the centre - note that it gets corrected to the actual location
// (see examples below)
Real get_bh_radius(RealVect &loc_bh, const Real bh_x_offset)
{
    loc_bh[0] -= bh_x_offset;
    return sqrt(loc_bh[0] * loc_bh[0] + loc_bh[1] * loc_bh[1] +
                loc_bh[2] * loc_bh[2]);
}

// Set Aij Bowen York data
// see Alcubierre pg 110 eqn (3.4.22)
Real get_Aij(const int i, const int j, const Real &rbh1, const Real &rbh2,
             const RealVect &n1, const RealVect &n2, const RealVect &J1,
             const RealVect &J2, const RealVect &P1, const RealVect &P2,
             const PoissonParameters &a_params)
{
    std::array<std::array<std::array<double, 3>, 3>, 3> epsilon = {0.};
    epsilon[0][1][2] = 1.0;
    epsilon[1][2][0] = 1.0;
    epsilon[2][0][1] = 1.0;
    epsilon[0][2][1] = -1.0;
    epsilon[2][1][0] = -1.0;
    epsilon[1][0][2] = -1.0;

    Real Aij = 1.5 / rbh1 / rbh1 * (n1[i] * P1[j] + n1[j] * P1[i]) +
               1.5 / rbh2 / rbh2 * (n2[i] * P2[j] + n2[j] * P2[i]);
    for (int k = 0; k < SpaceDim; k++)
    {
        Aij +=
            1.5 / rbh1 / rbh1 * (n1[i] * n1[j] - Real(i == j)) * P1[k] * n1[k] +
            1.5 / rbh2 / rbh2 * (n2[i] * n2[j] - Real(i == j)) * P2[k] * n2[k];

        for (int l = 0; l < SpaceDim; l++)
        {
            Aij += -3.0 / rbh1 / rbh1 / rbh1 *
                       (epsilon[i][l][k] * n1[j] + epsilon[j][l][k] * n1[i]) *
                       n1[l] * J1[k] -
                   3.0 / rbh2 / rbh2 / rbh2 *
                       (epsilon[i][l][k] * n2[j] + epsilon[j][l][k] * n2[i]) *
                       n2[l] * J2[k];
        }
    }

    return Aij;
}

void set_binary_bh_Aij(FArrayBox &multigrid_vars_box, const IntVect &iv,
                       const RealVect &loc, const PoissonParameters &a_params)
{
    // set the BH values - location
    RealVect loc_bh1 = loc;
    Real rbh1 = get_bh_radius(loc_bh1, a_params.bh1_offset);

    RealVect loc_bh2 = loc;
    Real rbh2 = get_bh_radius(loc_bh2, a_params.bh2_offset);

    RealVect n1 = {loc_bh1[0] / rbh1, loc_bh1[1] / rbh1, loc_bh1[2] / rbh1};
    RealVect n2 = {loc_bh2[0] / rbh2, loc_bh2[1] / rbh2, loc_bh2[2] / rbh2};

    // the Bowen York params
    RealVect J1 = {0.0, 0.0, a_params.bh1_spin};
    RealVect J2 = {0.0, 0.0, a_params.bh2_spin};
    RealVect P1 = {0.0, a_params.bh1_momentum, 0.0};
    RealVect P2 = {0.0, a_params.bh2_momentum, 0.0};

    // set the Aij vars
    multigrid_vars_box(iv, c_A11_0) +=
        get_Aij(0, 0, rbh1, rbh2, n1, n2, J1, J2, P1, P2, a_params);
    multigrid_vars_box(iv, c_A22_0) +=
        get_Aij(1, 1, rbh1, rbh2, n1, n2, J1, J2, P1, P2, a_params);
    multigrid_vars_box(iv, c_A33_0) +=
        get_Aij(2, 2, rbh1, rbh2, n1, n2, J1, J2, P1, P2, a_params);
    multigrid_vars_box(iv, c_A12_0) +=
        get_Aij(0, 1, rbh1, rbh2, n1, n2, J1, J2, P1, P2, a_params);
    multigrid_vars_box(iv, c_A13_0) +=
        get_Aij(0, 2, rbh1, rbh2, n1, n2, J1, J2, P1, P2, a_params);
    multigrid_vars_box(iv, c_A23_0) +=
        get_Aij(1, 2, rbh1, rbh2, n1, n2, J1, J2, P1, P2, a_params);
}

Real set_binary_bh_psi(const RealVect &loc, const PoissonParameters &a_params)
{
    // the Bowen York params
    Real m1 = a_params.bh1_bare_mass;
    Real m2 = a_params.bh2_bare_mass;

    // set the BH values - location
    RealVect loc_bh1 = loc;
    Real rbh1 = get_bh_radius(loc_bh1, a_params.bh1_offset);

    RealVect loc_bh2 = loc;
    Real rbh2 = get_bh_radius(loc_bh2, a_params.bh2_offset);

    return m1 / rbh1 + m2 / rbh2;
}

void set_Aij_0(FArrayBox &multigrid_vars_box, const IntVect &iv,
               const RealVect &loc, const RealVect &a_dx,
               const PoissonParameters &a_params,
               const FArrayBox &a_grad_multigrid)
{

    Real trace = (a_grad_multigrid(iv, 3 * c_V0_0 + 0) +
                  a_grad_multigrid(iv, 3 * c_V1_0 + 1) +
                  a_grad_multigrid(iv, 3 * c_V2_0 + 2));

    multigrid_vars_box(iv, c_A11_0) =
        3. / 4. *
            (a_grad_multigrid(iv, 3 * c_V0_0 + 0) +
             a_grad_multigrid(iv, 3 * c_V0_0 + 0)) -
        1. / 2. * trace;

    multigrid_vars_box(iv, c_A22_0) =
        3. / 4. *
            (a_grad_multigrid(iv, 3 * c_V1_0 + 1) +
             a_grad_multigrid(iv, 3 * c_V1_0 + 1)) -
        1. / 2. * trace;

    multigrid_vars_box(iv, c_A33_0) =
        3. / 4. *
            (a_grad_multigrid(iv, 3 * c_V2_0 + 2) +
             a_grad_multigrid(iv, 3 * c_V2_0 + 2)) -
        1. / 2. * trace;

    multigrid_vars_box(iv, c_A12_0) = 3. / 4. *
                                      (a_grad_multigrid(iv, 3 * c_V0_0 + 1) +
                                       a_grad_multigrid(iv, 3 * c_V1_0 + 0));

    multigrid_vars_box(iv, c_A13_0) = 3. / 4. *
                                      (a_grad_multigrid(iv, 3 * c_V0_0 + 2) +
                                       a_grad_multigrid(iv, 3 * c_V2_0 + 0));

    multigrid_vars_box(iv, c_A23_0) = 3. / 4. *
                                      (a_grad_multigrid(iv, 3 * c_V2_0 + 1) +
                                       a_grad_multigrid(iv, 3 * c_V1_0 + 2));
}

Real get_ricci(FArrayBox &multigrid_vars_box, const IntVect &iv,
               const RealVect &loc, const RealVect &a_dx,
               const PoissonParameters &a_params,
               const FArrayBox &a_grad_multigrid,
               const FArrayBox &a_grad2_multigrid,
               const FArrayBox &a_mixed_grad2_multigrid,
               const FArrayBox &a_grad_h_UU)
{

    Real ricci_tensor[3][3] = {0};
    Real h[3][3];
    h[0][0] = multigrid_vars_box(iv, c_h11_0);
    h[1][1] = multigrid_vars_box(iv, c_h22_0);
    h[2][2] = multigrid_vars_box(iv, c_h33_0);
    h[0][1] = multigrid_vars_box(iv, c_h12_0);
    h[0][2] = multigrid_vars_box(iv, c_h13_0);
    h[1][2] = multigrid_vars_box(iv, c_h23_0);
    h[0][1] = h[1][0];
    h[0][2] = h[2][0];
    h[1][2] = h[2][1];

    Real det = h[0][0] * (h[1][1] * h[2][2] - h[1][2] * h[2][1]) -
               h[0][1] * (h[2][2] * h[1][0] - h[1][2] * h[2][0]) +
               h[0][2] * (h[1][0] * h[2][1] - h[1][1] * h[2][0]);

    Real h_UU[3][3];
    Real det_inverse = 1. / det;
    h_UU[0][0] = (h[1][1] * h[2][2] - h[1][2] * h[2][1]) * det_inverse;
    h_UU[1][1] = (h[0][0] * h[2][2] - h[0][2] * h[2][0]) * det_inverse;
    h_UU[2][2] = (h[0][0] * h[1][1] - h[1][0] * h[0][1]) * det_inverse;
    h_UU[0][1] = (h[2][0] * h[1][2] - h[1][0] * h[2][2]) * det_inverse;
    h_UU[1][0] = (h[0][2] * h[2][1] - h[0][1] * h[2][2]) * det_inverse;
    h_UU[0][2] = (h[1][0] * h[2][1] - h[1][1] * h[2][0]) * det_inverse;
    h_UU[2][0] = (h[0][1] * h[1][2] - h[1][1] * h[0][2]) * det_inverse;
    h_UU[1][2] = (h[0][1] * h[2][0] - h[0][0] * h[2][1]) * det_inverse;
    h_UU[2][1] = (h[1][0] * h[0][2] - h[0][0] * h[1][2]) * det_inverse;

    Real d1_h[3][3][3]= {0};
    Real d1_hUU[3][3][3] = {0};
    Real d2_h[3][3][3][3] = {0};

    for (int i = 0; i < SpaceDim; i++)
    {
        d1_h[0][0][i] = a_grad_multigrid(iv, 3 * c_h11_0 + i);
        d1_h[1][1][i] = a_grad_multigrid(iv, 3 * c_h22_0 + i);
        d1_h[2][2][i] = a_grad_multigrid(iv, 3 * c_h33_0 + i);
        d1_h[0][1][i] = a_grad_multigrid(iv, 3 * c_h12_0 + i);
        d1_h[0][2][i] = a_grad_multigrid(iv, 3 * c_h13_0 + i);
        d1_h[1][2][i] = a_grad_multigrid(iv, 3 * c_h23_0 + i);
        d1_h[1][0][i] = d1_h[0][1][i];
        d1_h[2][0][i] = d1_h[0][2][i];
        d1_h[2][1][i] = d1_h[1][2][i];

        d2_h[0][0][i][i] = a_grad2_multigrid(iv, 3 * c_h11_0 + i);
        d2_h[1][1][i][i] = a_grad2_multigrid(iv, 3 * c_h22_0 + i);
        d2_h[2][2][i][i] = a_grad2_multigrid(iv, 3 * c_h33_0 + i);
        d2_h[0][1][i][i] = a_grad2_multigrid(iv, 3 * c_h12_0 + i);
        d2_h[0][2][i][i] = a_grad2_multigrid(iv, 3 * c_h13_0 + i);
        d2_h[1][2][i][i] = a_grad2_multigrid(iv, 3 * c_h23_0 + i);
        d2_h[1][0][i][i] = d2_h[0][1][i][i];
        d2_h[2][0][i][i] = d2_h[0][2][i][i];
        d2_h[2][1][i][i] = d2_h[1][2][i][i];

        d1_hUU[0][0][i] = a_grad_h_UU(iv, 3 * 0 + i);
        d1_hUU[1][1][i] = a_grad_h_UU(iv, 3 * 1 + i);
        d1_hUU[2][2][i] = a_grad_h_UU(iv, 3 * 2 + i);
        d1_hUU[0][1][i] = a_grad_h_UU(iv, 3 * 3 + i);
        d1_hUU[0][2][i] = a_grad_h_UU(iv, 3 * 4 + i);
        d1_hUU[1][2][i] = a_grad_h_UU(iv, 3 * 5 + i);
        d1_hUU[1][0][i] = d1_hUU[0][1][i];
        d1_hUU[2][0][i] = d1_hUU[0][2][i];
        d1_hUU[2][1][i] = d1_hUU[1][2][i];

        for (int j = i+1; i < SpaceDim; i++)
        {
            d2_h[0][0][i][j] = a_mixed_grad2_multigrid(iv, 3 * c_h11_0 + i + j);
            d2_h[1][1][i][j] = a_mixed_grad2_multigrid(iv, 3 * c_h22_0 + i + j);
            d2_h[2][2][i][j] = a_mixed_grad2_multigrid(iv, 3 * c_h33_0 + i + j);
            d2_h[0][1][i][j] = a_mixed_grad2_multigrid(iv, 3 * c_h12_0 + i + j);
            d2_h[0][2][i][j] = a_mixed_grad2_multigrid(iv, 3 * c_h13_0 + i + j);
            d2_h[1][2][i][j] = a_mixed_grad2_multigrid(iv, 3 * c_h23_0 + i + j);
            d2_h[1][0][i][j] = d2_h[0][1][i][j];
            d2_h[2][0][i][j] = d2_h[0][2][i][j];
            d2_h[2][1][i][j] = d2_h[1][2][i][j];

            d2_h[0][0][j][i] = d2_h[0][0][i][j];
            d2_h[1][1][j][i] = d2_h[1][1][i][j];
            d2_h[2][2][j][i] = d2_h[2][2][i][j];
            d2_h[0][1][j][i] = d2_h[0][1][i][j];
            d2_h[0][2][j][i] = d2_h[0][2][i][j];
            d2_h[1][2][j][i] = d2_h[1][2][i][j];
            d2_h[1][0][j][i] = d2_h[0][1][i][j];
            d2_h[2][0][j][i] = d2_h[0][2][i][j];
            d2_h[2][1][j][i] = d2_h[1][2][i][j];
        }
    }


    for (int i = 0; i < SpaceDim; i++)
    {
        for (int j = 0; j < SpaceDim; j++)
        {
            for (int k = 0; k < SpaceDim; k++)
            {
                for (int l = 0; l < SpaceDim; l++)
                {
                    ricci_tensor[i][j] +=
                        0.5 * d1_hUU[k][l][k] *
                            (d1_h[j][l][i] + d1_h[i][k][j] - d1_h[i][j][l]) +
                        0.5 * h_UU[k][l] *
                            (d2_h[i][l][j][k] + d2_h[j][l][i][k] -
                             d2_h[i][j][k][l] - d2_h[k][l][i][j]) -
                        0.5 * d1_hUU[k][l][j] * d1_h[k][l][i];

                    for (int m = 0; m < SpaceDim; m++)
                    {
                        for (int n = 0; n < SpaceDim; n++)
                        {
                            ricci_tensor[i][j] +=
                                0.25 * h_UU[m][n] * h_UU[k][l] *
                                (d1_h[k][l][n] *
                                     (d1_h[i][m][j] + d1_h[j][m][i] -
                                      d1_h[i][j][m]) -
                                 d1_h[m][k][i] * d1_h[n][l][j] -
                                 2 * d1_h[i][k][m] * d1_h[j][n][l] +
                                 2 * d1_h[i][k][m] * d1_h[j][l][n]);
                        }
                    }
                }
            }
        }
    }
    Real out_ricci = 0;

    for (int i = 0; i < SpaceDim; i++)
    {
        for (int j = 0; j < SpaceDim; j++)
        {
            out_ricci += h_UU[i][j] * ricci_tensor[i][j];
        }
    }

    return out_ricci = 0;
}
#endif /* SETBINARYBH_HPP */
