/*
 * Copyright 2025 Emiliano Augusto Gonzalez (egonzalez . hiperion @ gmail . com))
 * * Project Site: https://github.com/hiperiondev/Interdigital_Bandpass_Filter *
 *
 * This is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 *
 */

#ifndef CALCULATIONS_H_
#define CALCULATIONS_H_

#include <complex.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define LINT_FACTOR   0.889   // interior rod length factor relative to lambda/4
#define LEND_FACTOR   1.005   // end rod length multiplier relative to interior
#define K_QU_DEFAULT  0.06    // empirical Qu constant (conservative default)

double toinch(double mm);
void butterworth_g(int n, double *g);
void chebyshev_g(int n, double ripple_dB, double *g);
double clamp(double x, double lo, double hi);
double interdigital_round_rod_length_corrected(double f0, double d, double h, double g, int keep_legacy_factor, double legacy_factor_value);
double rod_length(double f0, double h, double d, double length_factor);
void build_coupling_matrix_from_k(const double *k, int N, double **M);
void compute_sparams_from_M(int N, double **M, double f0, double BW, double Qe1, double QeN, double Qu, const double *freqs, int F, double complex *S11,
        double complex *S21);

#endif /* CALCULATIONS_H_ */
