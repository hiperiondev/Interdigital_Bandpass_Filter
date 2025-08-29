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

#ifndef COUPLING_MATRIX_SPARAMS_H_
#define COUPLING_MATRIX_SPARAMS_H_

#include <complex.h>

/* Build coupling matrix M (NxN) from nearest-neighbour couplings */
void build_coupling_matrix_from_k(const double *k, int N, double **M);

void compute_sparams_from_M(int N, double **M, double f0, double BW, double Qe1, double QeN, double Qu, const double *freqs, int F, double complex *S11,
        double complex *S21);

#endif /* COUPLING_MATRIX_SPARAMS_H_ */
