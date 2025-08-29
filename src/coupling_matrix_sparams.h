/*
 * coupling_matrix_sparams.h
 *
 *  Created on: 27 ago 2025
 *      Author: egonzalez
 */

#ifndef COUPLING_MATRIX_SPARAMS_H_
#define COUPLING_MATRIX_SPARAMS_H_

#include <complex.h>

/* Build coupling matrix M (NxN) from nearest-neighbour couplings */
void build_coupling_matrix_from_k(const double *k, int N, double **M);

void compute_sparams_from_M(int N, double **M, double f0, double BW, double Qe1, double QeN, double Qu, const double *freqs, int F, double complex *S11,
        double complex *S21);

#endif /* COUPLING_MATRIX_SPARAMS_H_ */
