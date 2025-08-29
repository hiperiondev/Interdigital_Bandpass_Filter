/*
 * coupling_matrix_sparams.c
 *
 *  Created on: 27 ago 2025
 *      Author: egonzalez
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

/* Build coupling matrix M (NxN) from nearest-neighbour couplings */
void build_coupling_matrix_from_k(const double *k, int N, double **M) {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            M[i][j] = 0.0;
        }
    }

    for (int i = 0; i < N - 1; ++i) {
        double kij = fabs(k[i]); /* keep magnitude input safe */
        double sgn = (i % 2 == 0) ? +1.0 : -1.0; /* alternate sign for interdigital */
        double v = sgn * kij; /* apply sign pattern */
        M[i][i + 1] = v; /* signed coupling i <-> i+1 */
        M[i + 1][i] = v; /* enforce symmetry */
    }
}

/* Compute S-parameters from coupling matrix */
void compute_sparams_from_M(int N, double **M, double f0, double BW, /* added BW argument */
double Qe1, double QeN, double Qu, const double *freqs, int F, double complex *S11, double complex *S21) {
    int NA = N + 2; /* augmented size: source + N + load */
    double FBW = BW / f0; /* fractional bandwidth */

    /* allocate dynamic arrays (complex) */
    double complex *A = (double complex*) calloc((size_t) NA * (size_t) NA, sizeof(double complex));
    double complex *Inv = (double complex*) calloc((size_t) NA * (size_t) NA, sizeof(double complex));
    double complex *Work = (double complex*) calloc((size_t) NA * (size_t) (2 * NA), sizeof(double complex));

    if (!A || !Inv || !Work) {
        fprintf(stderr, "Allocation failed in compute_sparams_from_M\n");
        free(A);
        free(Inv);
        free(Work);
        return;
    }

    /* Build frequency-independent part: M_aug */
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            A[(i + 1) * NA + (j + 1)] = M[i][j];
        }
    }

    /* External couplings */
    double ms = (Qe1 > 0.0) ? (1.0 / sqrt(Qe1)) : 0.0;
    double ml = (QeN > 0.0) ? (1.0 / sqrt(QeN)) : 0.0;
    A[0 * NA + 1] = ms;
    A[1 * NA + 0] = ms;
    A[N * NA + (N + 1)] = ml;
    A[(N + 1) * NA + N] = ml;

    /* Loss factor */
    double dloss = (Qu > 0.0) ? (1.0 / (2.0 * Qu)) : 0.0;

    for (int t = 0; t < F; ++t) {
        double f = freqs[t];

        /* correct normalized low-pass variable using BW */
        double Omega = (f * f - f0 * f0) / (f * f0 * FBW);

        /* Form A(ω) as augmented matrix [A(ω) | I] */
        for (int r = 0; r < NA; ++r) {
            for (int c = 0; c < NA; ++c) {
                Work[r * (2 * NA) + c] = A[r * NA + c];
            }
        }

        for (int r = 0; r < NA; ++r) {
            for (int c = 0; c < NA; ++c) {
                Work[r * (2 * NA) + (NA + c)] = (r == c) ? 1.0 + 0.0 * I : 0.0 + 0.0 * I;
            }

            if (r == 0 || r == NA - 1) {
                Work[r * (2 * NA) + r] += -I; /* -j*G on ports */
            } else {
                Work[r * (2 * NA) + r] += (Omega + I * dloss); /* Ω + j*loss on resonators */
            }
        }

        /* -------- Gauss-Jordan inversion -------- */
        for (int col = 0; col < NA; ++col) {
            int piv = col;
            double best = cabs(Work[piv * (2 * NA) + col]);
            for (int r = col + 1; r < NA; ++r) {
                double mag = cabs(Work[r * (2 * NA) + col]);
                if (mag > best) {
                    best = mag;
                    piv = r;
                }
            }
            if (best == 0.0) {
                S11[t] = NAN + NAN * I;
                S21[t] = NAN + NAN * I;
                goto next_frequency;
            }
            if (piv != col) {
                for (int c = 0; c < 2 * NA; ++c) {
                    double complex tmp = Work[col * (2 * NA) + c];
                    Work[col * (2 * NA) + c] = Work[piv * (2 * NA) + c];
                    Work[piv * (2 * NA) + c] = tmp;
                }
            }
            double complex pivval = Work[col * (2 * NA) + col];
            for (int c = 0; c < 2 * NA; ++c) {
                Work[col * (2 * NA) + c] /= pivval;
            }
            for (int r = 0; r < NA; ++r) {
                if (r == col)
                    continue;
                double complex factor = Work[r * (2 * NA) + col];
                if (factor != 0.0 + 0.0 * I) {
                    for (int c = 0; c < 2 * NA; ++c) {
                        Work[r * (2 * NA) + c] -= factor * Work[col * (2 * NA) + c];
                    }
                }
            }
        }

        /* Extract inverse */
        for (int r = 0; r < NA; ++r) {
            for (int c = 0; c < NA; ++c) {
                Inv[r * NA + c] = Work[r * (2 * NA) + (NA + c)];
            }
        }

        /* Compute S11 and S21 */
        S11[t] = 1.0 + 2.0 * I * Inv[0 * NA + 0];
        S21[t] = -2.0 * I * Inv[(NA - 1) * NA + 0];

        next_frequency:
        ;
    }

    free(A);
    free(Inv);
    free(Work);
}
