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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

/* Physical constants */
#ifndef C0
#define C0 299792458.0
#endif

double toinch(double mm) {
    return mm * 0.03937007874015748;
}

/* Butterworth prototype g-values (ripple = 0). g0=g_{n+1}=1. */
void butterworth_g(int n, double *g) {
    g[0] = 1.0;
    for (int k = 1; k <= n; k++) {
        g[k] = 2.0 * sin((2.0 * k - 1.0) * M_PI / (2.0 * n));
    }
    g[n + 1] = 1.0;
}

/* Chebyshev prototype g-values (equal ripple r dB).
 Reference: standard closed-form from many texts/Wikipedia. */
void chebyshev_g(int n, double ripple_dB, double *g) {
    /* Compute epsilon and beta (standard definition) */
    double eps = sqrt(pow(10.0, ripple_dB / 10.0) - 1.0);
    double beta = asinh(1.0 / eps) / n;

    /* Base term */
    g[0] = 1.0;

    /* Recurrence arrays (sized to safely handle n <= 30) */
    double sinhb = sinh(beta);
    double a[64] = { 0 }, b[64];

    /* compute a_k */
    for (int k = 1; k <= n; k++) {
        a[k] = sin((2.0 * k - 1.0) * M_PI / (2.0 * n));
    }

    /* compute b_k = sinh^2(beta) + sin^2(k*pi/n) for k=1..n-1 */
    for (int k = 1; k <= n - 1; k++) {
        double s = sin(k * M_PI / n);
        b[k] = sinhb * sinhb + s * s;
    }

    /* g1 and recurrence for gk */
    g[1] = (2.0 * a[1]) / sinhb;
    for (int k = 2; k <= n; k++) {
        g[k] = (4.0 * a[k - 1] * a[k]) / (b[k - 1] * g[k - 1]);
    }

    if ((n % 2) == 0) {
        g[n + 1] = 1.0;
    } else {
        double tb2 = tanh(0.5 * beta);
        g[n + 1] = 1.0 / (tb2 * tb2);
    }
}

double clamp(double x, double lo, double hi) {
    return (x < lo) ? lo : (x > hi) ? hi : x;
}

/* --------------------------------------------------------------------------
 Equivalent strip width for a round rod in a two-lid (stripline-like) cavity.

 We need a width W_eq so we can use an open-end length extension that is
 tabulated for strips/bars. A pragmatic mapping is to set W_eq to roughly
 half the rod circumference: W_eq ≈ (π * d) / 2. This keeps the local
 fringing “aperture” of a rod comparable to a strip of width W_eq.

 This is an engineering approximation to feed the end-effect model; it
 consistently trends correctly with d/h and produces good first-pass
 results before EM/tuning. If you prefer a different mapping, swap it here.
 -------------------------------------------------------------------------- */
static inline double round_rod_equivalent_strip_width(double d) {
    return 0.5 * M_PI * d; /* W_eq ≈ π d / 2 */
}

/* --------------------------------------------------------------------------
 Open-end length extension for a (wide) line between two lids (air-filled).
 We reuse the widely-used Silvester–Benedek/Hammerstad-style microstrip
 open-end expression with ε_eff≈1, but we feed it with W_eq and h from our
 stripline-like cavity. This captures the main dependence on W/h and h.

 Δℓ/h ≈ 0.412 * ((ε_eff + 0.3)/(ε_eff - 0.258)) * ((W/h + 0.262)/(W/h + 0.813))

 For air-filled between two lids, use ε_eff ≈ 1 (fields largely in air).
 Source of the formula pattern: Silvester & Benedek; Hammerstad refinements.
 (We cite a teaching note that reproduces the equation.)
 -------------------------------------------------------------------------- */
static double open_end_delta_l_from_w_h(double w, double h) {
    const double er_eff = 1.0; /* air-filled */
    const double wh = clamp(w / h, 1e-6, 1e6);
    const double num1 = (er_eff + 0.3);
    const double den1 = (er_eff - 0.258);
    /* Guard against division by near-zero if someone forces er_eff≈0.258 */
    double factor1 = (fabs(den1) < 1e-6) ? 1.0 : (num1 / den1);
    double factor2 = (wh + 0.262) / (wh + 0.813);
    double dlh = 0.412 * factor1 * factor2;
    if (dlh < 0.0)
        dlh = 0.0;
    return dlh * h; /* meters */
}

/* Convenience wrapper for round rods */
static double open_end_delta_l_round_rod(double d, double h) {
    const double w_eq = round_rod_equivalent_strip_width(d);
    return open_end_delta_l_from_w_h(w_eq, h);
}

/* --------------------------------------------------------------------------
 Extra end-loading from the adjacent interdigital end gap (heuristic).

 The adjacent grounded (at the far end) rod forms an *electric* coupling
 capacitor at the open ends. This additional shunt-C increases electrical
 length, so we translate it into an extra length extension Δℓ_gap.

 We scale Δℓ_gap as a fraction of Δℓ_open, stronger when the end gap is
 small compared to the rod diameter. This tracks the well-known increase
 in end loading as g shrinks.

 Δℓ_gap ≈ α * Δℓ_open * 1/(1 + g/d)
 with α in [0.2, 0.5] depending on how “interdigital” the ends really are.
 We pick α = 0.35 as a conservative default for round-rod interdigital.

 Rationale:
 - Matthaei/Hong handle this rigorously via fringing/mutual-capacitance
 charts; our heuristic just biases Δℓ upward when g is small.
 - Safe as a *first-pass*; refine α once you calibrate to your hardware/EM.
 -------------------------------------------------------------------------- */
static double end_gap_delta_l_round_rod(double d, double g, double h) {
    const double alpha = 0.35; /* tune 0.2–0.5 if you have measurements */
    const double d_open = open_end_delta_l_round_rod(d, h);
    const double ratio = d > 0.0 ? (g / d) : 1e9;
    const double scale = 1.0 / (1.0 + ratio);
    return alpha * d_open * scale;
}

/* --------------------------------------------------------------------------
 Master: compute physical rod length for target f0 with geometry awareness.
 Inputs:
 f0  : resonant frequency [Hz]
 d   : rod diameter [m]
 h   : lid spacing (ground-to-ground) [m]
 g   : end gap to adjacent rod [m] (use the *smallest* one if asymmetric)
 keep_legacy_factor : if non-zero, multiply the final result by your old single length factor (e.g., 0.95) for compatibility.
 legacy_factor_value: that legacy factor (ignored if keep_legacy_factor=0)

 Returns:
 Physical length to cut before fine tuning.
 -------------------------------------------------------------------------- */
double interdigital_round_rod_length_corrected(double f0, double d, double h, double g, int keep_legacy_factor, double legacy_factor_value) {
    /* Quarter-wave in air (TEM-like), εeff≈1 */
    const double Lq = C0 / (4.0 * f0);

    const double dL_open = open_end_delta_l_round_rod(d, h);
    const double dL_gap = end_gap_delta_l_round_rod(d, g, h);

    /* Subtract the electrical lengthening from the geometric quarter-wave */
    double L = Lq - (dL_open + dL_gap);

    /* Safety clamps: don’t let us go crazy if inputs are odd */
    const double minL = 0.75 * Lq; /* generous lower bound */
    const double maxL = 1.05 * Lq; /* and upper bound */
    L = clamp(L, minL, maxL);

    if (keep_legacy_factor) {
        L *= legacy_factor_value; /* optional compatibility layer */
    }

    /* diagnostics */
    //printf("[rod-length] f0=%.6f MHz  Lq=%.6f mm  ΔL_open=%.3f mm  ΔL_gap=%.3f mm  -> L=%.6f mm%s\n", f0 * 1e-6, Lq * 1e3, dL_open * 1e3, dL_gap * 1e3, L * 1e3,
    //        keep_legacy_factor ? " (legacy factor applied)" : "");
    return L;
}

double rod_length(double f0, double h, double d, double length_factor) {
    /* If you don’t know g yet at this call site, pass a representative
     design gap (e.g., the smallest end gap in the passband section). */
    const double assumed_g = 0.5 * d; /* conservative default; override if you can */

    int keep = 1;
    double L = interdigital_round_rod_length_corrected(f0, d, h, assumed_g, keep, length_factor);

    return L;
}

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

