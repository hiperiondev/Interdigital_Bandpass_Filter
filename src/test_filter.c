/* interdigital_bpf.c
 Interdigital Bandpass Filter quick synthesizer (air-filled, round-rod cavity)
 - Quarter-wave resonators, alternate-short (interdigital)
 - k, Qe from low-pass prototype g-values (Butterworth or Chebyshev)
 - Coupling gaps from GM3SEK empirical formula
 - Tap from GM3SEK tapped-line formula
 - Rough Qu estimate and midband IL budget (IL_dB ≈ 4.343 * (QL/Qu) * sum(gk))
 - Simple narrowband insertion-loss template (Butterworth/Chebyshev mapping)

 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "coupling_matrix_sparams.h"
#include "rod.h"
#include "generate_openEMS.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define LINT_FACTOR   0.889   /* interior rod length factor relative to lambda/4 */
#define LEND_FACTOR   1.005   /* end rod length multiplier relative to interior */
#define K_QU_DEFAULT  0.06    /* empirical Qu constant (conservative default) */

static double inch(double mm) {
    return mm * 0.03937007874015748;
}
static double clamp(double x, double lo, double hi) {
    if (x < lo)
        return lo;
    if (x > hi)
        return hi;
    return x;
}

/* Butterworth prototype g-values (ripple = 0). g0=g_{n+1}=1. */
static void butterworth_g(int n, double *g) {
    g[0] = 1.0;
    for (int k = 1; k <= n; k++) {
        g[k] = 2.0 * sin((2.0 * k - 1.0) * M_PI / (2.0 * n));
    }
    g[n + 1] = 1.0;
}

/* Chebyshev prototype g-values (equal ripple r dB).
 Reference: standard closed-form from many texts/Wikipedia. */
static void chebyshev_g(int n, double ripple_dB, double *g) {
    /* Compute epsilon and beta (standard definition) */
    double eps = sqrt(pow(10.0, ripple_dB / 10.0) - 1.0);
    double beta = asinh(1.0 / eps) / n;

    /* Base term */
    g[0] = 1.0;

    /* Recurrence arrays (sized to safely handle n <= 30) */
    double sinhb = sinh(beta);
    double a[64], b[64];

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

int main(void) {
    /* Inputs */
    double f0_MHz, BW_MHz, R_ohm, H_mm, D_mm, E_mm;
    int NFR, ele;
    double ripple_dB;

    printf("Center frequency f0 (MHz): ");
    if (scanf("%lf", &f0_MHz) != 1)
        return 1;
    printf("Bandwidth BW (MHz): ");
    if (scanf("%lf", &BW_MHz) != 1)
        return 1;
    printf("System impedance R (Ohm): ");
    if (scanf("%lf", &R_ohm) != 1)
        return 1;
    printf("Ground-plane height H (mm) (recommended lambda/4: %.3f) :", (double) (299792458.0 / (4.0 * (f0_MHz * 1e6))) * 1000.0);
    if (scanf("%lf", &H_mm) != 1)
        return 1;

    double E_target_mm = 0.7 * H_mm;

    printf("Rod diameter D (mm)  (recommended H/3:  %.3f) :", (double) (H_mm / 3));
    if (scanf("%lf", &D_mm) != 1)
        return 1;
    printf("End spacing E (mm) (GM3SEK recommended:  %.3f mm): ", E_target_mm);
    if (scanf("%lf", &E_mm) != 1)
        return 1;
    printf("Number of freq points NFR: ");
    if (scanf("%d", &NFR) != 1)
        return 1;
    printf("Number of elements ele [2..19]: ");
    if (scanf("%d", &ele) != 1)
        return 1;
    printf("Passband ripple rip (dB) [0=Butterworth]: ");
    if (scanf("%lf", &ripple_dB) != 1)
        return 1;

    ele = (int) clamp(ele, 2, 19);
    NFR = (int) clamp(NFR, 11, 20001); /* keep reasonable */

    /* Enforce or default end-wall spacing E: GM3SEK requires end-resonator centers
     to be around 0.7 * H to realize the assumed impedance and coupling.
     If user provides E <= 0 we set a conservative default of 0.7*H and
     warn if the supplied E deviates strongly from 0.7*H. */
    if (E_mm <= 0.0) {
        E_mm = E_target_mm;
        printf("Note: End spacing E not provided or non-positive. Using default E = 0.7*H = %.3f mm\n", E_mm);
    }
    /* Warn if E differs from recommended 0.7*H by more than 10% of H */
    if (fabs(E_mm - E_target_mm) > 0.10 * H_mm) {
        double diff_pct = 100.0 * (E_mm - E_target_mm) / E_target_mm;
        printf("WARNING: End spacing E = %.3f mm differs from GM3SEK recommended 0.7*H = %.3f mm (%.1f%%).\n", E_mm, E_target_mm, diff_pct);
        printf("         This may invalidate simple spacing/impedance assumptions; consider using E ≈ 0.7*H or retune lengths.\n");
    }

    const double f0 = f0_MHz * 1e6;
    const double BW = BW_MHz * 1e6;
    const double FBW = BW / f0;
    const double QL = f0 / BW;

    double g[32] = { 0 }; /* g0..g_{n+1} */
    if (ripple_dB <= 0.0) {
        butterworth_g(ele, g);
    } else {
        chebyshev_g(ele, ripple_dB, g);
    }

    /* Couplings and external Q */
    double *k = (double*) calloc(ele, sizeof(double)); /* k[1..ele-1] */
    for (int i = 1; i <= ele - 1; i++) {
        k[i] = FBW / sqrt(g[i] * g[i + 1]);
    }
    double Qe1 = (g[0] * g[1]) / FBW;
    double QeN = (g[ele] * g[ele + 1]) / FBW;

    /* Spacing (edge-edge) from k, then center-center = gap + D */
    double *gap = (double*) calloc(ele, sizeof(double)); /* c[i]=gap between i and i+1, edge-edge */
    double *cc = (double*) calloc(ele, sizeof(double)); /* center-center distances */
    double d_over_h = D_mm / H_mm;
    for (int i = 1; i <= ele - 1; i++) {
        double Ki = k[i];
        if (Ki <= 0)
            Ki = 1e-6;
        double c_over_h = (log10(1.0 / Ki) + 0.91 * (d_over_h) - 0.048) / 1.37; /* GM3SEK */
        gap[i] = c_over_h * H_mm; /* edge-edge */
        if (gap[i] < 0.0)
            gap[i] = 0.0;
        cc[i] = gap[i] + D_mm; /* center-center */
    }

    /* Box coordinates along length: wall(0), center1=E, center2=E+cc1, ... , right wall = last center + E */
    double *pos = (double*) calloc(ele + 2, sizeof(double)); /* position of: 0: left wall, 1..ele: centers, ele+1: right wall */
    pos[0] = 0.0;
    pos[1] = E_mm;
    for (int i = 2; i <= ele; i++) {
        pos[i] = pos[i - 1] + cc[i - 1];
    }
    pos[ele + 1] = pos[ele] + E_mm;
    double box_length = pos[ele + 1];

    /* Quarter-wave height (air) */
    const double c0 = 299792458.0; /* m/s */
    double lambda_over4_mm = (c0 / (4.0 * f0)) * 1000.0;

    /* empirical mechanical length factor vs d/h:
     factor = 0.97 - 0.0666666667 * (d/h)  (clamped to [0.90,0.96])
     This yields ≈0.95 for typical d/h and drifts towards 0.90 for large rods or 0.96 for very small rods.
     (See GM3SEK: recommended ~0.95·λ/4 at UHF; treat this as a starting point.) */
    double Lcorr_factor = 0.97 - 0.0666666667 * d_over_h;
    if (Lcorr_factor > 0.96)
        Lcorr_factor = 0.96;
    if (Lcorr_factor < 0.90)
        Lcorr_factor = 0.90;
    double L_int_mm = Lcorr_factor * lambda_over4_mm;
    /* End resonators are usually slightly longer due to end-wall loading and tap screws.
     Use a modest default extension (+1%) relative to interior resonators; tune as needed. */
    double L_end_mm = L_int_mm * 1.01;
    printf("Note: Initial mechanical lengths: L_int factor=%.4f -> L_int=%.3f mm, L_end=%.3f mm (d/h=%.3f)\n", Lcorr_factor, L_int_mm, L_end_mm, d_over_h);

    /* Characteristic impedance of end resonator (trough-line, round rod) */
    double Z_end = 138.0 * log10((1.25 * H_mm) / D_mm);
    if (Z_end < 1.0)
        Z_end = 1.0; /* sanity */
    double Z_inner = Z_end + 5.0; /* simple uplift; typically a few ohms higher */

    /* Tap distance from shorted end using GM3SEK formula:
     t = (L/90deg) * asin( sqrt( (pi*R)/(4*Z*Qe1) ) )
     Implement with radians: t = L * (2/pi) * asin( sqrt( (M_PI*R)/(4*Z*Qe1) ) )
     Use end-resonator Z for the source. */
    double arg = (M_PI * R_ohm) / (4.0 * Z_end * Qe1);
    arg = sqrt(fabs(arg));
    arg = clamp(arg, 0.0, 0.999999);
    double tap_mm = L_end_mm * (2.0 / M_PI) * asin(arg);

    /* Estimated unloaded Q (very rough; tweak coefficient to match your workshop).
     Qu ≈ K * f_MHz * H_mm (empirical). */
    const double K_qu = 0.06; /* conservative air-cavity value; polish as you gain data */
    double Qu_est = K_qu * f0_MHz * H_mm;
    if (Qu_est < 100.0)
        Qu_est = 100.0;

    /* Midband IL estimate */
    double sum_g = 0.0;
    for (int i = 1; i <= ele; i++)
        sum_g += g[i];
    double IL_dB = 4.343 * (QL / Qu_est) * sum_g;

    /* Group delay (approx, midband, narrowband) ~ sum(g)/ (2π·BW) */
    double tau_s = sum_g / (2.0 * M_PI * BW);
    double tau_ns = tau_s * 1e9;

    printf("\n-------------------------------------------------------------------------\n");
    printf("Design data for a %d section interdigital bandpass filter.\n", ele);
    printf("Center Frequency           : %10.4f MHz\n", f0_MHz);
    printf("Passband Ripple            : %10.3f dB\n", ripple_dB <= 0.0 ? 0.0 : ripple_dB);
    printf("System Impedance           : %10.2f Ohm\n", R_ohm);
    double f1 = f0 - BW / 2.0, f2 = f0 + BW / 2.0;
    printf("Cutoff Frequency           : %10.4f MHz and %10.4f MHz\n", f1 * 1e-6, f2 * 1e-6);
    printf("Bandwidth (3dB)            : %10.4f MHz\n", BW_MHz);
    printf("Fractional Bandwidth (3dB) : %10.6f\n", FBW);
    printf("Filter Q                   : %10.3f\n", QL);
    printf("Estimated Qu               : %10.2f\n", Qu_est);
    printf("Loss, based on this Qu     : %10.3f dB\n", IL_dB);
    printf("Passband Delay             : %10.3f ns\n", tau_ns);

    printf("-------------------------------------------------------------------------\n");
    printf("Quarter Wavelength           : %6.2f mm\n", lambda_over4_mm);
    printf("Length interior Element      : %6.2f mm\n", L_int_mm);
    printf("Length of end Element        : %6.2f mm\n", L_end_mm);
    printf("Ground plane space           : %6.2f mm\n", H_mm);
    printf("Rod Diameter                 : %6.2f mm\n", D_mm);
    printf("End plate to center of Rod   : %6.2f mm\n", E_mm);
    printf("Tap to shorted End           : %6.2f mm\n", tap_mm);
    printf("Impedance end Rod            : %6.2f Ohm\n", Z_end);
    printf("Impedance inner Rod          : %6.2f Ohm\n", Z_inner);
    printf("Impedance ext. line          : %6.2f Ohm\n", R_ohm);

    printf("-------------------------------------------------------------------------\n");
    printf("                     **** Dimensions, mm (inch) ****\n");
    printf(" #    End to Center      Center-Center        G[k]         Q/Coup\n");
    printf(" 0  %8.2f (%6.3f)\n", 0.0, 0.0);
    for (int i = 1; i <= ele; i++) {
        double ec_mm = pos[i];
        double cc_mm = (i < ele) ? cc[i] : 0.0;
        printf("%2d  %8.2f (%6.3f)   ", i, ec_mm, inch(ec_mm));
        if (i < ele)
            printf("%8.2f (%6.3f)   ", cc_mm, inch(cc_mm));
        else
            printf("    %6s %8s   ", "", "");
        printf("%8.3f   ", g[i]);
        if (i < ele)
            printf("%10.6f", k[i]);
        else
            printf("%10.3f", QeN);
        printf("\n");
    }
    printf("%2d  %8.2f (%6.3f)\n", ele + 1, pos[ele + 1], inch(pos[ele + 1]));

    printf("-------------------------------------------------------------------------\n");
    printf("\n=== Rod length corrections (geometry-aware) ===\n");
    double f0_Hz = f0; /* f0 is in Hz in this scope */
    double h_m = H_mm / 1000.0;
    double d_m = D_mm / 1000.0;
    double *rod_lengths = (double*) calloc(ele, sizeof(double));
    for (int i = 1; i <= ele; i++) {
        /* edge-edge gaps (mm):
         - interior: gap[i-1] is between rod (i-1) and rod i
         - for rod i, left neighbor gap is gap[i-1] (if i>1), right neighbor is gap[i] (if i<ele)
         - end gaps to the cavity wall (edge-edge) are approximated as E_mm - D_mm/2
         */
        double g_left_mm = (i == 1) ? (pos[1] - pos[0] - (D_mm / 2.0)) : gap[i - 1];
        double g_right_mm = (i == ele) ? (pos[ele + 1] - pos[ele] - (D_mm / 2.0)) : gap[i];

        /* choose the larger gap as the 'open' end gap heuristic (alternating short/open not explicit here) */
        double g_end_mm = (g_right_mm > g_left_mm) ? g_right_mm : g_left_mm;
        if (g_end_mm < 0.0)
            g_end_mm = 0.0;

        double g_m = g_end_mm / 1000.0;
        int keep_legacy = 0;
        double legacy_val = (i == 1 || i == ele) ? LEND_FACTOR : LINT_FACTOR;

        /* call geometry-aware corrected length (returns meters) */
        double L_corr_m = interdigital_round_rod_length_corrected(f0_Hz, d_m, h_m, g_m, keep_legacy, legacy_val);
        double L_corr_mm = L_corr_m * 1000.0;
        rod_lengths[i - 1] = L_corr_mm;

        /* call legacy-style helper (returns meters) and convert to mm */
        double L_legacy_m = rod_length(f0_Hz, h_m, d_m, legacy_val);
        double L_legacy_mm = L_legacy_m * 1000.0;

        printf("Rod %2d: g_left=%.3f mm g_right=%.3f mm | corrected L = %.3f mm | legacy L = %.3f mm\n", i, g_left_mm, g_right_mm, L_corr_mm, L_legacy_mm);
    }

    printf("-------------------------------------------------------------------------\n");
    printf("                      **** Box inside dimensions ****\n");
    printf("Height : %5.2f mm or %5.3f inch\n", lambda_over4_mm, inch(lambda_over4_mm));
    printf("Length : %5.2f mm or %5.3f inch\n", box_length, inch(box_length));
    printf("Depth  : %5.2f mm or %5.3f inch\n", H_mm, inch(H_mm));
    printf("-------------------------------------------------------------------------\n");

    generate_openEMS_script("bpf_sim.m", f0_MHz, BW_MHz, R_ohm, H_mm, D_mm, E_mm, ele, ripple_dB, pos, gap, rod_lengths);
    free(rod_lengths);

    /*
     double fstart = f0 * (1.0 - 2.5 * FBW);
     double fstop = f0 * (1.0 + 2.5 * FBW);
     if (fstart < 1.0)
     fstart = 1.0;
     printf("\n");
     printf("Insertion-loss template (narrowband %s) over [%.3f..%.3f MHz]\n", (ripple_dB <= 0.0 ? "Butterworth" : "Chebyshev"), fstart * 1e-6, fstop * 1e-6);
     printf("\n");
     printf("   f (MHz)    |S21| (dB)\n");
     for (int i = 0; i < NFR; i++) {
     double frac = (double) i / (double) (NFR - 1);
     double f = fstart + frac * (fstop - fstart);
     double omega = ((f / f0) - (f0 / f)) / FBW; // low-pass mapping
     double H2;
     if (ripple_dB <= 0.0) {
     H2 = 1.0 / (1.0 + pow(fabs(omega), 2.0 * ele));
     } else {
     double eps = sqrt(pow(10.0, ripple_dB / 10.0) - 1.0);
     double T = Tn(ele, omega);
     H2 = 1.0 / (1.0 + (eps * eps * T * T));
     }
     double S21_dB = -10.0 * log10(H2);
     printf("%10.4f   %8.3f\n", f * 1e-6, S21_dB);
     }
     */
    /*
     // --- Compute S-parameters from coupling matrix using computed k and Q estimates --- integration with coupling_matrix_sparams
     int N = ele; // number of resonators
     // Build center frequencies array spanning f0 +/- 2*BW (sample)
     double *freqs = (double*) malloc(NFR * sizeof(double));
     if (freqs) {
     for (int i = 0; i < NFR; i++) {
     freqs[i] = f0 - 2.0 * BW + i * (4.0 * BW / (double) (NFR - 1));
     }
     }
     // Build coupling matrix M (N x N)
     double **M = (double**) malloc(N * sizeof(double*));
     if (M) {
     for (int i = 0; i < N; i++) {
     M[i] = (double*) calloc(N, sizeof(double));
     }
     // Note: k[] in this program is indexed from 1..ele-1, while build_coupling_matrix_from_k expects 0..N-2
     double *k0 = k + 1; //align to 0-based expected by builder
     build_coupling_matrix_from_k(k0, N, M); // build matrix from k

     // Allocate S-parameter arrays and compute
     double complex *S11 = (double complex*) malloc(NFR * sizeof(double complex));
     double complex *S21 = (double complex*) malloc(NFR * sizeof(double complex));
     if (S11 && S21) {
     compute_sparams_from_M(N, M, f0, BW, Qe1, QeN, Qu_est, freqs, NFR, S11, S21); //call into coupling module

     // Print results (preserve original printf formats)
     for (int ii = 0; ii < NFR; ii++) {
     double magS21 = 20.0 * log10(cabs(S21[ii]));
     double magS11 = 20.0 * log10(cabs(S11[ii]));
     printf("%.0f Hz: |S21|=%.2f dB, |S11|=%.2f dB\n", freqs[ii], magS21, magS11);
     }
     }
     // free S-params arrays
     if (S11)
     free(S11);
     if (S21)
     free(S21);

     // free matrix
     for (int i = 0; i < N; i++)
     free(M[i]);
     free(M);
     }

     if (freqs)
     free(freqs);
     */

    free(k);
    free(gap);
    free(cc);
    free(pos);
    return 0;
}
