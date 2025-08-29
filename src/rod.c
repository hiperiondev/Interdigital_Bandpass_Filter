/*
 * Copyright 2025 Emiliano Augusto Gonzalez (egonzalez . hiperion @ gmail . com))
 * * Project Site: https://github.com/hiperiondev/Interdigital_Bandpass_Filter *
 * * This code is based on: https://www.changpuak.ch/electronics/interdigital_bandpass_filter_designer.php
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

/* --------------------------------------------------------------------------
 Interdigital round-rod resonator length — geometry-aware first-pass
 - Air-filled, two lids (stripline-like), quarter-wave rods
 - Accounts for: d/h via equivalent strip width, lid spacing h, end-gap g
 - Keeps your legacy "length factor" available if you still want it.
 References in comments; see analysis note.
 -------------------------------------------------------------------------- */

#include <math.h>
#include <stdio.h>

/* Physical constants */
#ifndef C0
#define C0 299792458.0
#endif

static inline double clampd(double x, double lo, double hi) {
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
    const double wh = clampd(w / h, 1e-6, 1e6);
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
    L = clampd(L, minL, maxL);

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

