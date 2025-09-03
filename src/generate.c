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
#include <math.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>

#include "calculations.h"
#include "generate.h"

const char *material_str[] = { "PEC", "COPPER", "BRASS", "ALUMINUM" };

void generate_openEMS_script(const char *filename, double f0_MHz, double BW_MHz, double R_ohm, double H_mm, double D_mm, double E_mm, int ele, double ripple_dB,
        double *pos, double *gap, double *rod_lengths, double box_height, double box_length, double tap_z1, double tap_zN, filter_material_t material) {

}

void generate_dxf(const char *filename, double f0_MHz, double BW_MHz, double R_ohm, double H_mm, double D_mm, double E_mm, int ele, double ripple_dB,
        double *pos, double *gap, double *rod_lengths, double box_height, double box_length, double tap_z1, double tap_zN) {
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        fprintf(stderr, "Error opening %s: %s\n", filename, strerror(errno));
        exit(1);
    }

    double box_width_front = H_mm;    // Set to H_mm to represent depth accurately
    double y_center = 0.0;            // Rod centers at y=0 for front view
    double y_min = -box_width_front / 2.0;
    double y_max = box_width_front / 2.0;

    // DXF Header
    fprintf(fp, "  0\nSECTION\n  2\nHEADER\n  9\n$ACADVER\n  1\nAC1006\n  0\nENDSEC\n");
    fprintf(fp, "  0\nSECTION\n  2\nENTITIES\n");

    // *** Front View (Top-Down): Box outline, rods as circles, dimensions for E, gaps ***
    // Box outline (rectangle)
    fprintf(fp, "  0\nLINE\n  8\n0\n 10\n0.0\n 20\n%.2f\n 30\n0.0\n 11\n%.2f\n 21\n%.2f\n 31\n0.0\n", y_min, box_length, y_min);  // Bottom
    fprintf(fp, "  0\nLINE\n  8\n0\n 10\n0.0\n 20\n%.2f\n 30\n0.0\n 11\n%.2f\n 21\n%.2f\n 31\n0.0\n", y_max, box_length, y_max);  // Top
    fprintf(fp, "  0\nLINE\n  8\n0\n 10\n0.0\n 20\n%.2f\n 30\n0.0\n 11\n0.0\n 21\n%.2f\n 31\n0.0\n", y_min, y_max);  // Left
    fprintf(fp, "  0\nLINE\n  8\n0\n 10\n%.2f\n 20\n%.2f\n 30\n0.0\n 11\n%.2f\n 21\n%.2f\n 31\n0.0\n", box_length, y_min, box_length, y_max);  // Right

    // Left wall text
    fprintf(fp, "  0\nTEXT\n  8\nLABELS\n 10\n-20.0\n 20\n0.0\n 30\n0.0\n 40\n2.0\n  1\nLeft Wall\n");

    // Rod circles (diameter D_mm)
    for (int i = 1; i <= ele; i++) {
        fprintf(fp, "  0\nCIRCLE\n  8\n0\n 10\n%.2f\n 20\n0.0\n 30\n0.0\n 40\n%.2f\n", pos[i], D_mm / 2.0);
        // Label rod
        fprintf(fp, "  0\nTEXT\n  8\nLABELS\n 10\n%.2f\n 20\n5.0\n 30\n0.0\n 40\n2.0\n  1\nRod %d\n", pos[i], i);
    }

    // Add D dimension to first rod in Front View (horizontal dimension, aligned with rod1 center, label to the right)
    fprintf(fp, "  0\nDIMENSION\n  8\nDIMENSIONS\n  2\n*\n 70\n32\n 71\n5\n");
    fprintf(fp, " 10\n%.2f\n 20\n%.2f\n 30\n0.0\n", pos[1] + D_mm / 2.0 + 10.0, y_center);
    fprintf(fp, " 13\n%.2f\n 23\n%.2f\n 33\n0.0\n", pos[1] - D_mm / 2.0, y_center);
    fprintf(fp, " 14\n%.2f\n 24\n%.2f\n 34\n0.0\n", pos[1] + D_mm / 2.0, y_center);
    fprintf(fp, "  1\nD=%.2f mm\n", D_mm);

    // Right wall text (adjusted to +10.0 for symmetry)
    fprintf(fp, "  0\nTEXT\n  8\nLABELS\n 10\n%.2f\n 20\n0.0\n 30\n0.0\n 40\n2.0\n  1\nRight Wall\n", box_length + 10.0);

    // Add H dimension indicating length of right line (vertical dimension near Right Wall)
    fprintf(fp, "  0\nDIMENSION\n  8\nDIMENSIONS\n  2\n*\n 70\n33\n 71\n5\n");
    fprintf(fp, " 10\n%.2f\n 20\n%.2f\n 30\n0.0\n", box_length + 5.0, (y_min + y_max) / 2.0);
    fprintf(fp, " 13\n%.2f\n 23\n%.2f\n 33\n0.0\n", box_length, y_min);
    fprintf(fp, " 14\n%.2f\n 24\n%.2f\n 34\n0.0\n", box_length, y_max);
    fprintf(fp, "  1\nH=%.2f mm\n", H_mm);

    // Dimensions: Left E
    fprintf(fp, "  0\nDIMENSION\n  8\nDIMENSIONS\n  2\n*\n 70\n32\n 71\n5\n");
    fprintf(fp, " 10\n%.2f\n 20\n%.2f\n 30\n0.0\n", (pos[1] + 0) / 2.0, y_min - 5.0);  // Adjusted to y_min - 5.0
    fprintf(fp, " 13\n0.0\n 23\n%.2f\n 33\n0.0\n", y_center);
    fprintf(fp, " 14\n%.2f\n 24\n%.2f\n 34\n0.0\n", pos[1], y_center);
    fprintf(fp, "  1\nE=%.2f mm\n", E_mm);

    // Gaps between rods
    for (int i = 0; i < ele - 1; i++) {  // gap[1..ele-1], assuming gap is 1-indexed
        fprintf(fp, "  0\nDIMENSION\n  8\nDIMENSIONS\n  2\n*\n 70\n32\n 71\n5\n");
        fprintf(fp, " 10\n%.2f\n 20\n%.2f\n 30\n0.0\n", (pos[i + 1] + pos[i + 2]) / 2.0, y_min - 5.0);
        fprintf(fp, " 13\n%.2f\n 23\n%.2f\n 33\n0.0\n", pos[i + 1], y_center);
        fprintf(fp, " 14\n%.2f\n 24\n%.2f\n 34\n0.0\n", pos[i + 2], y_center);
        fprintf(fp, "  1\nGap=%.2f mm\n", gap[i + 1]);
    }

    // Right E
    fprintf(fp, "  0\nDIMENSION\n  8\nDIMENSIONS\n  2\n*\n 70\n32\n 71\n5\n");
    fprintf(fp, " 10\n%.2f\n 20\n%.2f\n 30\n0.0\n", (pos[ele] + box_length) / 2.0, y_min - 5.0);
    fprintf(fp, " 13\n%.2f\n 23\n%.2f\n 33\n0.0\n", pos[ele], y_center);
    fprintf(fp, " 14\n%.2f\n 24\n%.2f\n 34\n0.0\n", box_length, y_center);
    fprintf(fp, "  1\nE=%.2f mm\n", E_mm);

    // Overall length
    fprintf(fp, "  0\nDIMENSION\n  8\nDIMENSIONS\n  2\n*\n 70\n32\n 71\n5\n");
    fprintf(fp, " 10\n%.2f\n 20\n%.2f\n 30\n0.0\n", box_length / 2.0, y_min - 10.0);  // Adjusted to y_min - 10.0
    fprintf(fp, " 13\n0.0\n 23\n%.2f\n 33\n0.0\n", y_center);
    fprintf(fp, " 14\n%.2f\n 24\n%.2f\n 34\n0.0\n", box_length, y_center);
    fprintf(fp, "  1\nLength=%.2f mm\n", box_length);

    // Label front view
    fprintf(fp, "  0\nTEXT\n  8\nLABELS\n 10\n%.2f\n 20\n%.2f\n 30\n0.0\n 40\n3.0\n  1\nFront View (Top-Down)\n", box_length / 2.0, y_max + 5.0);

    // *** Height View (Side): Aligned to the right of the Front View ***
    double offset_x_height = box_length + 50.0;
    double offset_y_height = -H_mm / 2.0;  // Align base (bottom) with Front View's bottom

    // Ground planes (horizontal lines for bottom/top)
    fprintf(fp, "  0\nLINE\n  8\n0\n 10\n%.2f\n 20\n%.2f\n 30\n0.0\n 11\n%.2f\n 21\n%.2f\n 31\n0.0\n", offset_x_height, offset_y_height + 0.0, offset_x_height + box_length, offset_y_height + 0.0);  // Bottom ground
    fprintf(fp, "  0\nLINE\n  8\n0\n 10\n%.2f\n 20\n%.2f\n 30\n0.0\n 11\n%.2f\n 21\n%.2f\n 31\n0.0\n", offset_x_height, offset_y_height + box_height,
            offset_x_height + box_length, offset_y_height + box_height);  // Top ground

    // Vertical lines to close the box (left and right)
    fprintf(fp, "  0\nLINE\n  8\n0\n 10\n%.2f\n 20\n%.2f\n 30\n0.0\n 11\n%.2f\n 21\n%.2f\n 31\n0.0\n", offset_x_height, offset_y_height + 0.0, offset_x_height, offset_y_height + box_height);  // Left vertical
    fprintf(fp, "  0\nLINE\n  8\n0\n 10\n%.2f\n 20\n%.2f\n 30\n0.0\n 11\n%.2f\n 21\n%.2f\n 31\n0.0\n", offset_x_height + box_length, offset_y_height + 0.0, offset_x_height + box_length, offset_y_height + box_height);  // Right vertical

    // Rods as rectangles with width D_mm
    double r = D_mm / 2.0;
    for (int i = 1; i <= ele; i++) {
        double x = offset_x_height + pos[i];
        double rod_len = rod_lengths[i - 1];
        double z_start = ((i % 2) == 1) ? 0.0 : (box_height - rod_len);  // Alternating: odd short bottom, even short top
        double z_end = z_start + rod_len;
        // Left vertical
        fprintf(fp, "  0\nLINE\n  8\n0\n 10\n%.2f\n 20\n%.2f\n 30\n0.0\n 11\n%.2f\n 21\n%.2f\n 31\n0.0\n", x - r, offset_y_height + z_start, x - r, offset_y_height + z_end);
        // Right vertical
        fprintf(fp, "  0\nLINE\n  8\n0\n 10\n%.2f\n 20\n%.2f\n 30\n0.0\n 11\n%.2f\n 21\n%.2f\n 31\n0.0\n", x + r, offset_y_height + z_start, x + r, offset_y_height + z_end);
        // Bottom horizontal
        fprintf(fp, "  0\nLINE\n  8\n0\n 10\n%.2f\n 20\n%.2f\n 30\n0.0\n 11\n%.2f\n 21\n%.2f\n 31\n0.0\n", x - r, offset_y_height + z_start, x + r, offset_y_height + z_start);
        // Top horizontal
        fprintf(fp, "  0\nLINE\n  8\n0\n 10\n%.2f\n 20\n%.2f\n 30\n0.0\n 11\n%.2f\n 21\n%.2f\n 31\n0.0\n", x - r, offset_y_height + z_end, x + r, offset_y_height + z_end);
        // Rod length label
        fprintf(fp, "  0\nTEXT\n  8\nLABELS\n 10\n%.2f\n 20\n%.2f\n 30\n0.0\n 40\n2.0\n  1\nL%d=%.2f\n", x + r + 1.0, offset_y_height + (z_start + z_end) / 2.0, i, rod_len);
    }

    // Dimensions: Box height (vertical)
    fprintf(fp, "  0\nDIMENSION\n  8\nDIMENSIONS\n  2\n*\n 70\n33\n 71\n5\n");
    fprintf(fp, " 10\n%.2f\n 20\n%.2f\n 30\n0.0\n", offset_x_height - 5.0, offset_y_height + box_height / 2.0);
    fprintf(fp, " 13\n%.2f\n 23\n%.2f\n 33\n0.0\n", offset_x_height, offset_y_height + 0.0);
    fprintf(fp, "  14\n%.2f\n 24\n%.2f\n 34\n0.0\n", offset_x_height, offset_y_height + box_height);
    fprintf(fp, "  1\nHeight=%.2f mm\n", box_height);

    // Tap positions (first and last rods)
    double x_first = offset_x_height + pos[1];
    double z_tap1_start = ((1 % 2) == 1) ? 0.0 : (box_height - tap_z1);  // Adjust based on short
    fprintf(fp, "  0\nDIMENSION\n  8\nDIMENSIONS\n  2\n*\n 70\n33\n 71\n5\n");
    fprintf(fp, " 10\n%.2f\n 20\n%.2f\n 30\n0.0\n", x_first - r - 5.0, offset_y_height + (z_tap1_start + tap_z1) / 2.0);
    fprintf(fp, "  13\n%.2f\n 23\n%.2f\n 33\n0.0\n", x_first, offset_y_height + z_tap1_start);
    fprintf(fp, "  14\n%.2f\n 24\n%.2f\n 34\n0.0\n", x_first, offset_y_height + tap_z1 + z_tap1_start);
    fprintf(fp, "  1\nTap1=%.2f mm\n", tap_z1);

    double x_last = offset_x_height + pos[ele];
    double z_tapN_start = ((ele % 2) == 1) ? 0.0 : (box_height - tap_zN);
    fprintf(fp, "  0\nDIMENSION\n  8\nDIMENSIONS\n  2\n*\n 70\n33\n 71\n5\n");
    fprintf(fp, " 10\n%.2f\n 20\n%.2f\n 30\n0.0\n", x_last + r + 5.0, offset_y_height + (z_tapN_start + tap_zN) / 2.0);
    fprintf(fp, "  13\n%.2f\n 23\n%.2f\n 33\n0.0\n", x_last, offset_y_height + z_tapN_start);
    fprintf(fp, "  14\n%.2f\n 24\n%.2f\n 34\n0.0\n", x_last, offset_y_height + tap_zN + z_tapN_start);
    fprintf(fp, "  1\nTapN=%.2f mm\n", tap_zN);

    // Label height view
    fprintf(fp, "  0\nTEXT\n  8\nLABELS\n 10\n%.2f\n 20\n%.2f\n 30\n0.0\n 40\n3.0\n  1\nHeight View (Side)\n", offset_x_height + box_length / 2.0,
            offset_y_height + box_height + 5.0);

    // Title and parameters
    fprintf(fp,
            "  0\nTEXT\n  8\nLABELS\n 10\n0.0\n 20\n%.2f\n 30\n0.0\n 40\n3.0\n  1\nInterdigital BPF: f0=%.2f MHz, BW=%.2f MHz, Poles=%d, Ripple=%.2f dB, R=%.2f Ohm\n",
            y_max + 10.0, f0_MHz, BW_MHz, ele, ripple_dB, R_ohm);

    // End DXF
    fprintf(fp, "  0\nENDSEC\n  0\nEOF\n");
    fclose(fp);
    printf("DXF file written to %s\n", filename);
}
