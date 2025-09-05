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
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        fprintf(stderr, "Error opening %s: %s\n", filename, strerror(errno));
        exit(1);
    }

    // Header and setup
    fprintf(fp, "%% Auto-generated openEMS simulation script for interdigital BPF\n");
    fprintf(fp, "%% Parameters: f0=%.2f MHz, BW=%.2f MHz, R=%.2f Ohm, H=%.2f mm, D=%.2f mm, E=%.2f mm, ele=%d, ripple=%.2f dB, material=%s\n\n", f0_MHz, BW_MHz,
            R_ohm, H_mm, D_mm, E_mm, ele, ripple_dB, material_str[material]);
    fprintf(fp, "clear; close all; clc;\n");
    fprintf(fp, "unit = 1e-3; %% mm\n\n");

    // Paths - note to user: ensure OPENEMS_PATH points to openEMS/matlab, CSXCAD_PATH to CSXCAD/matlab (swap if needed)
    fprintf(fp, "openEMS_path = getenv('OPENEMS_PATH');\n");
    fprintf(fp, "if isempty(openEMS_path)\n");
    fprintf(fp, "    error('OPENEMS_PATH environment variable not set. Please set it to the openEMS/matlab directory.');\n");
    fprintf(fp, "end\n");
    fprintf(fp, "addpath(openEMS_path);\n\n");
    fprintf(fp, "csxcad_path = getenv('CSXCAD_PATH');\n");
    fprintf(fp, "if isempty(csxcad_path)\n");
    fprintf(fp, "    error('CSXCAD_PATH environment variable not set. Please set it to the CSXCAD/matlab directory.');\n");
    fprintf(fp, "end\n");
    fprintf(fp, "addpath(csxcad_path);\n\n");
    fprintf(fp, "physical_constants;\n");
    fprintf(fp, "%% Debug: Print versions\n");
    fprintf(fp, "fprintf('CSXCAD path: %%s\\n', csxcad_path);\n");
    fprintf(fp, "try\n");
    fprintf(fp, "    ver = CSXCAD_version;\n");
    fprintf(fp, "    fprintf('CSXCAD version: %%s\\n', ver);\n");
    fprintf(fp, "catch\n");
    fprintf(fp, "    fprintf('CSXCAD version: unknown\\n');\n");
    fprintf(fp, "end\n\n");

    // Custom DefineRectGrid (kept for compatibility)
    fprintf(fp, "function CSX = DefineRectGrid(CSX, mesh)\n");
    fprintf(fp, "if ~isfield(mesh, 'lines') || length(mesh.lines) < 3\n");
    fprintf(fp, "    error('mesh.lines must have at least 3 components');\n");
    fprintf(fp, "end\n");
    fprintf(fp, "if any(cellfun(@isempty, mesh.lines))\n");
    fprintf(fp, "    error('mesh.lines is empty in one or more directions');\n");
    fprintf(fp, "end\n");
    fprintf(fp, "if any(cellfun(@(x) any(diff(x) <= 0), mesh.lines))\n");
    fprintf(fp, "    error('mesh.lines is not strictly increasing');\n");
    fprintf(fp, "end\n");
    fprintf(fp, "if any(cellfun(@length, mesh.lines) < 15)\n");
    fprintf(fp, "    error('each mesh.lines must have at least 15 points');\n");
    fprintf(fp, "end\n");
    fprintf(fp, "fprintf('DefineRectGrid: lines{1} first 5: %%s\\n', sprintf('%%.3e ', mesh.lines{1}(1:min(5,end))));\n");
    fprintf(fp, "fprintf('DefineRectGrid: lines{2} first 5: %%s\\n', sprintf('%%.3e ', mesh.lines{2}(1:min(5,end))));\n");
    fprintf(fp, "fprintf('DefineRectGrid: lines{3} first 5: %%s\\n', sprintf('%%.3e ', mesh.lines{3}(1:min(5,end))));\n");
    fprintf(fp, "CSX.RectilinearGrid.XLines = mesh.lines{1};\n");
    fprintf(fp, "CSX.RectilinearGrid.YLines = mesh.lines{2};\n");
    fprintf(fp, "CSX.RectilinearGrid.ZLines = mesh.lines{3};\n");
    fprintf(fp, "CSX.RectilinearGrid.Type = 'Rectilinear';\n");
    fprintf(fp, "CSX.RectilinearGrid.DeltaUnit = 1; %% meters\n");
    fprintf(fp, "CSX.RectilinearGrid.CoordSystem = 0;\n");
    fprintf(fp, "CSX.RectilinearGrid.GridType = 'Cartesian';\n");
    fprintf(fp, "CSX.RectilinearGrid.Attribs.ID = 1;\n");
    fprintf(fp, "fprintf('Cartesian grid set with %%d x %%d x %%d points\\n', ...\n");
    fprintf(fp, "    length(mesh.lines{1}), length(mesh.lines{2}), length(mesh.lines{3}));\n");
    fprintf(fp, "end\n\n");

    // Parameters
    fprintf(fp, "f0 = %g * 1e6;\n", f0_MHz);
    fprintf(fp, "bw = %g * 1e6;\n", BW_MHz);
    fprintf(fp, "R = %g;\n", R_ohm);
    fprintf(fp, "h_mm = %g;\n", box_height);
    fprintf(fp, "w_mm = %g;\n", H_mm);
    fprintf(fp, "d_mm = %g;\n", D_mm);
    fprintf(fp, "e_mm = %g;\n", E_mm);
    fprintf(fp, "l_mm = %g;\n", box_length);
    fprintf(fp, "tap1_mm = %g;\n", tap_z1);
    fprintf(fp, "tapN_mm = %g;\n", tap_zN);
    fprintf(fp, "n_ele = %d;\n", ele);
    fprintf(fp, "pos_mm = [");
    for (int i = 0; i <= ele + 1; i++) {
        fprintf(fp, " %g", pos[i]);
    }
    fprintf(fp, "];\n");
    fprintf(fp, "rod_len_mm = [");
    for (int i = 0; i < ele; i++) {
        fprintf(fp, " %g", rod_lengths[i]);
    }
    fprintf(fp, "];\n\n");

    fprintf(fp, "h = h_mm * unit;\n");
    fprintf(fp, "w = w_mm * unit;\n");
    fprintf(fp, "l = l_mm * unit;\n");
    fprintf(fp, "d = d_mm * unit;\n");
    fprintf(fp, "r = d / 2;\n");
    fprintf(fp, "pos = pos_mm * unit;\n");
    fprintf(fp, "rod_len = rod_len_mm * unit;\n");
    fprintf(fp, "tap1 = tap1_mm * unit;\n");
    fprintf(fp, "tapN = tapN_mm * unit;\n");
    fprintf(fp, "port_half = 1.5 * unit;\n");
    fprintf(fp, "mesh_offset = 0.333 * unit; %% Adjusted for 1/3 rule approximation\n\n");

    // FDTD setup (increased NrTS to satisfy >3x excitation)
    fprintf(fp, "FDTD = InitFDTD('NrTS', 10000000, 'EndCriteria', 1e-5);\n");
    fprintf(fp, "FDTD = SetGaussExcite(FDTD, f0, bw/2);\n");
    fprintf(fp, "BC = {'PEC' 'PEC' 'PEC' 'PEC' 'PEC' 'PEC'};\n");
    fprintf(fp, "FDTD = SetBoundaryCond(FDTD, BC);\n\n");

    // CSX setup
    fprintf(fp, "CSX = InitCSX();\n\n");

    // Material (corrected for PEC as metal, others as conducting with kappa)
    if (material == FM_PEC) {
        fprintf(fp, "CSX = AddMetal(CSX, '%s');\n\n", material_str[material]);
    } else {
        double kappa = 0.0;
        if (material == FM_COPPER)
            kappa = 5.8e7;
        else if (material == FM_BRASS)
            kappa = 1.5e7;
        else if (material == FM_ALUMINUM)
            kappa = 3.8e7;
        fprintf(fp, "CSX = AddMaterial(CSX, '%s', 'Kappa', %g);\n\n", material_str[material], kappa);
    }

    // Add rods (alternating shorts)
    for (int i = 1; i <= ele; i++) {
        double start_z = ((i % 2) == 1) ? 0.0 : (box_height - rod_lengths[i - 1]);
        double stop_z = start_z + rod_lengths[i - 1];
        fprintf(fp, "start = [%.4f, 0, %.4f] * unit;\n", pos[i], start_z);
        fprintf(fp, "stop = [%.4f, 0, %.4f] * unit;\n", pos[i], stop_z);
        fprintf(fp, "CSX = AddCylinder(CSX, '%s', 10, start, stop, r);\n\n", material_str[material]);
    }

    // Add ports (with explicit half-width)
    double port_half_mm = 1.5;
    // Port 1
    double start_z1 = ((1 % 2) == 1) ? tap_z1 - port_half_mm : (box_height - rod_lengths[0] + tap_z1 - port_half_mm);
    double stop_z1 = start_z1 + 2 * port_half_mm;
    fprintf(fp, "start = [%.4f, 0, %.5f] * unit;\n", pos[1], start_z1);
    fprintf(fp, "stop = [%.4f, 0, %.5f] * unit;\n", pos[1], stop_z1);
    fprintf(fp, "dir = [0 0 1];\n");
    fprintf(fp, "CSX = AddLumpedPort(CSX, 5, 1, R, start, stop, dir, 1);\n\n");

    // Port N
    double start_zN = ((ele % 2) == 1) ? tap_zN - port_half_mm : (box_height - rod_lengths[ele - 1] + tap_zN - port_half_mm);
    double stop_zN = start_zN + 2 * port_half_mm;
    fprintf(fp, "start = [%.4f, 0, %.5f] * unit;\n", pos[ele], start_zN);
    fprintf(fp, "stop = [%.4f, 0, %.5f] * unit;\n", pos[ele], stop_zN);
    fprintf(fp, "dir = [0 0 1];\n");
    fprintf(fp, "CSX = AddLumpedPort(CSX, 5, 2, R, start, stop, dir, 0);\n\n");

    // Mesh (corrected: increased min_cell to reduce small cells, adjusted offset for thirds rule)
    fprintf(fp, "mesh = [];\n");
    fprintf(fp, "x_points = unique(sort([pos - r, pos, pos + r, pos - r + mesh_offset, pos + r - mesh_offset]));\n");
    fprintf(fp, "mesh.lines{1} = x_points;\n");
    fprintf(fp, "mesh.lines{2} = linspace(-w/2, w/2, 15);\n");
    fprintf(fp, "z_points = [0, h, tap1 - port_half, tap1 + port_half, tapN - port_half, tapN + port_half, ...\n");
    fprintf(fp, "            h - tap1 - port_half, h - tap1 + port_half, h - tapN - port_half, h - tapN + port_half];\n");
    for (int i = 0; i < ele; i++) {
        fprintf(fp, "z_points = [z_points, rod_len(%d), h - rod_len(%d), rod_len(%d) + mesh_offset, h - rod_len(%d) - mesh_offset];\n", i + 1, i + 1, i + 1,
                i + 1);
    }
    fprintf(fp, "mesh.lines{3} = unique(sort(z_points));\n");
    fprintf(fp, "%% Smooth mesh to prevent tiny cells and small timesteps\n");
    fprintf(fp, "min_cell = unit * 2; %% Increased to 2 mm to avoid small timesteps\n");
    fprintf(fp, "max_cell = unit * 30; %% ~lambda/10 to avoid max_res warning\n");
    fprintf(fp, "mesh.lines{1} = SmoothMeshLines(mesh.lines{1}, min_cell, 1.5, max_cell); %% Reduced ratio for smoother transitions\n");
    fprintf(fp, "mesh.lines{2} = SmoothMeshLines(mesh.lines{2}, min_cell, 1.5, max_cell);\n");
    fprintf(fp, "mesh.lines{3} = SmoothMeshLines(mesh.lines{3}, min_cell, 1.5, max_cell);\n");
    fprintf(fp, "CSX = DefineRectGrid(CSX, mesh);\n\n");

    // Run simulation
    fprintf(fp, "Sim_Path = '.';\n");
    fprintf(fp, "if ~exist(Sim_Path, 'dir'); mkdir(Sim_Path); end\n");
    fprintf(fp, "WriteOpenEMS([Sim_Path '/geometry.xml'], FDTD, CSX);\n");
    fprintf(fp, "RunOpenEMS(Sim_Path, 'geometry.xml');\n");

    fclose(fp);
    printf("openEMS script written to %s\n", filename);
}

void generate_dxf(const char *filename, double f0_MHz, double BW_MHz, double R_ohm, double H_mm, double D_mm, double E_mm, int ele, double ripple_dB,
        double *pos, double *gap, double *rod_lengths, double box_height, double box_length, double tap_z1, double tap_zN) {
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        fprintf(stderr, "Error opening %s: %s\n", filename, strerror(errno));
        exit(1);
    }

    double box_width_front = H_mm;
    double y_center = 0.0;
    double y_min = -box_width_front / 2.0;
    double y_max = box_width_front / 2.0;

    // DXF Header
    fprintf(fp, "  0\nSECTION\n  2\nHEADER\n  9\n$ACADVER\n  1\nAC1006\n  0\nENDSEC\n");
    fprintf(fp, "  0\nSECTION\n  2\nENTITIES\n");

    // Front View (Top-Down): Box outline, rods as circles, dimensions for E, gaps ***
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
    fprintf(fp, "  0\nTEXT\n  8\nLABELS\n 10\n%.2f\n 20\n%.2f\n 30\n0.0\n 40\n3.0\n  1\nFront View (Top-Down)\n", box_length / 2.0 - 20.0, y_max + 5.0);

    // Height View (Side)
    double offset_x_height = box_length + 50.0;
    double offset_y_height = -H_mm / 2.0;  // Align base (bottom) with Front View's bottom

    // Ground planes (horizontal lines for bottom/top)
    fprintf(fp, "  0\nLINE\n  8\n0\n 10\n%.2f\n 20\n%.2f\n 30\n0.0\n 11\n%.2f\n 21\n%.2f\n 31\n0.0\n", offset_x_height, offset_y_height + 0.0,
            offset_x_height + box_length, offset_y_height + 0.0);  // Bottom ground
    fprintf(fp, "  0\nLINE\n  8\n0\n 10\n%.2f\n 20\n%.2f\n 30\n0.0\n 11\n%.2f\n 21\n%.2f\n 31\n0.0\n", offset_x_height, offset_y_height + box_height,
            offset_x_height + box_length, offset_y_height + box_height);  // Top ground

    // Vertical lines to close the box (left and right)
    fprintf(fp, "  0\nLINE\n  8\n0\n 10\n%.2f\n 20\n%.2f\n 30\n0.0\n 11\n%.2f\n 21\n%.2f\n 31\n0.0\n", offset_x_height, offset_y_height + 0.0, offset_x_height,
            offset_y_height + box_height);  // Left vertical
    fprintf(fp, "  0\nLINE\n  8\n0\n 10\n%.2f\n 20\n%.2f\n 30\n0.0\n 11\n%.2f\n 21\n%.2f\n 31\n0.0\n", offset_x_height + box_length, offset_y_height + 0.0,
            offset_x_height + box_length, offset_y_height + box_height);  // Right vertical

    // Rods as rectangles with width D_mm
    double r = D_mm / 2.0;
    for (int i = 1; i <= ele; i++) {
        double x = offset_x_height + pos[i];
        double rod_len = rod_lengths[i - 1];
        double z_start = ((i % 2) == 1) ? 0.0 : (box_height - rod_len);  // Alternating: odd short bottom, even short top
        double z_end = z_start + rod_len;
        // Left vertical
        fprintf(fp, "  0\nLINE\n  8\n0\n 10\n%.2f\n 20\n%.2f\n 30\n0.0\n 11\n%.2f\n 21\n%.2f\n 31\n0.0\n", x - r, offset_y_height + z_start, x - r,
                offset_y_height + z_end);
        // Right vertical
        fprintf(fp, "  0\nLINE\n  8\n0\n 10\n%.2f\n 20\n%.2f\n 30\n0.0\n 11\n%.2f\n 21\n%.2f\n 31\n0.0\n", x + r, offset_y_height + z_start, x + r,
                offset_y_height + z_end);
        // Bottom horizontal
        fprintf(fp, "  0\nLINE\n  8\n0\n 10\n%.2f\n 20\n%.2f\n 30\n0.0\n 11\n%.2f\n 21\n%.2f\n 31\n0.0\n", x - r, offset_y_height + z_start, x + r,
                offset_y_height + z_start);
        // Top horizontal
        fprintf(fp, "  0\nLINE\n  8\n0\n 10\n%.2f\n 20\n%.2f\n 30\n0.0\n 11\n%.2f\n 21\n%.2f\n 31\n0.0\n", x - r, offset_y_height + z_end, x + r,
                offset_y_height + z_end);
        // Rod length label
        fprintf(fp, "  0\nTEXT\n  8\nLABELS\n 10\n%.2f\n 20\n%.2f\n 30\n0.0\n 40\n2.0\n  1\nL%d=%.2f\n", x + r + 1.0, offset_y_height + (z_start + z_end) / 2.0,
                i, rod_len);
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
    fprintf(fp, "  0\nTEXT\n  8\nLABELS\n 10\n%.2f\n 20\n%.2f\n 30\n0.0\n 40\n3.0\n  1\nHeight View (Side)\n", offset_x_height + box_length / 2.0 - 17.0,
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
