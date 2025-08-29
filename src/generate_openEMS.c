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

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>

#include "rod.h"

void generate_openEMS_script(const char *filename, double f0_MHz, double BW_MHz, double R_ohm, double H_mm, double D_mm, double E_mm, int ele, double ripple_dB,
        double *pos, double *gap, double *rod_lengths) {
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        fprintf(stderr, "Error opening %s: %s\n", filename, strerror(errno));
        exit(1);
    }

    double f0 = f0_MHz * 1e6;  // Hz
    double lambda4_mm = 299792458.0 / (4 * f0) * 1000.0;  // Precise quarter-wavelength in mm
    double box_length = pos[ele + 1];  // Cavity length from positions
    double f_min = f0 - 2 * BW_MHz * 1e6;
    double f_max = f0 + 2 * BW_MHz * 1e6;

    // Validate geometry
    if (box_length <= 0 || H_mm <= 0 || lambda4_mm <= 0) {
        fprintf(stderr, "Error: Invalid box dimensions (length=%.2f, ground spacing=%.2f, height=%.2f)\n", box_length, H_mm, lambda4_mm);
        fclose(fp);
        exit(1);
    }

    // Header and setup
    fprintf(fp, "%% openEMS script for %d-pole interdigital BPF at %.2f MHz\n", ele, f0_MHz);
    fprintf(fp, "close all; clear; clc;\n");
    fprintf(fp, "addpath(getenv('OPENEMS_PATH')); %% Set OPENEMS_PATH to openEMS/matlab directory\n");
    fprintf(fp, "addpath(getenv('CSXCAD_PATH')); %% Set CSXCAD_PATH to CSXCAD/matlab directory\n");
    fprintf(fp, "if ~exist('InitFDTD', 'file') || ~exist('DefineRectGrid', 'file')\n");
    fprintf(fp, "    error('openEMS or CSXCAD not found. Check OPENEMS_PATH and CSXCAD_PATH.');\n");
    fprintf(fp, "end\n");
    fprintf(fp, "physical_constants;\n");
    fprintf(fp, "unit = 1e-3; %% mm\n");
    fprintf(fp, "f0 = %.2e;\n", f0);
    fprintf(fp, "f_min = %.2e; f_max = %.2e;\n", f_min, f_max);
    fprintf(fp, "box_l = %.2f; %% box length in x\n", box_length);
    fprintf(fp, "ground_spacing = %.2f; %% ground spacing in y\n", H_mm);
    fprintf(fp, "box_height = %.2f; %% box height in z (rod direction)\n", lambda4_mm);
    fprintf(fp, "Sim_Path = 'sim_bpf';\n");
    fprintf(fp, "mkdir(Sim_Path);\n");
    fprintf(fp, "disp(['box_l = ', num2str(box_l), ' mm']);\n");
    fprintf(fp, "disp(['ground_spacing = ', num2str(ground_spacing), ' mm']);\n");
    fprintf(fp, "disp(['box_height = ', num2str(box_height), ' mm']);\n");
    fprintf(fp, "FDTD = InitFDTD('NrTS', 30000000, 'EndCriteria', 1e-5);\n");
    fprintf(fp, "FDTD = SetGaussExcite(FDTD, f0, (f_max - f_min)/2);\n");
    fprintf(fp, "BC = {'PEC' 'PEC' 'PEC' 'PEC' 'PEC' 'PEC'}; %% All PEC for closed metal cavity\n");
    fprintf(fp, "FDTD = SetBoundaryCond(FDTD, BC);\n");

    // CSX and mesh
    fprintf(fp, "CSX = InitCSX();\n");
    fprintf(fp, "CSX = AddMetal(CSX, 'PEC'); %% Define Perfect Electric Conductor (PEC)\n");
    fprintf(fp, "resolution = c0 / (f_max) / unit / 20; %% lambda/20 for finer mesh\n");
    fprintf(fp, "mesh.x = linspace(-box_l/2, box_l/2, max(3, ceil(box_l / resolution)));\n");
    fprintf(fp, "mesh.y = linspace(-ground_spacing/2, ground_spacing/2, max(3, ceil(ground_spacing / resolution)));\n");
    fprintf(fp, "mesh.z = linspace(0, box_height, max(3, ceil(box_height / resolution)));\n");

    // Add refinement around ports and open ends
    double QL = f0 / (BW_MHz * 1e6);  // Loaded Q
    double FBW = BW_MHz / f0_MHz;  // Fractional BW
    double Qe = 1.0 / FBW;  // Approximate symmetric Qe for Butterworth
    double tap_norm = sqrt(R_ohm / (M_PI * QL * Qe));  // Normalized tap from grounded end

    double tap_z1 = 0.0, tap_zN = 0.0;
    double half_delta1 = 0.0, half_deltaN = 0.0;
    double delta_z = 1.0;  // Small z-span for tap gaps

    for (int i = 1; i <= ele; i++) {
        double len = rod_lengths[i - 1];
        double open_z = 0.0;
        if (i % 2 == 1) {  // Odd
            open_z = len;
        } else {  // Even
            open_z = lambda4_mm - len;
        }
        // Refine around open with clipped range to avoid out-of-bounds
        fprintf(fp, "mesh.z = sort(unique([mesh.z, linspace(max(0,%.2f -3), min(%.2f,%.2f +3), 5)]));\n", open_z, lambda4_mm, open_z);

        if (i == 1 || i == ele) {
            double tap_z = (i % 2 == 1) ? tap_norm * len : lambda4_mm - tap_norm * len;
            double z_start = (i % 2 == 1) ? 0.0 : lambda4_mm - len;
            double half_delta = delta_z / 2.0;
            half_delta = fmin(half_delta, tap_z - z_start);
            if (i == 1) {
                tap_z1 = tap_z;
                half_delta1 = half_delta;
            }
            if (i == ele) {
                tap_zN = tap_z;
                half_deltaN = half_delta;
            }
        }
    }

    // Calculate port positions
    double start_x1 = pos[1] - box_length / 2.0;
    double start_xN = pos[ele] - box_length / 2.0;
    double start_z1 = tap_z1 - half_delta1;
    double end_z1 = tap_z1 + half_delta1;
    double start_zN = tap_zN - half_deltaN;
    double end_zN = tap_zN + half_deltaN;

    // Validate port bounds
    if (start_x1 < -box_length / 2 || start_x1 > box_length / 2 || start_z1 < 0 || end_z1 > lambda4_mm ||
        start_xN < -box_length / 2 || start_xN > box_length / 2 || start_zN < 0 || end_zN > lambda4_mm) {
        fprintf(stderr, "Error: Port out of bounds (input: x=%.2f, z=[%.2f, %.2f], output: x=%.2f, z=[%.2f, %.2f])\n", start_x1, start_z1, end_z1, start_xN,
                start_zN, end_zN);
        fclose(fp);
        exit(1);
    }

    fprintf(fp, "mesh.z = sort(unique([mesh.z, linspace(%.2f, %.2f, 3), linspace(%.2f, %.2f, 3)]));\n", start_z1, end_z1, start_zN, end_zN);

    // Manual merge of close points to prevent tiny cells
    fprintf(fp, "mesh.z = sort(unique(mesh.z));\n");
    fprintf(fp, "temp = [];\n");
    fprintf(fp, "for ii=1:length(mesh.z)-1\n");
    fprintf(fp, "  if mesh.z(ii+1) - mesh.z(ii) > 0.1\n");
    fprintf(fp, "    temp = [temp, mesh.z(ii)];\n");
    fprintf(fp, "  end\n");
    fprintf(fp, "end\n");
    fprintf(fp, "temp = [temp, mesh.z(end)];\n");
    fprintf(fp, "mesh.z = temp;\n");

    // Smooth mesh with higher ratio
    fprintf(fp, "mesh.x = SmoothMeshLines(mesh.x, resolution / 5, 5.0);\n");
    fprintf(fp, "mesh.y = SmoothMeshLines(mesh.y, resolution / 5, 5.0);\n");
    fprintf(fp, "mesh.z = SmoothMeshLines(mesh.z, resolution / 5, 5.0);\n");

    // Define the rectangular grid
    fprintf(fp, "CSX = DefineRectGrid(CSX, unit, mesh);\n");

    // Add cavity box (PEC walls)
    fprintf(fp, "start = [-box_l/2, -ground_spacing/2, 0]; stop = [box_l/2, ground_spacing/2, box_height];\n");
    fprintf(fp, "CSX = AddBox(CSX, 'PEC', 10, start, stop);\n"); // Priority 10 for cavity

    // Add rods (cylinders)
    double delta_perp = 0.5;  // Assumed small perpendicular span for ports; adjust if needed
    for (int i = 1; i <= ele; i++) {
        double x = pos[i] - box_length / 2;
        double len = rod_lengths[i - 1];
        double z_start = (i % 2 == 1) ? 0.0 : lambda4_mm - len;
        double z_end = (i % 2 == 1) ? len : lambda4_mm;
        if (i == 1 || i == ele) {
            double tap_z = (i == 1) ? tap_z1 : tap_zN;
            double half_delta = (i == 1) ? half_delta1 : half_deltaN;
            // Lower part
            double lower_end = tap_z - half_delta;
            if (lower_end > z_start) {
                fprintf(fp, "start = [%.2f, 0, %.2f]; stop = [%.2f, 0, %.2f];\n", x, z_start, x, lower_end);
                fprintf(fp, "CSX = AddCylinder(CSX, 'PEC', 20, start, stop, %.2f);\n", D_mm / 2);
            }
            // Upper part
            double upper_start = tap_z + half_delta;
            if (upper_start < z_end) {
                fprintf(fp, "start = [%.2f, 0, %.2f]; stop = [%.2f, 0, %.2f];\n", x, upper_start, x, z_end);
                fprintf(fp, "CSX = AddCylinder(CSX, 'PEC', 20, start, stop, %.2f);\n", D_mm / 2);
            }
        } else {
            // Full cylinder for middle rods
            fprintf(fp, "start = [%.2f, 0, %.2f]; stop = [%.2f, 0, %.2f];\n", x, z_start, x, z_end);
            fprintf(fp, "CSX = AddCylinder(CSX, 'PEC', 20, start, stop, %.2f);\n", D_mm / 2);
        }
    }

    // Ports (across gaps)
    fprintf(fp, "delta_perp = %.2f;\n", delta_perp);
    fprintf(fp, "start_x1 = %.2f; start_z1 = %.2f; end_x1 = %.2f; end_z1 = %.2f;\n", start_x1, start_z1, start_x1, end_z1);
    fprintf(fp, "start_xN = %.2f; start_zN = %.2f; end_xN = %.2f; end_zN = %.2f;\n", start_xN, start_zN, start_xN, end_zN);
    fprintf(fp, "disp(['Input port: x = ', num2str(start_x1), ', z = [', num2str(start_z1), ', ', num2str(end_z1), ']']);\n");
    fprintf(fp, "disp(['Output port: x = ', num2str(start_xN), ', z = [', num2str(start_zN), ', ', num2str(end_zN), ']']);\n");
    fprintf(fp, "[CSX, port{1}] = AddLumpedPort(CSX, 30, 1, %.2f, [start_x1 - delta_perp, -delta_perp, start_z1], [start_x1 + delta_perp, delta_perp, end_z1], [0 0 1], 1);\n", R_ohm);  // Input active
    fprintf(fp, "[CSX, port{2}] = AddLumpedPort(CSX, 30, 2, %.2f, [start_xN - delta_perp, -delta_perp, start_zN], [start_xN + delta_perp, delta_perp, end_zN], [0 0 1]);\n", R_ohm);  // Output passive

    // Excitation and dumps (fields)
    fprintf(fp, "CSX = AddDump(CSX, 'Et', 'DumpType', 0, 'DumpMode', 0);\n");  // E-field time-domain
    fprintf(fp, "start = [-box_l/2, -ground_spacing/2, 0]; stop = [box_l/2, ground_spacing/2, box_height];\n");
    fprintf(fp, "CSX = AddBox(CSX, 'Et', 0, start, stop);\n");

    // Run simulation
    fprintf(fp, "WriteOpenEMS([Sim_Path, '/bpf.xml'], FDTD, CSX);\n");
    fprintf(fp, "disp('XML written, verifying geometry...'); pause(1);\n");  // Debug pause
    fprintf(fp, "RunOpenEMS(Sim_Path, 'bpf.xml');\n");

    // Post-process S-params
    fprintf(fp, "f = linspace(f_min, f_max, 2001);\n");
    fprintf(fp, "port = calcPort(port, Sim_Path, f);\n");
    fprintf(fp, "s11 = port{1}.uf.ref ./ port{1}.uf.inc;\n");
    fprintf(fp, "s21 = port{2}.uf.ref ./ port{1}.uf.inc;\n");
    fprintf(fp, "figure; plot(f/1e9, 20*log10(abs(s11)), 'k--', 'DisplayName', 'S11');\n");
    fprintf(fp, "hold on; plot(f/1e9, 20*log10(abs(s21)), 'b-', 'DisplayName', 'S21');\n");
    fprintf(fp, "xlabel('Frequency (GHz)'); ylabel('Magnitude (dB)'); legend; grid on;\n");
    fprintf(fp, "%% Save S-parameters to separate file\n");
    fprintf(fp, "data = [f(:)/1e9, real(s11(:)), imag(s11(:)), real(s21(:)), imag(s21(:))];\n");
    fprintf(fp, "fid = fopen('s_params.csv', 'w');\n");
    fprintf(fp, "fprintf(fid, 'Frequency (GHz),Re(S11),Im(S11),Re(S21),Im(S21)\\n');\n");
    fprintf(fp, "fclose(fid);\n");
    fprintf(fp, "dlmwrite('s_params.csv', data, '-append', 'delimiter', ',', 'precision', '%%.6f');\n");
    fprintf(fp, "disp('S-parameters saved to s_params.csv');\n");

    fclose(fp);
    printf("openEMS script written to %s\n", filename);
}
