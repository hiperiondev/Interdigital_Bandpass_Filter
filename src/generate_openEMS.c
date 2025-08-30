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
    fprintf(fp, "FDTD = InitFDTD('NrTS', 5000000, 'EndCriteria', 1e-4);\n");
    fprintf(fp, "FDTD = SetGaussExcite(FDTD, f0, (f_max - f_min)/2);\n");
    fprintf(fp, "BC = {'PEC' 'PEC' 'PEC' 'PEC' 'PEC' 'PEC'}; %% All PEC for closed metal cavity\n");
    fprintf(fp, "FDTD = SetBoundaryCond(FDTD, BC);\n");

    // CSX and mesh
    fprintf(fp, "CSX = InitCSX();\n");
    fprintf(fp, "CSX = AddMetal(CSX, 'PEC'); %% Define Perfect Electric Conductor (PEC)\n");
    fprintf(fp, "resolution = c0 / (f_max) / unit / 12; %% Base mesh lambda/12 (~2.5mm) for reduced cell count\n");
    fprintf(fp, "mesh.x = linspace(-box_l/2, box_l/2, max(3, ceil(box_l / resolution)));\n");
    fprintf(fp, "mesh.y = linspace(-ground_spacing/2, ground_spacing/2, max(3, ceil(ground_spacing / resolution)));\n");
    fprintf(fp, "mesh.z = linspace(0, box_height, max(3, ceil(box_height / resolution)));\n");

    // Add refinement around ports and open ends
    double QL = f0 / (BW_MHz * 1e6);  // Loaded Q
    double FBW = BW_MHz / f0_MHz;  // Fractional BW
    double Qe = 1.0 / FBW;  // Approximate symmetric Qe for Butterworth
    double tap_norm = sqrt(R_ohm / (M_PI * QL * Qe));  // Normalized tap from grounded end

    double delta_z = 1.5;  // Z-span for ports (unchanged, sufficient for mesh)
    double delta_perp = 0.5;  // X/Y span for 1D z-directed ports

    // Define port variables
    double start_x1 = 0.0, start_z1 = 0.0, end_z1 = 0.0;
    double start_xN = 0.0, start_zN = 0.0, end_zN = 0.0;

    for (int i = 1; i <= ele; i++) {
        double x = pos[i] - box_length / 2;
        double len = rod_lengths[i - 1];
        double open_z = 0.0;
        if (i % 2 == 1) {  // Odd
            open_z = len;
        } else {  // Even
            open_z = lambda4_mm - len;
        }
        // Refine around open end (tighter, fewer points)
        fprintf(fp, "mesh.z = sort(unique([mesh.z, linspace(max(0,%.2f -1.0), min(%.2f,%.2f +1.0), 2)]));\n", open_z, lambda4_mm, open_z);

        if (i == 1 || i == ele) {
            double tap_z = (i % 2 == 1) ? tap_norm * len : lambda4_mm - tap_norm * len;
            double z_start = (i % 2 == 1) ? 0.0 : lambda4_mm - len;
            double z_end = (i % 2 == 1) ? len : lambda4_mm;
            double half_delta = delta_z / 2.0;
            double lower_end = tap_z - half_delta;
            double upper_start = tap_z + half_delta;

            // Set port positions
            if (i == 1) {
                start_x1 = x;
                start_z1 = lower_end;
                end_z1 = upper_start;
            }
            if (i == ele) {
                start_xN = x;
                start_zN = lower_end;
                end_zN = upper_start;
            }

            // Add refinements for tap_z, port z-edges, and x-position (tighter, fewer points)
            fprintf(fp, "mesh.z = sort(unique([mesh.z, linspace(max(0,%.2f -1.0), min(%.2f,%.2f +1.0), 2)]));\n", tap_z, lambda4_mm, tap_z);
            fprintf(fp, "mesh.z = sort(unique([mesh.z, %.2f, %.2f]));\n", lower_end, upper_start);
            fprintf(fp, "mesh.x = sort(unique([mesh.x, linspace(%.2f -1.0, %.2f +1.0, 2)]));\n", x, x);
            // Debug: Log x-mesh after each refinement
            fprintf(fp, "fid = fopen([Sim_Path, '/mesh_debug.txt'], 'a'); fprintf(fid, 'x-mesh after rod %d: %%s\\n', num2str(mesh.x)); fclose(fid);\n", i);

            // Lower part
            if (lower_end > z_start) {
                fprintf(fp, "start = [%.2f, 0, %.2f]; stop = [%.2f, 0, %.2f];\n", x, z_start, x, lower_end);
                fprintf(fp, "CSX = AddCylinder(CSX, 'PEC', 20, start, stop, %.2f);\n", D_mm / 2);
            }
            // Upper part
            if (upper_start < z_end) {
                fprintf(fp, "start = [%.2f, 0, %.2f]; stop = [%.2f, 0, %.2f];\n", x, upper_start, x, z_end);
                fprintf(fp, "CSX = AddCylinder(CSX, 'PEC', 20, start, stop, %.2f);\n", D_mm / 2);
            }
        } else {
            // Full cylinder for middle rods
            double z_start = (i % 2 == 1) ? 0.0 : lambda4_mm - len;
            double z_end = (i % 2 == 1) ? len : lambda4_mm;
            fprintf(fp, "start = [%.2f, 0, %.2f]; stop = [%.2f, 0, %.2f];\n", x, z_start, x, z_end);
            fprintf(fp, "CSX = AddCylinder(CSX, 'PEC', 20, start, stop, %.2f);\n", D_mm / 2);
            // Refine middle rod ends (optional, minimal)
            fprintf(fp, "mesh.z = sort(unique([mesh.z, linspace(max(0,%.2f -1.0), min(%.2f,%.2f +1.0), 2)]));\n", open_z, lambda4_mm, open_z);
            fprintf(fp, "fid = fopen([Sim_Path, '/mesh_debug.txt'], 'a'); fprintf(fid, 'x-mesh after rod %d: %%s\\n', num2str(mesh.x)); fclose(fid);\n", i);
        }
    }

    // Consolidate mesh and log final sizes
    fprintf(fp, "mesh.x = sort(unique(mesh.x));\n");
    fprintf(fp, "mesh.y = sort(unique(mesh.y));\n");
    fprintf(fp, "mesh.z = sort(unique(mesh.z));\n");
    fprintf(fp, "fid = fopen([Sim_Path, '/mesh_debug.txt'], 'a');\n");
    fprintf(fp, "fprintf(fid, 'Final mesh.x (%%d points): %%s\\n', length(mesh.x), num2str(mesh.x));\n");
    fprintf(fp, "fprintf(fid, 'Final mesh.y (%%d points): %%s\\n', length(mesh.y), num2str(mesh.y));\n");
    fprintf(fp, "fprintf(fid, 'Final mesh.z (%%d points): %%s\\n', length(mesh.z), num2str(mesh.z));\n");
    fprintf(fp, "fprintf(fid, 'Estimated cell count: %%d\\n', length(mesh.x) * length(mesh.y) * length(mesh.z));\n");
    fprintf(fp, "fclose(fid);\n");

    // Smooth mesh with higher ratio
    fprintf(fp, "mesh.x = SmoothMeshLines(mesh.x, resolution / 4, 15.0);\n");
    fprintf(fp, "mesh.y = SmoothMeshLines(mesh.y, resolution / 4, 15.0);\n");
    fprintf(fp, "mesh.z = SmoothMeshLines(mesh.z, resolution / 4, 15.0);\n");

    // Define the rectangular grid
    fprintf(fp, "CSX = DefineRectGrid(CSX, unit, mesh);\n");

    // Add cavity box (PEC walls)
    fprintf(fp, "start = [-box_l/2, -ground_spacing/2, 0]; stop = [box_l/2, ground_spacing/2, box_height];\n");
    fprintf(fp, "CSX = AddBox(CSX, 'PEC', 10, start, stop);\n");  // Priority 10 for cavity

    // Ports (across gaps)
    fprintf(fp, "delta_perp = %.2f;\n", delta_perp);
    fprintf(fp, "start_x1 = %.2f; start_z1 = %.2f; end_x1 = %.2f; end_z1 = %.2f;\n", start_x1, start_z1, start_x1, end_z1);
    fprintf(fp, "start_xN = %.2f; start_zN = %.2f; end_xN = %.2f; end_zN = %.2f;\n", start_xN, start_zN, start_xN, end_zN);
    fprintf(fp, "disp(['Input port: x = ', num2str(start_x1), ', z = [', num2str(start_z1), ', ', num2str(end_z1), ']']);\n");
    fprintf(fp, "disp(['Output port: x = ', num2str(start_xN), ', z = [', num2str(start_zN), ', ', num2str(end_zN), ']']);\n");
    fprintf(fp,
            "[CSX, port{1}] = AddLumpedPort(CSX, 30, 1, %.2f, [start_x1 - delta_perp, -delta_perp, start_z1], [start_x1 + delta_perp, delta_perp, end_z1], [0 0 1], 1);\n",
            R_ohm);  // Input active
    fprintf(fp,
            "[CSX, port{2}] = AddLumpedPort(CSX, 30, 2, %.2f, [start_xN - delta_perp, -delta_perp, start_zN], [start_xN + delta_perp, delta_perp, end_zN], [0 0 1]);\n",
            R_ohm);  // Output passive

    // Excitation and dumps (fields)
    fprintf(fp, "CSX = AddDump(CSX, 'Et', 'DumpType', 0, 'DumpMode', 0);\n");  // E-field time-domain
    fprintf(fp, "start = [-box_l/2, -ground_spacing/2, 0]; stop = [box_l/2, ground_spacing/2, box_height];\n");
    fprintf(fp, "CSX = AddBox(CSX, 'Et', 0, start, stop);\n");

    // Run simulation
    fprintf(fp, "WriteOpenEMS([Sim_Path, '/bpf.xml'], FDTD, CSX);\n");
    fprintf(fp, "disp('XML written, verifying geometry...'); pause(1);\n");  // Debug pause
    fprintf(fp, "RunOpenEMS(Sim_Path, 'bpf.xml');\n");

    // Post-process S-params (with checks for division by zero to avoid NaNs)
    fprintf(fp, "f = linspace(f_min, f_max, 2001);\n");
    fprintf(fp, "port = calcPort(port, Sim_Path, f);\n");
    fprintf(fp, "for ii=1:length(f)\n");
    fprintf(fp, "  if abs(port{1}.uf.inc(ii)) > 1e-10\n");
    fprintf(fp, "    s11(ii) = port{1}.uf.ref(ii) ./ port{1}.uf.inc(ii);\n");
    fprintf(fp, "    s21(ii) = port{2}.uf.ref(ii) ./ port{1}.uf.inc(ii);\n");
    fprintf(fp, "  else\n");
    fprintf(fp, "    s11(ii) = 1; s21(ii) = 0; %% Fallback for zero incident\n");
    fprintf(fp, "  end\n");
    fprintf(fp, "end\n");
    fprintf(fp, "figure; plot(f/1e9, 20*log10(abs(s11)), 'k--', 'DisplayName', 'S11');\n");
    fprintf(fp, "hold on; plot(f/1e9, 20*log10(abs(s21)), 'b-', 'DisplayName', 'S21');\n");
    fprintf(fp, "xlabel('Frequency (GHz)'); ylabel('Magnitude (dB)'); legend; grid on;\n");
    fprintf(fp, "%% Save S-parameters to separate file\n");
    fprintf(fp, "data = [f(:)/1e9, real(s11(:)), imag(s11(:)), abs(s11(:)), real(s21(:)), imag(s21(:)), abs(s21(:))];\n");
    fprintf(fp, "fid = fopen('s_params.csv', 'w');\n");
    fprintf(fp, "fprintf(fid, 'Frequency (GHz),Re(S11),Im(S11),|S11|,Re(S21),Im(S21),|S21|\\n');\n");
    fprintf(fp, "fclose(fid);\n");
    fprintf(fp, "dlmwrite('s_params.csv', data, '-append', 'delimiter', ',', 'precision', '%%.6f');\n");
    fprintf(fp, "disp('S-parameters saved to s_params.csv');\n");

    fclose(fp);
    printf("openEMS script written to %s\n", filename);
}
