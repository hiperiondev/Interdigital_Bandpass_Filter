/*
 * generate_openEMS.c
 *
 *  Created on: 28 ago 2025
 *      Author: egonzalez
 */

/*
 * generate_openEMS.c
 * Generates MATLAB script for openEMS simulation of interdigital BPF.
 * Call after computing parameters in test_filter.c.
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
    fprintf(fp, "FDTD = InitFDTD('NrTS', 1e7, 'EndCriteria', 1e-5);\n");
    fprintf(fp, "FDTD = SetGaussExcite(FDTD, f0, (f_max - f_min)/2);\n");
    fprintf(fp, "BC = {'PEC' 'PEC' 'PEC' 'PEC' 'PEC' 'PEC'}; %% All PEC for closed metal cavity\n");
    fprintf(fp, "FDTD = SetBoundaryCond(FDTD, BC);\n");

    // CSX and mesh
    fprintf(fp, "CSX = InitCSX();\n");
    fprintf(fp, "resolution = c0 / (f_max) / unit / 15; %% lambda/15 for stability\n");
    fprintf(fp, "mesh.x = linspace(-box_l/2, box_l/2, ceil(box_l / resolution));\n");
    fprintf(fp, "mesh.y = linspace(-ground_spacing/2, ground_spacing/2, ceil(ground_spacing / resolution));\n");
    fprintf(fp, "mesh.z = linspace(0, box_height, ceil(box_height / resolution));\n");

    // Add refinement around ports and open ends
    double QL = f0 / (BW_MHz * 1e6);  // Loaded Q
    double FBW = BW_MHz / f0_MHz;  // Fractional BW
    double Qe = 1.0 / FBW;  // Approximate symmetric Qe for Butterworth
    double tap_norm = sqrt(R_ohm / (M_PI * QL * Qe));  // Normalized tap from grounded end

    // Collect refinement points: taps and opens
    for (int i = 1; i <= ele; i++) {
        double len = rod_lengths[i - 1];
        double tap_z = 0.0;
        double open_z = 0.0;
        if (i % 2 == 1) {  // Odd
            tap_z = tap_norm * len;
            open_z = len;
        } else {  // Even
            tap_z = lambda4_mm - tap_norm * len;
            open_z = lambda4_mm - len;
        }
        // Refine around open
        fprintf(fp, "mesh.z = sort(unique([mesh.z, linspace(%.2f -2, %.2f +2, 5)]));\n", open_z, open_z);
        if (i == 1 || i == ele) {
            // Refine around tap
            fprintf(fp, "mesh.z = sort(unique([mesh.z, linspace(%.2f -2, %.2f +2, 5)]));\n", tap_z, tap_z);
        }
    }

    // Add mesh lines for rods and ports
    fprintf(fp, "mesh.y = sort(unique([mesh.y, 0]));\n");
    for (int i = 1; i <= ele; i++) {
        double x = pos[i] - box_length / 2;
        fprintf(fp, "mesh.x = sort(unique([mesh.x, %.2f]));\n", x);
    }

    fprintf(fp, "mesh = SmoothMesh(mesh, resolution, 1.3, 0.5);\n");
    fprintf(fp, "CSX = DefineRectGrid(CSX, unit, mesh);\n");

    // Add metal for rods
    fprintf(fp, "CSX = AddMetal(CSX, 'PEC');\n");

    double delta_z = 1.0;  // Small gap for port (mm)
    double half_delta = delta_z / 2.0;

    for (int i = 1; i <= ele; i++) {
        double x = pos[i] - box_length / 2;
        double len = rod_lengths[i - 1];
        double z_ground, z_open, tap_z;
        if (i % 2 == 1) {  // Odd: grounded at bottom (z=0), open at top
            z_ground = 0.0;
            z_open = len;
            tap_z = z_ground + tap_norm * len;
        } else {  // Even: grounded at top (z=box_height), open at bottom
            z_ground = lambda4_mm;
            z_open = lambda4_mm - len;
            tap_z = z_ground - tap_norm * len;  // Tap from grounded downward
        }
        double z_start = fmin(z_ground, z_open);
        double z_end = fmax(z_ground, z_open);

        // Clamp to [0, lambda4_mm]
        z_start = fmax(0.0, z_start);
        z_end = fmin(lambda4_mm, z_end);

        // Validate rod
        if (x < -box_length / 2 || x > box_length / 2 || z_start >= z_end || len <= 0 || len > lambda4_mm) {
            fprintf(stderr, "Error: Rod %d out of bounds or invalid (x=%.2f, z=[%.2f, %.2f], len=%.2f > height=%.2f?)\n", i, x, z_start, z_end, len,
                    lambda4_mm);
            fclose(fp);
            exit(1);
        }

        fprintf(fp, "disp(['Rod %d: x = ', num2str(%.2f), ', z_start = ', num2str(%.2f), ', length = ', num2str(%.2f)]);\n", i, x, z_start, len);

        if (i == 1 || i == ele) {  // Split for input/output rods at tap
        // Lower part: from z_start to tap_z - half_delta
            double lower_end = tap_z - half_delta;
            if (lower_end > z_start) {
                fprintf(fp, "start = [%.2f, 0, %.2f]; stop = [%.2f, 0, %.2f];\n", x, z_start, x, lower_end);
                fprintf(fp, "CSX = AddCylinder(CSX, 'PEC', 20, start, stop, %.2f);\n", D_mm / 2);
            }

            // Upper part: from tap_z + half_delta to z_end
            double upper_start = tap_z + half_delta;
            if (upper_start < z_end) {
                fprintf(fp, "start = [%.2f, 0, %.2f]; stop = [%.2f, 0, %.2f];\n", x, upper_start, x, z_end);
                fprintf(fp, "CSX = AddCylinder(CSX, 'PEC', 20, start, stop, %.2f);\n", D_mm / 2);
            }
        } else {  // Full cylinder for middle rods
            fprintf(fp, "start = [%.2f, 0, %.2f]; stop = [%.2f, 0, %.2f];\n", x, z_start, x, z_end);
            fprintf(fp, "CSX = AddCylinder(CSX, 'PEC', 20, start, stop, %.2f);\n", D_mm / 2);
        }
    }

    // Ports (across gaps)
    double start_x1 = pos[1] - box_length / 2;
    double len1 = rod_lengths[0];
    double tap_z1 = (1 % 2 == 1) ? tap_norm * len1 : lambda4_mm - tap_norm * len1;
    double start_z1 = tap_z1 - half_delta;
    double end_z1 = tap_z1 + half_delta;

    double start_xN = pos[ele] - box_length / 2;
    double lenN = rod_lengths[ele - 1];
    double tap_zN = (ele % 2 == 1) ? tap_norm * lenN : lambda4_mm - tap_norm * lenN;
    double start_zN = tap_zN - half_delta;
    double end_zN = tap_zN + half_delta;

    // Ensure start < end
    if (start_z1 > end_z1) {
        double tmp = start_z1;
        start_z1 = end_z1;
        end_z1 = tmp;
    }
    if (start_zN > end_zN) {
        double tmp = start_zN;
        start_zN = end_zN;
        end_zN = tmp;
    }

    // Validate ports
    if (start_x1 < -box_length / 2 || start_x1 > box_length / 2 || start_z1 < 0 || end_z1 > lambda4_mm || start_xN < -box_length / 2
            || start_xN > box_length / 2 || start_zN < 0 || end_zN > lambda4_mm) {
        fprintf(stderr, "Error: Port out of bounds (input: x=%.2f, z=[%.2f, %.2f], output: x=%.2f, z=[%.2f, %.2f])\n", start_x1, start_z1, end_z1, start_xN,
                start_zN, end_zN);
        fclose(fp);
        exit(1);
    }

    // Add exact mesh lines for port z positions
    fprintf(fp, "mesh.z = sort(unique([mesh.z, %.2f, %.2f, %.2f, %.2f]));\n", start_z1, end_z1, start_zN, end_zN);

    fprintf(fp, "start_x1 = %.2f; start_z1 = %.2f; end_x1 = %.2f; end_z1 = %.2f;\n", start_x1, start_z1, start_x1, end_z1);
    fprintf(fp, "start_xN = %.2f; start_zN = %.2f; end_xN = %.2f; end_zN = %.2f;\n", start_xN, start_zN, start_xN, end_zN);
    fprintf(fp, "disp(['Input port: x = ', num2str(start_x1), ', z = [', num2str(start_z1), ', ', num2str(end_z1), ']']);\n");
    fprintf(fp, "disp(['Output port: x = ', num2str(start_xN), ', z = [', num2str(start_zN), ', ', num2str(end_zN), ']']);\n");
    fprintf(fp, "[CSX, port{1}] = AddLumpedPort(CSX, 30, 1, %.2f, [start_x1, 0, start_z1], [end_x1, 0, end_z1], [0 0 1], 1);\n", R_ohm);  // Input active
    fprintf(fp, "[CSX, port{2}] = AddLumpedPort(CSX, 30, 2, %.2f, [start_xN, 0, start_zN], [end_xN, 0, end_zN], [0 0 1]);\n", R_ohm);  // Output passive

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
