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

    // Find max rod length for box_height with margin
    double max_l = 0.0;
    for (int i = 0; i < ele; i++) {
        if (rod_lengths[i] > max_l)
            max_l = rod_lengths[i];
    }
    double box_height = max_l * 1.1;  // 10% margin for open ends

    // Compute g values
    double *g = (double*) malloc((ele + 2) * sizeof(double));
    if (ripple_dB <= 0.0) {
        g[0] = 1.0;
        for (int k = 1; k <= ele; k++) {
            g[k] = 2.0 * sin((2.0 * k - 1.0) * M_PI / (2.0 * ele));
        }
        g[ele + 1] = 1.0;
    } else {
        double eps = sqrt(pow(10.0, ripple_dB / 10.0) - 1.0);
        double beta = asinh(1.0 / eps) / ele;
        g[0] = 1.0;
        double sinhb = sinh(beta);
        double a[64] = { 0 }, b[64];
        for (int k = 1; k <= ele; k++) {
            a[k - 1] = sin((2.0 * k - 1.0) * M_PI / (2.0 * ele));  // 0-based
        }
        for (int k = 1; k <= ele - 1; k++) {
            double s = sin(k * M_PI / ele);
            b[k - 1] = sinhb * sinhb + s * s;
        }
        g[1] = (2.0 * a[0]) / sinhb;
        for (int k = 2; k <= ele; k++) {
            g[k] = (4.0 * a[k - 2] * a[k - 1]) / (b[k - 2] * g[k - 1]);
        }
        if ((ele % 2) == 0) {
            g[ele + 1] = 1.0;
        } else {
            double tb2 = tanh(0.5 * beta);
            g[ele + 1] = 1.0 / (tb2 * tb2);
        }
    }

    double QL = f0_MHz / BW_MHz;
    double Qe1 = g[0] * QL / g[1];
    double QeN = g[ele + 1] * QL / g[ele];

    // Zr for round rod
    double Zr = 60.0 * log(4.0 * H_mm / (M_PI * D_mm));

    // Tap distances
    double arg1 = M_PI * R_ohm / (4.0 * Zr * Qe1);
    if (arg1 > 1.0)
        arg1 = 1.0;
    double sin_theta1 = sqrt(arg1);
    double theta1 = asin(sin_theta1);
    double d1 = (2.0 * rod_lengths[0] / M_PI) * theta1;

    double argN = M_PI * R_ohm / (4.0 * Zr * QeN);
    if (argN > 1.0)
        argN = 1.0;
    double sin_thetaN = sqrt(argN);
    double thetaN = asin(sin_thetaN);
    double dN = (2.0 * rod_lengths[ele - 1] / M_PI) * thetaN;

    // Tap positions based on shorted end
    double tap_z1 = d1;  // First rod (i=0 even) shorted at z=0
    double tap_zN;
    if ((ele - 1) % 2 == 0) {
        tap_zN = dN;  // Shorted at z=0
    } else {
        tap_zN = box_height - dN;  // Shorted at box_height
    }

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
    fprintf(fp, "box_height = %.2f; %% box height in z (rod direction)\n", box_height);
    fprintf(fp, "Sim_Path = 'sim_bpf';\n");
    fprintf(fp, "mkdir(Sim_Path);\n");
    fprintf(fp, "disp(['box_l = ', num2str(box_l), ' mm']);\n");
    fprintf(fp, "disp(['ground_spacing = ', num2str(ground_spacing), ' mm']);\n");
    fprintf(fp, "disp(['box_height = ', num2str(box_height), ' mm']);\n");
    fprintf(fp, "FDTD = InitFDTD('NrTS', 10000000, 'EndCriteria', 1e-5);\n");
    fprintf(fp, "FDTD = SetGaussExcite(FDTD, f0, (f_max - f_min)/2);\n");
    fprintf(fp, "BC = {'PEC' 'PEC' 'PEC' 'PEC' 'PEC' 'PEC'}; %% All PEC for closed metal cavity\n");
    fprintf(fp, "FDTD = SetBoundaryCond(FDTD, BC);\n");

    // CSX and mesh
    fprintf(fp, "CSX = InitCSX();\n");
    fprintf(fp, "CSX = AddMetal(CSX, 'PEC'); %% Define Perfect Electric Conductor (PEC)\n");
    fprintf(fp, "resolution = c0 / (f_max) / unit / 30; %% Finer mesh lambda/30 for better accuracy\n");
    fprintf(fp, "mesh.x = linspace(-box_l/2, box_l/2, max(3, ceil(box_l / resolution)));\n");
    fprintf(fp, "mesh.y = linspace(-ground_spacing/2, ground_spacing/2, max(3, ceil(ground_spacing / resolution)));\n");
    fprintf(fp, "mesh.z = linspace(0, box_height, max(3, ceil(box_height / resolution)));\n");
    fprintf(fp,
            "mesh.y = sort(unique([mesh.y, linspace(-ground_spacing/2, -ground_spacing/2 + 2, 5), linspace(ground_spacing/2 - 2, ground_spacing/2, 5), linspace(-%.2f/2 - 1, -%.2f/2 + 1, 5), linspace(%.2f/2 - 1, %.2f/2 + 1, 5)]));\n",
            D_mm, D_mm, D_mm, D_mm);

    // Add rods
    double x_shift = box_length / 2.0;
    for (int i = 0; i < ele; i++) {
        double x = pos[i + 1] - x_shift;
        double radius = D_mm / 2.0;
        double l = rod_lengths[i];
        double z_start, z_end;
        if (i % 2 == 0) {
            z_start = 0.0;
            z_end = l;
        } else {
            z_start = box_height - l;
            z_end = box_height;
        }
        fprintf(fp, "CSX = AddCylinder(CSX, 'PEC', 10, [%.2f 0 %.2f], [%.2f 0 %.2f], %.2f);\n", x, z_start, x, z_end, radius);
    }

    // Ports (tapped feed from side wall to rod surface in y direction)
    double delta_x = 0.5;  // Small extent in x
    double delta_z = 0.5;  // Small extent in z
    double y_wall_in = -H_mm / 2.0;
    double y_rod_in = -D_mm / 2.0;
    double y_wall_out = H_mm / 2.0;
    double y_rod_out = D_mm / 2.0;
    double start_x1 = pos[1] - x_shift;
    double start_xN = pos[ele] - x_shift;

    fprintf(fp, "delta_x = %.2f;\n", delta_x);
    fprintf(fp, "delta_z = %.2f;\n", delta_z);
    fprintf(fp, "start_x1 = %.2f; tap_z1 = %.2f;\n", start_x1, tap_z1);
    fprintf(fp, "start_xN = %.2f; tap_zN = %.2f;\n", start_xN, tap_zN);
    fprintf(fp, "disp(['Input port: x = ', num2str(start_x1), ', y = [', num2str(%.2f), ', ', num2str(%.2f), '], z = ', num2str(tap_z1)]);\n", y_wall_in,
            y_rod_in);
    fprintf(fp, "disp(['Output port: x = ', num2str(start_xN), ', y = [', num2str(%.2f), ', ', num2str(%.2f), '], z = ', num2str(tap_zN)]);\n", y_wall_out,
            y_rod_out);
    fprintf(fp,
            "[CSX, port{1}] = AddLumpedPort(CSX, 30, 1, %.2f, [start_x1 - delta_x, %.2f, tap_z1 - delta_z], [start_x1 + delta_x, %.2f, tap_z1 + delta_z], [0 1 0], 1);\n",
            R_ohm, y_wall_in, y_rod_in);  // Input active, dir [0 1 0]
    fprintf(fp,
            "[CSX, port{2}] = AddLumpedPort(CSX, 30, 2, %.2f, [start_xN - delta_x, %.2f, tap_zN - delta_z], [start_xN + delta_x, %.2f, tap_zN + delta_z], [0 1 0]);\n",
            R_ohm, y_rod_out, y_wall_out);  // Output passive, dir [0 1 0]

    // Excitation and dumps (fields)
    fprintf(fp, "CSX = AddDump(CSX, 'Et', 'DumpType', 0, 'DumpMode', 0);\n");  // E-field time-domain
    fprintf(fp, "start = [-box_l/2, -ground_spacing/2, 0]; stop = [box_l/2, ground_spacing/2, box_height];\n");
    fprintf(fp, "CSX = AddBox(CSX, 'Et', 0, start, stop);\n");

    // Mesh smoothing and definition
    fprintf(fp, "mesh = SmoothMesh(mesh, resolution / 4, 1.5);\n");  // Adjust as needed
    fprintf(fp, "CSX = DefineRectGrid(CSX, unit, mesh);\n");

    // Run simulation
    fprintf(fp, "WriteOpenEMS([Sim_Path, '/bpf.xml'], FDTD, CSX);\n");
    fprintf(fp, "disp('XML written, verifying geometry...'); pause(1);\n");
    fprintf(fp, "RunOpenEMS(Sim_Path, 'bpf.xml');\n");

    // Post-process S-params
    fprintf(fp, "f = linspace(f_min, f_max, 2001);\n");
    fprintf(fp, "port = calcPort(port, Sim_Path, f);\n");
    fprintf(fp, "for ii=1:length(f)\n");
    fprintf(fp, "  if abs(port{1}.uf.inc(ii)) > 1e-10\n");
    fprintf(fp, "    s11(ii) = port{1}.uf.ref(ii) ./ port{1}.uf.inc(ii);\n");
    fprintf(fp, "    s21(ii) = port{2}.uf.tot(ii) ./ port{1}.uf.inc(ii);\n");  // Correct to uf.tot for transmitted
    fprintf(fp, "  else\n");
    fprintf(fp, "    s11(ii) = 1; s21(ii) = 0;\n");
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

    free(g);
    fclose(fp);
    printf("openEMS script written to %s\n", filename);
}
