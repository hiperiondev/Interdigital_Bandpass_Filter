/*
 * generate_openEMS.h
 *
 *  Created on: 28 ago 2025
 *      Author: egonzalez
 */

#ifndef GENERATE_OPENEMS_H_
#define GENERATE_OPENEMS_H_

void generate_openEMS_script(const char *filename, double f0_MHz, double BW_MHz, double R_ohm, double H_mm, double D_mm, double E_mm, int ele, double ripple_dB,
        double *pos, double *gap, double *rod_lengths);

#endif /* GENERATE_OPENEMS_H_ */
