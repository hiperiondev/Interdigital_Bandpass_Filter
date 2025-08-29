/*
 * rod.h
 *
 *  Created on: 27 ago 2025
 *      Author: egonzalez
 */

#ifndef ROD_H_
#define ROD_H_

double interdigital_round_rod_length_corrected(double f0, double d, double h, double g, int keep_legacy_factor, double legacy_factor_value);
double rod_length(double f0, double h, double d, double length_factor);

#endif /* ROD_H_ */
