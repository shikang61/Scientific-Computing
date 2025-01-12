#ifndef EULER1D_H
#define EULER1D_H

#include <array>
#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <cmath>
#include <math.h>
#include <algorithm>
#include <memory>
#include "EoS.H"

double c_s_star(double p_star, VectOfDouble &w, EoS_ptr &EoS);
double A(VectOfDouble &w, EoS_ptr &EoS);
double B(VectOfDouble &w, EoS_ptr &EoS);
double f_K(double p, VectOfDouble &w, EoS_ptr &EoS);
double df_K(double p, VectOfDouble &w, EoS_ptr &EoS);
double f(double p, VectOfDouble &w_L, VectOfDouble &w_R, EoS_ptr &EoS);
double df(double p, VectOfDouble &w_L, VectOfDouble &w_R, EoS_ptr &EoS);
double newtonRaphson(double initialGuess, double tolerance, int maxIterations, VectOfDouble &w_L, VectOfDouble &w_R, EoS_ptr &EoS);
VectOfDouble raref(double S, VectOfDouble &w, std::string side, EoS_ptr &EoS);
double rho_star_shock(double p_star, VectOfDouble w, EoS_ptr &EoS);
double S_shock(double p_star, VectOfDouble w, std::string side, EoS_ptr &EoS);
double S_H(VectOfDouble w, std::string side, EoS_ptr &EoS);
double S_T(VectOfDouble w_star, double p_star, std::string side, EoS_ptr &EoS);
double rho_star_raref(double p_star, VectOfDouble w, EoS_ptr &EoS);
VectOfDouble exactRiemannSolver(double x, double mid, double tStop, VectOfDouble &w_L_initial, VectOfDouble &w_R_initial, EoS_ptr &EoS);

std::pair<double, double> wave_speed_estimate(VectOfDouble &u_L, VectOfDouble &u_R, EoS_ptr &EoS);
#endif