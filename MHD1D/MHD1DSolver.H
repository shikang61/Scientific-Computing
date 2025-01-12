#ifndef MHD1DSOLVER_H
#define MHD1DSOLVER_H

#include "EoS.H"
#include "method.H"
#include "material.H"
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>

void solveEuler(std::string scheme, std::string initial, std::string geometry);
std::ostringstream generate_filename(method_ptr &numMethod, std::string initial, std::string &geometry, EoS_ptr &EoS);
void writeToFile(std::ofstream &output, double t, VectOfDouble xCells, VectofVectDouble &u, EoS_ptr &Eo);
void plot(std::ostringstream &filename, std::string initial);
void ConstructFVDomain(double &dx, VectOfDouble &xCells);
void setInitialConditions(const VectOfDouble &xCells, VectofVectDouble &u, double &tStop, std::string &initial, method_ptr &numMethod, EoS_ptr &EoS);
void computeDt(VectofVectDouble &u, double &dx, double &dt, int &step, EoS_ptr &EoS);
void applyBoundaryConditions(VectofVectDouble &u, std::array<std::string, nVar> boundary, int nCells);

#endif