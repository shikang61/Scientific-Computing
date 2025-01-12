#ifndef EULER1DSOLVER_H
#define EULER1DSOLVER_H

#include "EoS.H"
#include "method.H"
#include "material.H"
#include <iostream>
#include <fstream>
#include <sstream>

void solveEuler(std::string method, std::string initial, std::string geometry, Mat_ptr &material);
std::ostringstream generate_filename(Method_ptr &numMethod, std::string initial, std::string &geometry, EoS_ptr &EoS);
void writeToFile(std::ofstream &output, double t, VectOfDouble xCells, VectofVectDouble &u, VectofVectDouble &w_ExactRiemann, EoS_ptr &EoS);
void plot(std::ostringstream &filename, std::string initial);
void ConstructFVDomain(double &dx, VectOfDouble &xCells);
void setInitialConditions(const VectOfDouble &xCells, VectofVectDouble &u, double &tStop, VectofVectDouble &w_ExactRiemann, std::string &initial, Method_ptr &numMethod, bool exactSolution, EoS_ptr &EoS);
void computeDt(VectofVectDouble &u, double &dx, double &dt, int &step, EoS_ptr &EoS);
void applyBoundaryConditions(VectofVectDouble &u, std::array<std::string, nVar> boundary, int nCells);

#endif