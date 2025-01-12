#ifndef EULERGHOST1DSOLVER_H
#define EULERGHOST1DSOLVER_H

#include "euler1D.H"
#include "method.H"
#include "EoS.H"
#include "material.H"
#include <iostream>
#include <fstream>
#include <sstream>

void solveEuler(std::string scheme, std::string initial, std::string geometry);
std::ostringstream generate_filename(std::string scheme, std::string initial, std::string &geometry, EoS_ptr &EoS_1, EoS_ptr &EoS_2);
void writeToFile(std::ofstream &output, double t, VectOfDouble xCells, VectofVectDouble &u, VectofVectDouble &u2,
                 VectOfDouble &phi, VectofVectDouble &w_ExactRiemann, EoS_ptr &EoS_1, EoS_ptr &EoS_2);
void plot(std::ostringstream &filename, std::string initial);
void ConstructFVDomain(double &dx, VectOfDouble &xCells);
void setInitialConditions(const VectOfDouble &xCells, VectofVectDouble &u, VectofVectDouble &u2, VectOfDouble &phi, double &tStop,
                          VectofVectDouble &w_ExactRiemann, std::string &initial, std::string &scheme, bool exactSolution,
                          EoS_ptr &EoS_1, EoS_ptr &EoS_2);
void computeDt(VectofVectDouble &u, double &dx, double &dt, int &step, EoS_ptr &EoS);
void applyBoundaryConditions(VectofVectDouble &u, std::array<std::string, nVar> boundary, EoS_ptr &EoS);
void applyGhostBoundaryConditions(VectofVectDouble &u, VectofVectDouble &u2, VectOfDouble &phi, EoS_ptr &EoS_1, EoS_ptr &EoS_2);

#endif