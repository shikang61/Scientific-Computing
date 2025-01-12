#ifndef EULER2DSOLVER_H
#define EULER2DSOLVER_H

#include "method.H"
#include "EoS.H"
#include "material.H"

void solveEuler(std::string scheme, std::string initial, std::pair<std::string, int> &slice, std::string geometry);
std::ostringstream generate_filename(method_ptr &numMethod, std::string &initial, std::pair<std::string, int> &slice, std::string &geometry, EoS_ptr &EoS);
void writeToFile(std::ofstream &output, double &t, Vect2D &xyCells, Vect2D &u, std::pair<std::string, int> &slice, EoS_ptr &EoS);
void plot(std::ostringstream &filename, std::string &initial, std::pair<std::string, int> &slice);

void ConstructFVDomain(double &dx, double &dy, Vect2D &xyCells);
void setInitialConditions(const Vect2D &xyCells, Vect2D &u, double &tStop, std::string &initial, method_ptr &numMethod, std::string &geometry, EoS_ptr &EoS);
double findMax(const Vect1D &a);
void computeDt(Vect2D &u, double &dx, double &dy, double &dt, int step, EoS_ptr &EoS);

void applyBoundaryConditions(Vect2D &u, std::array<std::string, nVar> &boundary, int &Nx, int &Ny, const std::string &geometry);

#endif