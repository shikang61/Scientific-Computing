#ifndef EOS_H
#define EOS_H

#include <sstream>
#include <array>
#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <cmath>
#include <math.h>
#include <algorithm>

// Forward declaration
class EquationOfState;
class NumericalMethod;
class ODE_solver;
class Material;

typedef std::vector<double> VectOfDouble;
typedef std::vector<std::vector<double>> VectofVectDouble;
typedef std::shared_ptr<EquationOfState> EoS_ptr;
typedef std::shared_ptr<NumericalMethod> method_ptr;
typedef std::shared_ptr<ODE_solver> ODE_ptr;
typedef std::shared_ptr<Material> Mat_ptr;

const int nVar = 8;
enum PINDEX
{
    RHO = 0,
    VX,
    VY,
    VZ,
    P,
    BX,
    BY,
    BZ
};
enum CINDEX
{
    MOMX = 1,
    MOMY,
    MOMZ,
    ENE
};

const double atm = 101325.0;

/*-------------------------------------------------------------------------------------------------------------------------------------------*/

void createEoS(EoS_ptr &EoS, const std::string &type, const double &gamma, const double &p_ref, const double &e_ref = 0.0);

class EquationOfState
{
protected:
    double m_gamma;
    double m_p_ref;
    std::string m_name;

public:
    EquationOfState(double, double); // initialisation constructon
    EquationOfState(double);
    EquationOfState(); // default constructor
    virtual ~EquationOfState() = default;

    std::string name() const;
    double gamma() const;
    double p_ref() const;

    virtual std::string info() const = 0;
    virtual VectOfDouble conservativeToPrimitive(const VectOfDouble &u) const = 0;
    virtual VectOfDouble primitiveToConservative(const VectOfDouble &w) const = 0;
    virtual double c_s(const VectOfDouble &w) const = 0;
    double c_a(const VectOfDouble &w) const;
    double c_f(const VectOfDouble &w) const;
    double c_sl(const VectOfDouble &w) const;
    virtual double enthalpy(const VectOfDouble &u, const VectOfDouble &w) const;
};

// maybe inherit from Mie-gruneisen EoS is better
class IdealGasEoS : public EquationOfState
{
private:
    std::string m_name = "ideal";

public:
    IdealGasEoS(double);      // Constructor
    ~IdealGasEoS() = default; // Destructor

    std::string info() const override;

    VectOfDouble conservativeToPrimitive(const VectOfDouble &u) const override;
    VectOfDouble primitiveToConservative(const VectOfDouble &w) const override;
    double c_s(const VectOfDouble &w) const override;
};

class StiffenedGasEoS : public EquationOfState
{
private:
    std::string m_name = "stiffened";

public:
    StiffenedGasEoS(double, double); // Constructor
    ~StiffenedGasEoS() = default;    // Destructor

    std::string info() const override;

    VectOfDouble conservativeToPrimitive(const VectOfDouble &u) const override;
    VectOfDouble primitiveToConservative(const VectOfDouble &w) const override;
    double c_s(const VectOfDouble &w) const override;
    double enthalpy(const VectOfDouble &u, const VectOfDouble &w) const override;
};

class TabulatedEoS : public EquationOfState
{
private:
    std::string m_name = "tabulated";
    VectOfDouble densities;
    VectOfDouble pressures;
    VectofVectDouble soundSpeeds;
    VectofVectDouble specificEnergies;

public:
    TabulatedEoS();            // Constructor
    ~TabulatedEoS() = default; // Destructor

    std::string info() const override;

    VectOfDouble conservativeToPrimitive(const VectOfDouble &u) const override;
    VectOfDouble primitiveToConservative(const VectOfDouble &w) const override;
    double c_s(const VectOfDouble &w) const override;
    void readData(const std::string &filename);
    double search_pressure_index(const double &p) const;
    double search_density_index(const double &rho) const;
    double bilinear_interp(const double &frac_p_index, const double &frac_rho_index, const VectofVectDouble &table) const;
    std::tuple<double, double> guess_pressure_from_rho_and_E(const double &frac_rho_index, const double &target, const VectofVectDouble &table) const;
};
#endif