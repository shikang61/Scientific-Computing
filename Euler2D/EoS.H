#ifndef EOS_H
#define EOS_H

#include <sstream>
#include <array>
#include <vector>
#include <array>
#include <cmath>
#include <math.h>
#include <algorithm>
#include <memory>
#include <iostream>
#include <fstream>
#include <sstream>

// Forward declaration
class EquationOfState;
class NumericalMethod;
class ODE_solver;
class Material;

typedef std::vector<double> VectOfDouble;
typedef std::vector<std::vector<double>> Vect1D;
typedef std::vector<std::vector<std::vector<double>>> Vect2D;
typedef std::vector<std::vector<std::vector<std::vector<double>>>> Vect3D;

typedef std::shared_ptr<EquationOfState> EoS_ptr;
typedef std::shared_ptr<NumericalMethod> method_ptr;
typedef std::shared_ptr<ODE_solver> ODE_ptr;
typedef std::shared_ptr<Material> Mat_ptr;

const int nVar = 4;
enum CINDEX
{
    RHO = 0,
    MOMX = 1,
    MOMY = 2,
    ENE = 3
};
enum PINDEX
{
    VX = 1,
    VY = 2,
    P = 3
};

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
};

class JWLEoS : public EquationOfState
{
private:
    std::string m_name = "JWL";
    double JWL_rho0 = 1840;
    double JWL_nu0 = 1 / JWL_rho0;
    double JWL_A = 854.5e9;
    double JWL_B = 20.5e9;
    double JWL_R1 = 4.6;
    double JWL_R2 = 1.35;

public:
    JWLEoS(double);      // Constructor
    ~JWLEoS() = default; // Destructor

    std::string info() const override;

    VectOfDouble conservativeToPrimitive(const VectOfDouble &u) const override;
    VectOfDouble primitiveToConservative(const VectOfDouble &w) const override;
    double JWL_pref(const double &rho) const;
    double JWL_eref(const double &rho) const;
    double dJWL_pref(const double &rho) const;
    double dJWL_eref(const double &rho) const;

    double c_s(const VectOfDouble &w) const override;
};

#endif