#ifndef METHOD_H
#define METHOD_H

#include "EoS.H"

void createNumericalMethod(method_ptr &numMethod, const std::string &scheme, const std::string &update);

/* -------------------------------------------------------------------*/
/*                                                                    */
/*                      Numerical method (PDE)                        */
/*                                                                    */
/* -------------------------------------------------------------------*/
class NumericalMethod
{
protected:
    int m_ng;
    std::string m_name;
    std::string update;

public:
    NumericalMethod();
    NumericalMethod(const std::string &);
    virtual ~NumericalMethod() = default;

    std::string name() const;
    int ng() const;

    void UpdatewithFluxes(Vect2D &u, Vect2D &uPlus1, const double &dx, const double &dy, const double &dt, EoS_ptr &EoS);
    VectOfDouble calc_F(const VectOfDouble &u, const VectOfDouble &w);
    virtual VectOfDouble computeFlux(VectOfDouble &u_L, VectOfDouble &u_R, const double &dx, const double &dt, EoS_ptr &EoS);
    virtual void getFullFlux(Vect2D &F, Vect2D &u, const double &dx, const double &dy, const double &dt, EoS_ptr &EoS);
};

class FORCE : public NumericalMethod
{

public:
    FORCE();
    FORCE(const std::string &);
    ~FORCE() = default;
    VectOfDouble LF_Flux(VectOfDouble &u_L, VectOfDouble &u_R, const double &dx, const double &dt, EoS_ptr &EoS);
    VectOfDouble RI_Flux(VectOfDouble &u_L, VectOfDouble &u_R, const double &dx, const double &dt, EoS_ptr &EoS);
    VectOfDouble computeFlux(VectOfDouble &u_L, VectOfDouble &u_R, const double &dx, const double &dt, EoS_ptr &EoS) override;
};

class SLIC : public FORCE
{
private:
    std::string limiter;

public:
    SLIC();
    SLIC(const std::string &);
    ~SLIC() = default;
    virtual VectOfDouble slope_limiter(VectOfDouble &u_L, VectOfDouble &u_C, VectOfDouble &u_R, std::string limiter);
    virtual VectOfDouble data_reconstruction(VectOfDouble &u_L, VectOfDouble &u_C, VectOfDouble &u_R, std::string side, std::string limiter);
    void getFullFlux(Vect2D &F, Vect2D &u, const double &dx, const double &dy, const double &dt, EoS_ptr &EoS) override;
};

/* -------------------------------------------------------------------*/
/*                                                                    */
/*                    ODE Solver (Source term)                        */
/*                                                                    */
/* -------------------------------------------------------------------*/
class ODE_solver
{
protected:
    std::string m_method;
    std::string m_geometry;
    int m_ng;

public:
    ODE_solver(const std::string &, int &);
    virtual ~ODE_solver() = default;

    std::string method() const;

    VectOfDouble sourceTerm(const double &x, VectOfDouble &u, EoS_ptr &EoS);
    void sourceUpdate(Vect2D &uOld, Vect2D &uNew, Vect2D &xyCells, const double &dt, EoS_ptr &EoS);
    virtual VectOfDouble ODEUpdate(VectOfDouble &u, const double &x, const double &dt, EoS_ptr &EoS);
};

#endif