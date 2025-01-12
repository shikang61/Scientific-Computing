#ifndef METHOD_H
#define METHOD_H

#include "EoS.H"
#include "euler1D.H"
/*-------------------------------------------------------------------------------------------------------------------------------------------*/

void createNumericalMethod(Method_ptr &numMethod, std::string scheme);

/*-------------------------------------------------------------------------------------------------------------------------------------------*/

/* -------------------------------------------------------------------*/
/*                                                                    */
/*                 Numerical method base class (PDE)                  */
/*                                                                    */
/* -------------------------------------------------------------------*/
class NumericalMethod
{
protected:
    int m_ng;
    std::string m_name;
    bool m_MHD;

public:
    NumericalMethod();
    virtual ~NumericalMethod() = default;

    virtual std::string name() const;
    virtual int ng() const;

    void UpdatewithFluxes(VectofVectDouble &u, VectofVectDouble &uPlus1, const double &dx, const double &dt, EoS_ptr &EoS);
    virtual VectOfDouble calc_F(VectOfDouble &u, VectOfDouble &w);
    virtual VectOfDouble computeFlux(VectOfDouble &u_L, VectOfDouble &u_R, const double &dx, const double &dt, EoS_ptr &EoS);
    virtual void getFullFlux(VectofVectDouble &F, VectofVectDouble &u, const double &dx, const double &dt, EoS_ptr &EoS);
};

class FORCE : public NumericalMethod
{
public:
    FORCE();
    ~FORCE() = default;
    VectOfDouble LF_Flux(VectOfDouble &u_L, VectOfDouble &u_R, const double &dx, const double &dt, EoS_ptr &EoS);
    VectOfDouble RI_Flux(VectOfDouble &u_L, VectOfDouble &u_R, const double &dx, const double &dt, EoS_ptr &EoS);
    VectOfDouble computeFlux(VectOfDouble &u_L, VectOfDouble &u_R, const double &dx, const double &dt, EoS_ptr &EoS) override;
};

class Godunov : public NumericalMethod
{
public:
    Godunov();
    ~Godunov() = default;
    VectOfDouble computeFlux(VectOfDouble &u_L, VectOfDouble &u_R, const double &dx, const double &dt, EoS_ptr &EoS) override;
};

class SLIC : public FORCE
{
private:
    std::string limiter;

public:
    SLIC();
    ~SLIC() = default;
    virtual VectOfDouble slope_limiter(VectOfDouble &u_L, VectOfDouble &u_C, VectOfDouble &u_R, std::string limiter);
    virtual VectOfDouble data_reconstruction(VectOfDouble &u_L, VectOfDouble &u_C, VectOfDouble &u_R, std::string side, std::string limiter);
    void getFullFlux(VectofVectDouble &F, VectofVectDouble &u, const double &dx, const double &dt, EoS_ptr &EoS) override;
};

class HLL : public NumericalMethod
{
public:
    HLL();
    ~HLL() = default;
    VectOfDouble computeFlux(VectOfDouble &u_L, VectOfDouble &u_R, const double &dx, const double &dt, EoS_ptr &EoS) override;
};

class HLLC : public NumericalMethod
{
public:
    HLLC();
    ~HLLC() = default;
    VectOfDouble calc_u_HLLC(VectOfDouble &u, double &S, double &Sstar, EoS_ptr &EoS);
    VectOfDouble computeFlux(VectOfDouble &u_L, VectOfDouble &u_R, const double &dx, const double &dt, EoS_ptr &EoS) override;
};

class MUSCL_HLLC : public NumericalMethod
{
private:
    HLLC m_HLLC;
    SLIC m_SLIC;
    std::string limiter;

public:
    MUSCL_HLLC();
    ~MUSCL_HLLC() = default;
    void getFullFlux(VectofVectDouble &F, VectofVectDouble &u, const double &dx, const double &dt, EoS_ptr &EoS) override;
    VectOfDouble computeFlux(VectOfDouble &u_L, VectOfDouble &u_R, const double &dx, const double &dt, EoS_ptr &EoS) override;
};

class MUSCL_Godunov : public MUSCL_HLLC
{
private:
    Godunov m_Godunov;
    SLIC m_SLIC;
    std::string limiter;

public:
    MUSCL_Godunov();
    ~MUSCL_Godunov() = default;
    VectOfDouble computeFlux(VectOfDouble &u_L, VectOfDouble &u_R, const double &dx, const double &dt, EoS_ptr &EoS) override;
};

class FVTV : public NumericalMethod
{
public:
    FVTV();
    ~FVTV() = default;
    double TV_A_ideal(VectOfDouble &w, EoS_ptr &EoS);
    double TV_C(VectOfDouble &w, double &A, std::string side);
    std::pair<double, double> TV_uStar_pStar(VectOfDouble &u_L, VectOfDouble &u_R, EoS_ptr &EoS);
    VectOfDouble computeFlux(VectOfDouble &u_L, VectOfDouble &u_R, const double &dx, const double &dt, EoS_ptr &EoS) override;
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
    int m_alpha;
    int m_ng;

public:
    ODE_solver(const std::string &, int &);
    virtual ~ODE_solver() = default;

    std::string method() const;

    VectOfDouble sourceTerm(const double &x, VectOfDouble &u, EoS_ptr &EoS);
    void sourceUpdate(VectofVectDouble &uOld, VectofVectDouble &uNew, VectOfDouble &xCells, const double &dt, EoS_ptr &EoS);
    virtual VectOfDouble ODEUpdate(VectOfDouble &u, const double &x, const double &dt, EoS_ptr &EoS);
};

#endif