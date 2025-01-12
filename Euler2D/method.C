#include "method.H"

void createNumericalMethod(method_ptr &numMethod, const std::string &scheme, const std::string &update)
{
    if (scheme == "FORCE")
    {
        numMethod = std::make_shared<FORCE>(update);
    }
    else if (scheme == "SLIC")
    {
        numMethod = std::make_shared<SLIC>(update);
    }
    else
    {
        numMethod = std::make_shared<NumericalMethod>();
        std::cout << "Numerical method not implemented" << std::endl;
    }
}

/* -------------------------------------------------------------------*/
/*                                                                    */
/*                          Base class                                */
/*                                                                    */
/* -------------------------------------------------------------------*/
NumericalMethod::NumericalMethod() {}
NumericalMethod::NumericalMethod(const std::string &update) : update(update) {}
std::string NumericalMethod::name() const
{
    return m_name;
}

int NumericalMethod::ng() const
{
    return m_ng;
};

VectOfDouble NumericalMethod::calc_F(const VectOfDouble &u, const VectOfDouble &w)
{
    if (update == "x")
    {
        // return {u[MOMX], u[MOMX]*w[VX] + w[P], w[MOMX]*w[VY], (u[ENE] + w[P]) * w[VX]}; WTF THIS BUG TOOK 3 DAYS TO FIND
        return {u[MOMX], u[MOMX] * w[VX] + w[P], u[MOMX] * w[VY], (u[ENE] + w[P]) * w[VX]};
    }
    else if (update == "y")
    {
        return {u[MOMY], u[MOMY] * w[VX], u[MOMY] * w[VY] + w[P], (u[ENE] + w[P]) * w[VY]};
    };
    return VectOfDouble(nVar, 0.0);
}

VectOfDouble NumericalMethod::computeFlux(VectOfDouble &u_L, VectOfDouble &u_R, const double &dx, const double &dt, EoS_ptr &EoS)
{
    return VectOfDouble(nVar, 0.0);
}

void NumericalMethod::getFullFlux(Vect2D &F, Vect2D &u, const double &dx, const double &dy, const double &dt, EoS_ptr &EoS)
{
    int Nx = u.size() - 2 * m_ng;
    int Ny = u[0].size() - 2 * m_ng;
    for (int i = 0; i < Nx + 1; i++)
    {
        for (int j = 0; j < Ny + 1; j++)
        {
            if (update == "x")
            {
                F[i][j] = computeFlux(u[i][j + 1], u[i + 1][j + 1], dx, dt, EoS);
            }
            else if (update == "y")
            {
                F[i][j] = computeFlux(u[i + 1][j], u[i + 1][j + 1], dy, dt, EoS);
            };
        };
    };
}

void NumericalMethod::UpdatewithFluxes(Vect2D &u, Vect2D &uPlus1, const double &dx, const double &dy, const double &dt, EoS_ptr &EoS)
{
    int Nx = u.size() - 2 * m_ng;
    int Ny = u[0].size() - 2 * m_ng;
    Vect2D F(Nx + 1, Vect1D(Ny + 1, VectOfDouble(nVar)));
    getFullFlux(F, u, dx, dy, dt, EoS);

    for (int i = 1; i < Nx + 1; i++)
    {
        for (int j = 1; j < Ny + 1; j++)
        {
            for (int n = 0; n < nVar; n++)
            {
                if (update == "x")
                {
                    uPlus1[i + 1][j + 1][n] = u[i + 1][j + 1][n] - (dt / dx) * (F[i][j - 1][n] - F[i - 1][j - 1][n]);
                }
                else if (update == "y")
                {
                    uPlus1[i + 1][j + 1][n] = u[i + 1][j + 1][n] - (dt / dy) * (F[i - 1][j][n] - F[i - 1][j - 1][n]);
                };
            };
        };
    };
};

/* -------------------------------------------------------------------*/
/*                                                                    */
/*                          FORCE class                               */
/*                                                                    */
/* -------------------------------------------------------------------*/

FORCE::FORCE() {};
FORCE::FORCE(const std::string &update) : NumericalMethod(update)
{
    m_name = "FORCE";
    m_ng = 1;
}

VectOfDouble FORCE::LF_Flux(VectOfDouble &u_L, VectOfDouble &u_R, const double &dx, const double &dt, EoS_ptr &EoS)
{
    VectOfDouble w_L = EoS->conservativeToPrimitive(u_L);
    VectOfDouble w_R = EoS->conservativeToPrimitive(u_R);
    VectOfDouble F_LF(nVar);
    for (int i = 0; i < nVar; i++)
    {
        F_LF[i] = 0.5 * dx / dt * (u_L[i] - u_R[i]) + 0.5 * (calc_F(u_L, w_L)[i] + calc_F(u_R, w_R)[i]);
    }
    return F_LF;
}

VectOfDouble FORCE::RI_Flux(VectOfDouble &u_L, VectOfDouble &u_R, const double &dx, const double &dt, EoS_ptr &EoS)
{
    VectOfDouble w_L = EoS->conservativeToPrimitive(u_L);
    VectOfDouble w_R = EoS->conservativeToPrimitive(u_R);
    VectOfDouble F_RI(nVar);
    VectOfDouble u_RI(nVar);
    for (int i = 0; i < nVar; i++)
    {
        u_RI[i] = 0.5 * (u_L[i] + u_R[i]) - 0.5 * dt / dx * (calc_F(u_R, w_R)[i] - calc_F(u_L, w_L)[i]);
    }
    VectOfDouble w_RI = EoS->conservativeToPrimitive(u_RI);
    F_RI = calc_F(u_RI, w_RI);
    return F_RI;
}

VectOfDouble FORCE::computeFlux(VectOfDouble &u_L, VectOfDouble &u_R, const double &dx, const double &dt, EoS_ptr &EoS)
{
    VectOfDouble F_LF(nVar);
    VectOfDouble F_RI(nVar);
    VectOfDouble F_FORCE(nVar);
    F_LF = LF_Flux(u_L, u_R, dx, dt, EoS);
    F_RI = RI_Flux(u_L, u_R, dx, dt, EoS);
    for (int i = 0; i < nVar; i++)
    {
        F_FORCE[i] = 0.5 * (F_LF[i] + F_RI[i]);
    }
    return F_FORCE;
}

/* -------------------------------------------------------------------*/
/*                                                                    */
/*                              SLIC Class                            */
/*                                                                    */
/* -------------------------------------------------------------------*/
// SLIC is data reconstruction, half time-step update and and compute flux with FORCE method

SLIC::SLIC() {};
SLIC::SLIC(const std::string &update) : FORCE(update)
{
    m_name = "SLIC";
    m_ng = 2;
    limiter = "Minbee";
}

VectOfDouble SLIC::slope_limiter(VectOfDouble &u_L, VectOfDouble &u_C, VectOfDouble &u_R, std::string limiter)
{
    VectOfDouble xi(nVar);
    if (limiter == "Minbee")
    {
        // Limit on every variable
        for (int i = 0; i < nVar; i++)
        {
            if ((u_R[i] - u_C[i]) == 0)
            {
                xi[i] = 0;
            }
            else
            {
                double r = (u_C[i] - u_L[i]) / (u_R[i] - u_C[i]);
                double xi_R = 2 / (1 + r);
                if (r <= 0)
                {
                    xi[i] = 0;
                }
                else if ((r <= 1) && (r > 0))
                {
                    xi[i] = r;
                }
                else
                {
                    xi[i] = std::min(r, xi_R);
                };
            };
        };
    }
    return xi;
}

VectOfDouble SLIC::data_reconstruction(VectOfDouble &u_L, VectOfDouble &u_C, VectOfDouble &u_R, std::string side, std::string limiter)
{
    VectOfDouble uBar(nVar);
    VectOfDouble xi = slope_limiter(u_L, u_C, u_R, limiter);
    double omega = 0.0;
    VectOfDouble Delta(nVar);
    for (int i = 0; i < nVar; i++)
    {
        Delta[i] = 0.5 * (1 + omega) * (u_C[i] - u_L[i]) + 0.5 * (1 - omega) * (u_R[i] - u_C[i]);
        if (side == "L")
        {
            uBar[i] = u_C[i] - 0.5 * xi[i] * Delta[i];
        }
        else if (side == "R")
        {
            uBar[i] = u_C[i] + 0.5 * xi[i] * Delta[i];
        }
    }
    return uBar;
}

void SLIC::getFullFlux(Vect2D &F, Vect2D &u, const double &dx, const double &dy, const double &dt, EoS_ptr &EoS)
{
    int Nx = u.size() - 2 * m_ng;
    int Ny = u[0].size() - 2 * m_ng;

    Vect2D uBarL(Nx + m_ng, Vect1D(Ny + m_ng, VectOfDouble(nVar)));
    Vect2D uBarR(Nx + m_ng, Vect1D(Ny + m_ng, VectOfDouble(nVar)));
    Vect2D uBarLUpdate(Nx + m_ng, Vect1D(Ny + m_ng, VectOfDouble(nVar)));
    Vect2D uBarRUpdate(Nx + m_ng, Vect1D(Ny + m_ng, VectOfDouble(nVar)));

    for (int i = 0; i < Nx + m_ng; i++)
    {
        for (int j = 0; j < Ny + m_ng; j++)
        {
            if (update == "x")
            {
                uBarL[i][j] = data_reconstruction(u[i][j + 1], u[i + 1][j + 1], u[i + 2][j + 1], "L", limiter);
                uBarR[i][j] = data_reconstruction(u[i][j + 1], u[i + 1][j + 1], u[i + 2][j + 1], "R", limiter);
            }
            else if (update == "y")
            {
                uBarL[i][j] = data_reconstruction(u[i + 1][j], u[i + 1][j + 1], u[i + 1][j + 2], "L", limiter);
                uBarR[i][j] = data_reconstruction(u[i + 1][j], u[i + 1][j + 1], u[i + 1][j + 2], "R", limiter);
            };
        };
    };
    for (int i = 0; i < Nx + m_ng; i++)
    {
        for (int j = 0; j < Ny + m_ng; j++)
        {
            VectOfDouble wBarL = EoS->conservativeToPrimitive(uBarL[i][j]);
            VectOfDouble wBarR = EoS->conservativeToPrimitive(uBarR[i][j]);
            for (int n = 0; n < nVar; n++)
            {
                if (update == "x")
                {
                    uBarLUpdate[i][j][n] = uBarL[i][j][n] - 0.5 * dt / dx * (calc_F(uBarR[i][j], wBarR)[n] - calc_F(uBarL[i][j], wBarL)[n]);
                    uBarRUpdate[i][j][n] = uBarR[i][j][n] - 0.5 * dt / dx * (calc_F(uBarR[i][j], wBarR)[n] - calc_F(uBarL[i][j], wBarL)[n]);
                }
                else if (update == "y")
                {
                    uBarLUpdate[i][j][n] = uBarL[i][j][n] - 0.5 * dt / dy * (calc_F(uBarR[i][j], wBarR)[n] - calc_F(uBarL[i][j], wBarL)[n]);
                    uBarRUpdate[i][j][n] = uBarR[i][j][n] - 0.5 * dt / dy * (calc_F(uBarR[i][j], wBarR)[n] - calc_F(uBarL[i][j], wBarL)[n]);
                };
            };
        };
    };

    for (int i = 0; i < Nx + 1; i++)
    {
        for (int j = 0; j < Ny + 1; j++)
        {
            if (update == "x")
            {
                F[i][j] = computeFlux(uBarRUpdate[i][j + 1], uBarLUpdate[i + 1][j + 1], dx, dt, EoS);
            }
            else if (update == "y")
            {
                F[i][j] = computeFlux(uBarRUpdate[i + 1][j], uBarLUpdate[i + 1][j + 1], dy, dt, EoS);
            };
        };
    };
}

/* -------------------------------------------------------------------*/
/*                                                                    */
/*                    ODE Solver (Source term)                        */
/*                                                                    */
/* -------------------------------------------------------------------*/

ODE_solver::ODE_solver(const std::string &geometry, int &ng) : m_geometry(geometry),
                                                               m_ng(ng)
{
    m_method = "heun"; // default method, write derived class for other high order methods
}

std::string ODE_solver::method() const
{
    return m_method;
}

VectOfDouble ODE_solver::sourceTerm(const double &x, VectOfDouble &u, EoS_ptr &EoS)
{
    VectOfDouble w(nVar), S(nVar);
    w = EoS->conservativeToPrimitive(u);
    double rho = w[RHO];
    double vx = w[VX];
    double vy = w[VY];
    double p = w[P];
    double E = u[ENE];
    S[0] = -1 * rho * vx / x;
    S[1] = -1 * rho * vx * vx / x;
    S[2] = -1 * rho * vx * vy / x;
    S[3] = -1 * vx * (E + p) / x;
    return S;
};

VectOfDouble ODE_solver::ODEUpdate(VectOfDouble &u, const double &x, const double &dt, EoS_ptr &EoS)
{
    VectOfDouble uBar(nVar), source(nVar);
    if (m_geometry == "spherical")
    {
        source = sourceTerm(x, u, EoS);
        if (m_method == "heun")
        {
            VectOfDouble k_1(nVar), k_2(nVar), temp(nVar);
            for (int n = 0; n < nVar; n++)
            {
                k_1[n] = dt * source[n];
            };
            for (int n = 0; n < nVar; n++)
            {
                temp[n] = u[n] + k_1[n];
            }
            for (int n = 0; n < nVar; n++)
            {
                k_2[n] = dt * sourceTerm(x, temp, EoS)[n];
            }
            for (int n = 0; n < nVar; n++)
            {
                uBar[n] = u[n] + 0.5 * (k_1[n] + k_2[n]);
            }
        };
    }
    else
    {
        return u;
    }

    return uBar;
};

void ODE_solver::sourceUpdate(Vect2D &uOld, Vect2D &uNew, Vect2D &xyCells, const double &dt, EoS_ptr &EoS)
{
    int Nx = uOld.size() - 2 * m_ng;
    int Ny = uOld[0].size() - 2 * m_ng;
    for (int i = 0; i < Nx; i++)
    {
        for (int j = 0; j < Ny; j++)
        {
            int posx = i + m_ng;
            int posy = j + m_ng;
            double x = xyCells[posx][posy][0];
            uNew[posx][posy] = ODEUpdate(uOld[posx][posy], x, dt, EoS);
        };
    };
};