#include "method.H"

void createNumericalMethod(method_ptr &numMethod, std::string scheme)
{
    if (scheme == "FORCE")
    {
        numMethod = std::make_shared<FORCE>();
    }
    else if (scheme == "Godunov")
    {
        numMethod = std::make_shared<Godunov>();
    }
    else if (scheme == "SLIC")
    {
        numMethod = std::make_shared<SLIC>();
    }
    else if (scheme == "HLL")
    {
        numMethod = std::make_shared<HLL>();
    }
    else if (scheme == "HLLC")
    {
        numMethod = std::make_shared<HLLC>();
    }
    else if (scheme == "FVTV")
    {
        numMethod = std::make_shared<FVTV>();
    }
    else if (scheme == "MUSCL-Hancock")
    {
        std::cout << "To be implemented" << std::endl;
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
std::string NumericalMethod::name() const
{
    return m_name;
}

int NumericalMethod::ng() const
{
    return m_ng;
};

VectOfDouble NumericalMethod::calc_F(VectOfDouble &u, VectOfDouble &w)
{
    VectOfDouble F(nVar);
    F[0] = u[MOMX];
    F[1] = u[MOMX] * w[VX] + w[P] + 0.5 * (w[BY] * w[BY] + w[BZ] * w[BZ] + w[BX] * w[BX]) - w[BX] * w[BX];
    F[2] = u[MOMX] * w[VY] - w[BX] * w[BY];
    F[3] = u[MOMX] * w[VZ] - w[BX] * w[BZ];
    F[4] = (u[ENE] + w[P] + 0.5 * (w[BY] * w[BY] + w[BZ] * w[BZ] + w[BX] * w[BX])) * w[VX] - (w[VX] * w[BX] + w[VY] * w[BY] + w[VZ] * w[BZ]) * w[BX];
    F[5] = 0;
    F[6] = w[BY] * w[VX] - w[BX] * w[VY];
    F[7] = w[BZ] * w[VX] - w[BX] * w[VZ];
    return F;
}

VectOfDouble NumericalMethod::computeFlux(VectOfDouble &u_L, VectOfDouble &u_R, const double &dx, const double &dt, EoS_ptr &EoS)
{
    return VectOfDouble(nVar, 0.0);
}

// define full flux
void NumericalMethod::getFullFlux(VectofVectDouble &F, VectofVectDouble &u, const double &dx, const double &dt, EoS_ptr &EoS)
{
    for (int i = 0; i < F.size(); i++)
    {
        F[i] = computeFlux(u[i], u[i + 1], dx, dt, EoS);
    };
}

void NumericalMethod::UpdatewithFluxes(VectofVectDouble &u, VectofVectDouble &uPlus1, const double &dx, const double &dt, EoS_ptr &EoS)
{
    int nCells = u.size() - 2 * m_ng;
    VectofVectDouble F(nCells + 1, VectOfDouble(nVar));
    getFullFlux(F, u, dx, dt, EoS);

    for (int i = 1; i < nCells + 1; i++)
    {
        int index = i + m_ng - 1;
        for (int n = 0; n < nVar; n++)
        {
            uPlus1[index][n] = u[index][n] - (dt / dx) * (F[i][n] - F[i - 1][n]);
        };
    };
};

/* -------------------------------------------------------------------*/
/*                                                                    */
/*                          FORCE class                               */
/*                                                                    */
/* -------------------------------------------------------------------*/

FORCE::FORCE()
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
/*                          Godunov class                             */
/*                                                                    */
/* -------------------------------------------------------------------*/

Godunov::Godunov()
{
    m_name = "Godunov";
    m_ng = 1;
}

VectOfDouble Godunov::computeFlux(VectOfDouble &u_L, VectOfDouble &u_R, const double &dx, const double &dt, EoS_ptr &EoS)
{
    VectOfDouble w_state(nVar), u_state(nVar), w_L(nVar), w_R(nVar);
    w_L = EoS->conservativeToPrimitive(u_L);
    w_R = EoS->conservativeToPrimitive(u_R);
    w_state = exactRiemannSolver(0.0, 0.0, 1.0, w_L, w_R, EoS);
    u_state = EoS->primitiveToConservative(w_state);
    return calc_F(u_state, w_state);
}

/* -------------------------------------------------------------------*/
/*                                                                    */
/*                              SLIC Class                            */
/*                                                                    */
/* -------------------------------------------------------------------*/
// SLIC is data reconstruction, half time-step update and and compute flux with FORCE method

SLIC::SLIC()
{
    m_name = "SLIC";
    m_ng = 2;
    limiter = "Minbee";
}

VectOfDouble SLIC::slope_limiter(VectOfDouble &u_L, VectOfDouble &u_C, VectOfDouble &u_R, std::string limiter)
{
    double r, xi_R, xi_L;
    VectOfDouble xi(nVar);
    for (int i = 0; i < nVar; i++)
    {
        if (limiter == "Minbee")
        {
            // if ((u_R[2] - u_C[2]) == 0){
            //     xi = {0.0, 0.0, 0.0};
            // } else {
            //     double r = (u_C[2] - u_L[2]) / (u_R[2] - u_C[2]);
            //     double xi_R = 2/(1+r);
            //     if (r <= 0){
            //             xi = {0.0, 0.0, 0.0};
            //         } else if ((r <= 1) && (r>0)){
            //             xi = {r, r, r};
            //         } else{
            //             xi = {std::min(r, xi_R), std::min(r, xi_R), std::min(r, xi_R)};
            //         };
            // }
            if ((u_R[i] - u_C[i]) == 0)
            {
                xi[i] = 0.0;
            }
            else
            {
                r = (u_C[i] - u_L[i]) / (u_R[i] - u_C[i]);
                xi_R = 2.0 / (1.0 + r);
                xi_L = 2.0 * r / (2 + r);
                if (r <= 0)
                {
                    xi[i] = 0.0;
                }
                else if ((r <= 1.0) && (r > 0.0))
                {
                    xi[i] = r;
                }
                else
                {
                    xi[i] = std::min(r, xi_R);
                };
            };
        }
        else if (limiter == "Vanleer")
        {
            if ((u_R[i] - u_C[i]) == 0)
            {
                xi[i] = 0.0;
            }
            else
            {
                r = (u_C[i] - u_L[i]) / (u_R[i] - u_C[i]);
                xi_R = 2.0 / (1.0 + r);
                xi_L = 2.0 * r / (2 + r);
                if (r <= 0)
                {
                    xi[i] = 0.0;
                }
                else
                {
                    xi[i] = std::min(xi_L, xi_R);
                }
            }
        }
        else if (limiter == "VanAlbada")
        {
            if ((u_R[i] - u_C[i]) == 0)
            {
                xi[i] = 0.0;
            }
            else
            {
                r = (u_C[i] - u_L[i]) / (u_R[i] - u_C[i]);
                xi_R = 2.0 / (1.0 + r);
                xi_L = 2.0 * r / (2 + r);
                if (r <= 0)
                {
                    xi[i] = 0.0;
                }
                else
                {
                    xi[i] = std::min((r * (1 + r)) / (1 + r * r), xi_R);
                }
            }
        }
        else if (limiter == "Superbee")
        {
            if ((u_R[i] - u_C[i]) == 0)
            {
                xi[i] = 0.0;
            }
            else
            {
                r = (u_C[i] - u_L[i]) / (u_R[i] - u_C[i]);
                xi_R = 2.0 / (1.0 + r);
                xi_L = 2.0 * r / (2 + r);
                if (r <= 0)
                {
                    xi[i] = 0.0;
                }
                else if (r > 0 && r <= 0.5)
                {
                    xi[i] = 2 * r;
                }
                else if (r > 0.5 && r <= 1)
                {
                    xi[i] = 1;
                }
                else
                {
                    xi[i] = std::min({r, xi_R, 2.0});
                }
            }
        }
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

void SLIC::getFullFlux(VectofVectDouble &F, VectofVectDouble &u, const double &dx, const double &dt, EoS_ptr &EoS)
{
    int size = u.size() - m_ng;
    VectofVectDouble uBarLUpdate(size, VectOfDouble(nVar));
    VectofVectDouble uBarRUpdate(size, VectOfDouble(nVar));
    VectofVectDouble uBarL(size, VectOfDouble(nVar));
    VectofVectDouble uBarR(size, VectOfDouble(nVar));

    for (int i = 0; i < size; i++)
    {
        uBarL[i] = data_reconstruction(u[i], u[i + 1], u[i + 2], "L", limiter);
        uBarR[i] = data_reconstruction(u[i], u[i + 1], u[i + 2], "R", limiter);
    };

    for (int i = 0; i < size; i++)
    {
        VectOfDouble wBarL = EoS->conservativeToPrimitive(uBarL[i]);
        VectOfDouble wBarR = EoS->conservativeToPrimitive(uBarR[i]);
        for (int j = 0; j < nVar; j++)
        {
            uBarLUpdate[i][j] = uBarL[i][j] - 0.5 * dt / dx * (calc_F(uBarR[i], wBarR)[j] - calc_F(uBarL[i], wBarL)[j]);
            uBarRUpdate[i][j] = uBarR[i][j] - 0.5 * dt / dx * (calc_F(uBarR[i], wBarR)[j] - calc_F(uBarL[i], wBarL)[j]);
        };
    };

    for (int i = 0; i < F.size(); i++)
    {
        F[i] = computeFlux(uBarRUpdate[i], uBarLUpdate[i + 1], dx, dt, EoS);
    };
}

/* -------------------------------------------------------------------*/
/*                                                                    */
/*                          HLL Flux                                  */
/*                                                                    */
/* -------------------------------------------------------------------*/

HLL::HLL()
{
    m_name = "HLL";
    m_ng = 1;
}

VectOfDouble HLL::computeFlux(VectOfDouble &u_L, VectOfDouble &u_R, const double &dx, const double &dt, EoS_ptr &EoS)
{
    double S_L, S_R;
    VectOfDouble w_L(nVar), w_R(nVar), F_HLL(nVar), F_L(nVar), F_R(nVar);
    w_L = EoS->conservativeToPrimitive(u_L);
    w_R = EoS->conservativeToPrimitive(u_R);
    std::tie(S_L, S_R) = wave_speed_estimate(u_L, u_R, EoS);
    // S_plus = std::max(abs(w_L[1])+c_s(gamma, w_L), abs(w_R[1])+c_s(gamma, w_R));
    // S_L = -S_plus;
    // S_R = S_plus;
    F_L = calc_F(u_L, w_L);
    F_R = calc_F(u_R, w_R);

    if (S_L >= 0)
    {
        return F_L;
    }
    else if (S_R <= 0)
    {
        return F_R;
    }
    for (int i = 0; i < nVar; i++)
    {
        F_HLL[i] = (S_R * F_L[i] - S_L * F_R[i] + S_L * S_R * (u_R[i] - u_L[i])) / (S_R - S_L);
    };
    return F_HLL;
};

/* -------------------------------------------------------------------*/
/*                                                                    */
/*                          HLLC Flux                                 */
/*                                                                    */
/* -------------------------------------------------------------------*/

HLLC::HLLC()
{
    m_name = "HLLC";
    m_ng = 1;
}
VectOfDouble HLLC::calc_u_HLLC(VectOfDouble &u, double &S, double &Sstar, EoS_ptr &EoS)
{
    VectOfDouble w(nVar), u_HLLC(nVar);
    w = EoS->conservativeToPrimitive(u);
    double E = u[ENE];
    double rho = w[RHO];
    double v = w[VX];
    double p = w[P];

    u_HLLC[RHO] = rho * (S - v) / (S - Sstar);
    u_HLLC[MOMX] = u_HLLC[RHO] * Sstar;
    u_HLLC[ENE] = u_HLLC[RHO] * (E / rho + (Sstar - v) * (Sstar + p / (rho * (S - v))));
    return u_HLLC;
}

VectOfDouble HLLC::computeFlux(VectOfDouble &u_L, VectOfDouble &u_R, const double &dx, const double &dt, EoS_ptr &EoS)
{
    double S_L, Sstar, S_R;
    VectOfDouble w_L(nVar), w_R(nVar), F_HLLC(nVar), F_L(nVar), F_R(nVar), u_L_HLLC(nVar), u_R_HLLC(nVar);
    w_L = EoS->conservativeToPrimitive(u_L);
    w_R = EoS->conservativeToPrimitive(u_R);
    std::tie(S_L, S_R) = wave_speed_estimate(u_L, u_R, EoS);
    F_L = calc_F(u_L, w_L);
    F_R = calc_F(u_R, w_R);
    double rhol = w_L[RHO];
    double rhor = w_R[RHO];
    double vl = w_L[VX];
    double vr = w_R[VX];
    double pl = w_L[P];
    double pr = w_R[P];

    Sstar = (pr - pl + rhol * vl * (S_L - vl) - rhor * vr * (S_R - vr)) / (rhol * (S_L - vl) - rhor * (S_R - vr));
    u_L_HLLC = calc_u_HLLC(u_L, S_L, Sstar, EoS);
    u_R_HLLC = calc_u_HLLC(u_R, S_R, Sstar, EoS);

    if (S_L >= 0)
    {
        return F_L;
    }
    if (S_R < 0)
    {
        return F_R;
    }
    if ((S_L < 0) && (Sstar >= 0))
    {
        for (int i = 0; i < nVar; i++)
        {
            F_HLLC[i] = F_L[i] + S_L * (u_L_HLLC[i] - u_L[i]);
        }
    }
    if ((Sstar < 0) && (S_R >= 0))
    {
        for (int i = 0; i < nVar; i++)
        {
            F_HLLC[i] = F_R[i] + S_R * (u_R_HLLC[i] - u_R[i]);
        }
    }
    return F_HLLC;
};

/* -------------------------------------------------------------------*/
/*                                                                    */
/*                    Flux vector splitting Flux                      */
/*                        (Toro - Vazquez)                            */ // OK
/*                                                                    */
/* -------------------------------------------------------------------*/
FVTV::FVTV()
{
    m_name = "FVTV";
    m_ng = 1;
}
double FVTV::TV_A_ideal(VectOfDouble &w, EoS_ptr &EoS)
{
    double c_s = EoS->c_s(w);
    return sqrt(w[VX] * w[VX] + 4 * c_s * c_s);
}

double FVTV::TV_C(VectOfDouble &w, double &A, std::string side)
{
    if (side == "L")
    {
        return w[RHO] * (w[VX] - A);
    }
    else
    {
        return w[RHO] * (w[VX] + A);
    }
}

std::pair<double, double> FVTV::TV_uStar_pStar(VectOfDouble &u_L, VectOfDouble &u_R, EoS_ptr &EoS)
{
    double C_L, C_R, A_L, A_R, v_R, v_L, p_L, p_R, uStar, pStar;
    VectOfDouble w_L(nVar), w_R(nVar);
    w_L = EoS->conservativeToPrimitive(u_L);
    w_R = EoS->conservativeToPrimitive(u_R);
    v_L = w_L[VX];
    v_R = w_R[VX];
    p_L = w_L[P];
    p_R = w_R[P];
    A_L = TV_A_ideal(w_L, EoS);
    A_R = TV_A_ideal(w_R, EoS);
    C_L = TV_C(w_L, A_L, "L");
    C_R = TV_C(w_R, A_R, "R");
    uStar = (C_R * v_R - C_L * v_L) / (C_R - C_L) - 2 * (p_R - p_L) / (C_R - C_L);
    pStar = (C_R * p_L - C_L * p_R) / (C_R - C_L) + 0.5 * C_L * C_R * (v_R - v_L) / (C_R - C_L);
    return std::make_pair(uStar, pStar);
}

VectOfDouble FVTV::computeFlux(VectOfDouble &u_L, VectOfDouble &u_R, const double &dx, const double &dt, EoS_ptr &EoS)
{
    double uStar, pStar, rho_K;
    VectOfDouble w_L(nVar), w_R(nVar), A_flux(nVar), P_flux(nVar), FVTV_Flux(nVar);
    w_L = EoS->conservativeToPrimitive(u_L);
    w_R = EoS->conservativeToPrimitive(u_R);
    double rhol = w_L[RHO];
    double rhor = w_R[RHO];
    double vl = w_L[VX];
    double vr = w_R[VX];
    std::tie(uStar, pStar) = TV_uStar_pStar(u_L, u_R, EoS);
    if (uStar >= 0)
    {
        A_flux[0] = uStar * rhol;
        A_flux[1] = uStar * rhol * vl;
        A_flux[2] = uStar * 0.5 * rhol * vl * vl;

        rho_K = rhol;
    }
    else
    {
        A_flux[0] = uStar * rhor;
        A_flux[1] = uStar * rhor * vr;
        A_flux[2] = uStar * 0.5 * rhor * vr * vr;

        rho_K = rhor;
    };

    P_flux[0] = 0;
    P_flux[1] = pStar;
    double ideal_e = pStar / ((EoS->gamma() - 1) * rho_K);
    P_flux[2] = uStar * (rho_K * ideal_e + pStar);
    for (int i = 0; i < nVar; i++)
    {
        FVTV_Flux[i] = A_flux[i] + P_flux[i];
    }
    return FVTV_Flux;
};

/* -------------------------------------------------------------------*/
/*                                                                    */
/*                    ODE Solver (Source term)                        */
/*                                                                    */
/* -------------------------------------------------------------------*/

ODE_solver::ODE_solver(const std::string &geometry, int &ng) : m_geometry(geometry), m_ng(ng)
{
    m_method = "heun"; // default method, write derived class for other high order methods

    if (geometry == "3Dspherical")
    {
        m_alpha = 2;
    }
    else if (geometry == "2Dcylindrical")
    {
        m_alpha = 1;
    }
    else
    {
        m_alpha = 0;
    }
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
    double v = w[VX];
    double p = w[P];
    double E = u[ENE];
    S[0] = -m_alpha * rho * v / x;
    S[1] = -m_alpha * rho * v * v / x;
    S[2] = -m_alpha * v * (E + p) / x;
    return S;
};

VectOfDouble ODE_solver::ODEUpdate(VectOfDouble &u, const double &x, const double &dt, EoS_ptr &EoS)
{
    VectOfDouble uBar(nVar), source(nVar);
    source = sourceTerm(x, u, EoS);

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
    return uBar;
};

void ODE_solver::sourceUpdate(VectofVectDouble &uOld, VectofVectDouble &uNew, VectOfDouble &xCells, const double &dt, EoS_ptr &EoS)
{
    for (int i = m_ng; i < xCells.size() - m_ng; i++)
    {
        double x = xCells[i];
        uNew[i] = ODEUpdate(uOld[i], x, dt, EoS);
    };
};
