#include "euler1D.H"

/* -------------------------------------------------------------------*/
/*                                                                    */
/*                      Exact Riemann Solver                          */
/*                                                                    */
/* -------------------------------------------------------------------*/

double c_s_star(double p_star, VectOfDouble &w, EoS_ptr &EoS)
{
    double gamma = EoS->gamma();
    double c_s = EoS->c_s(w);
    double p_ref = EoS->p_ref();
    return c_s * pow((p_star + p_ref) / (w[P] + p_ref), 0.5 * (gamma - 1) / gamma);
};

double A(VectOfDouble &w, EoS_ptr &EoS)
{
    double gamma = EoS->gamma();
    return 2.0 / ((gamma + 1) * w[RHO]);
};

double B(VectOfDouble &w, EoS_ptr &EoS)
{
    double gamma = EoS->gamma();
    double p_ref = EoS->p_ref();
    return (gamma - 1) / (gamma + 1) * w[P] + 2.0 * gamma * p_ref / (gamma + 1);
};

double f_K(double p, VectOfDouble &w, EoS_ptr &EoS)
{
    double gamma = EoS->gamma();
    double c_s = EoS->c_s(w);
    double p_ref = EoS->p_ref();
    if (p > w[P])
    {
        return (p - w[P]) * sqrt(A(w, EoS) / (p + B(w, EoS)));
    }
    else
    {
        return 2.0 * c_s / (gamma - 1) * (pow((p + p_ref) / (w[P] + p_ref), (gamma - 1) / (2 * gamma)) - 1);
    }
};

double df_K(double p, VectOfDouble &w, EoS_ptr &EoS)
{
    double gamma = EoS->gamma();
    double c_s = EoS->c_s(w);
    double p_ref = EoS->p_ref();
    if (p > w[P])
    {
        return sqrt(A(w, EoS) / (p + B(w, EoS))) * (1 - 0.5 * (p - w[P]) / (p + B(w, EoS)));
    }
    else
    {
        return 1 / (c_s * w[RHO]) * pow((p + p_ref) / (w[P] + p_ref), -1 * (gamma + 1) / (2 * gamma));
    }
};

double f(double p, VectOfDouble &w_L, VectOfDouble &w_R, EoS_ptr &EoS)
{
    return f_K(p, w_L, EoS) + f_K(p, w_R, EoS) + (w_R[VX] - w_L[VX]);
};

double df(double p, VectOfDouble &w_L, VectOfDouble &w_R, EoS_ptr &EoS)
{
    return df_K(p, w_L, EoS) + df_K(p, w_R, EoS);
};

double newtonRaphson(double initialGuess, double tolerance, int maxIterations, VectOfDouble &w_L, VectOfDouble &w_R, EoS_ptr &EoS)
{
    double p_prev = initialGuess; // Start
    double p;
    double f_p;
    double df_p;
    double CHA;
    for (int i = 0; i < maxIterations; ++i)
    {
        f_p = f(p_prev, w_L, w_R, EoS);
        df_p = df(p_prev, w_L, w_R, EoS);
        p = p_prev - f_p / df_p; // Update p using Newton-Raphson formula
        if (p <= 0)
        {
            p = tolerance;
        }
        CHA = abs(abs(p - p_prev) / (0.5 * (p + p_prev)));
        if (CHA <= tolerance)
        {
            return p; // Return the root approximation
        }
        p_prev = p;
    }
    std::cout << "Max iterations reached without convergence.\n";
    return p;
};

VectOfDouble raref(double S, VectOfDouble &w, std::string side, EoS_ptr &EoS)
{
    double rho = 0.0;
    double u = 0.0;
    double p = 0.0;
    double gamma = EoS->gamma();
    double c_s = EoS->c_s(w);
    double p_ref = EoS->p_ref();
    VectOfDouble w_raref(nVar);

    if (side == "L")
    {
        rho = w[RHO] * pow(2 / (gamma + 1) + (gamma - 1) / ((gamma + 1) * c_s) * (w[VX] - S), 2 / (gamma - 1));
        u = 2 / (gamma + 1) * (c_s + 0.5 * (gamma - 1) * w[VX] + S);
        p = (w[P] + p_ref) * pow(2 / (gamma + 1) + (gamma - 1) / ((gamma + 1) * c_s) * (w[VX] - S), 2 * gamma / (gamma - 1));
    }
    else if (side == "R")
    {
        rho = w[RHO] * pow(2 / (gamma + 1) - (gamma - 1) / ((gamma + 1) * c_s) * (w[VX] - S), 2 / (gamma - 1));
        u = 2 / (gamma + 1) * (-1 * c_s + 0.5 * (gamma - 1) * w[VX] + S);
        p = (w[P] + p_ref) * pow(2 / (gamma + 1) - (gamma - 1) / ((gamma + 1) * c_s) * (w[VX] - S), 2 * gamma / (gamma - 1));
    }

    w_raref[RHO] = rho;
    w_raref[VX] = u;
    w_raref[P] = p - p_ref;
    return w_raref;
};

double rho_star_shock(double p_star, VectOfDouble w, EoS_ptr &EoS)
{
    double p_ref = EoS->p_ref();
    double gamma = EoS->gamma();
    return w[RHO] * ((p_star + p_ref) / (w[P] + p_ref) + (gamma - 1) / (gamma + 1)) / ((gamma - 1) / (gamma + 1) * (p_star + p_ref) / (w[P] + p_ref) + 1);
};

double S_shock(double p_star, VectOfDouble w, std::string side, EoS_ptr &EoS)
{
    double gamma = EoS->gamma();
    double c_s = EoS->c_s(w);
    double p_ref = EoS->p_ref();
    if (side == "L")
    {
        return w[VX] - c_s * sqrt((gamma + 1) / (2 * gamma) * (p_star + p_ref) / (w[P] + p_ref) + (gamma - 1) / (2 * gamma));
    }
    else if (side == "R")
    {
        return w[VX] + c_s * sqrt((gamma + 1) / (2 * gamma) * (p_star + p_ref) / (w[P] + p_ref) + (gamma - 1) / (2 * gamma));
    }
    return 0.0;
};

double S_H(VectOfDouble w, std::string side, EoS_ptr &EoS)
{
    double c_s = EoS->c_s(w);
    if (side == "L")
    {
        return (w[VX] - c_s);
    }
    else if (side == "R")
    {
        return (w[VX] + c_s);
    }
    return 0.0;
};

double S_T(VectOfDouble w_star, double p_star, std::string side, EoS_ptr &EoS)
{
    if (side == "L")
    {
        return (w_star[VX] - c_s_star(p_star, w_star, EoS));
    }
    else if (side == "R")
    {
        return (w_star[VX] + c_s_star(p_star, w_star, EoS));
    }
    return 0.0;
};

double rho_star_raref(double p_star, VectOfDouble w, EoS_ptr &EoS)
{
    double gamma = EoS->gamma();
    double p_ref = EoS->p_ref();
    return w[RHO] * pow((p_star + p_ref) / (w[P] + p_ref), 1 / gamma);
};

VectOfDouble exactRiemannSolver(double x, double mid, double tStop, VectOfDouble &w_L_initial, VectOfDouble &w_R_initial, EoS_ptr &EoS)
{
    double initialGuess = 0.5 * (w_L_initial[P] + w_R_initial[P]);
    double p_star = newtonRaphson(initialGuess, 1e-6, 10e6, w_L_initial, w_R_initial, EoS);
    double u_star = 0.5 * (w_L_initial[VX] + w_R_initial[VX]) + 0.5 * (f_K(p_star, w_R_initial, EoS) - f_K(p_star, w_L_initial, EoS));
    double S = (x - mid) / tStop;

    if (S <= u_star)
    {
        VectOfDouble w_LEFT(nVar);
        if (p_star > w_L_initial[P])
        { // Left shock wave
            double rho_L_star = rho_star_shock(p_star, w_L_initial, EoS);
            double S_L = S_shock(p_star, w_L_initial, "L", EoS);
            VectOfDouble w_sho_L = {rho_L_star, u_star, p_star};
            if (S >= S_L)
            {
                w_LEFT = w_sho_L;
            }
            else if (S <= S_L)
            {
                w_LEFT = w_L_initial;
            }
        }
        else if (p_star <= w_L_initial[P])
        { // Left rarefaction wave
            double rho_L_star = rho_star_raref(p_star, w_L_initial, EoS);
            VectOfDouble w_Lfan = raref(S, w_L_initial, "L", EoS);
            VectOfDouble w_Lstar = {rho_L_star, u_star, p_star};
            double S_HL = S_H(w_L_initial, "L", EoS);
            double S_TL = S_T(w_Lstar, p_star, "L", EoS);
            if (S <= S_HL)
            {
                w_LEFT = w_L_initial;
            }
            else if ((S >= S_HL) && (S <= S_TL))
            {
                w_LEFT = w_Lfan;
            }
            else if (S >= S_TL)
            {
                w_LEFT = w_Lstar;
            }
        };
        return w_LEFT;
    }
    else if (S >= u_star)
    {
        VectOfDouble w_RIGHT(nVar);
        if (p_star > w_R_initial[P])
        { // Right shock wave
            double rho_R_star = rho_star_shock(p_star, w_R_initial, EoS);
            double S_R = S_shock(p_star, w_R_initial, "R", EoS);
            VectOfDouble w_sho_R = {rho_R_star, u_star, p_star};
            if (S <= S_R)
            {
                w_RIGHT = w_sho_R;
            }
            else if (S >= S_R)
            {
                w_RIGHT = w_R_initial;
            }
        }
        else if (p_star <= w_R_initial[P])
        { // Right rarefaction wave
            VectOfDouble w_Rfan = raref(S, w_R_initial, "R", EoS);
            VectOfDouble w_Rstar = {rho_star_raref(p_star, w_R_initial, EoS), u_star, p_star};
            double S_HR = S_H(w_R_initial, "R", EoS);
            double S_TR = S_T(w_Rstar, p_star, "R", EoS);
            if (S >= S_HR)
            {
                w_RIGHT = w_R_initial;
            }
            else if ((S <= S_HR) && (S >= S_TR))
            {
                w_RIGHT = w_Rfan;
            }
            else if (S <= S_TR)
            {
                w_RIGHT = w_Rstar;
            }
        }
        return w_RIGHT;
    };
    return VectOfDouble(nVar, 0.0);
};

/* -------------------------------------------------------------------*/
/*                                                                    */
/*                            Additional                              */
/*                                                                    */
/* -------------------------------------------------------------------*/

std::pair<double, double> wave_speed_estimate(VectOfDouble &u_L, VectOfDouble &u_R, EoS_ptr &EoS)
{
    VectOfDouble w_L(nVar), w_R(nVar);
    double S_L, S_R;
    w_L = EoS->conservativeToPrimitive(u_L);
    w_R = EoS->conservativeToPrimitive(u_R);
    double initialGuess = 0.5 * (w_L[P] + w_R[P]);
    double p_star = newtonRaphson(initialGuess, 1e-6, 10e6, w_L, w_R, EoS);

    if (p_star <= w_L[P])
    {
        S_L = S_H(w_L, "L", EoS);
    }
    else
    {
        S_L = S_shock(p_star, w_L, "L", EoS);
    }

    if (p_star <= w_R[P])
    {
        S_R = S_H(w_R, "R", EoS);
    }
    else
    {
        S_R = S_shock(p_star, w_R, "R", EoS);
    }
    return std::make_pair(S_L, S_R);
}
