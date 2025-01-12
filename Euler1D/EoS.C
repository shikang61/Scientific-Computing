#include "EoS.H"

void createEoS(EoS_ptr &EoS, const std::string &type, const double &gamma, const double &p_ref, const double &e_ref)
{
    if (type == "ideal")
    {
        EoS = std::make_shared<IdealGasEoS>(gamma);
    }
    else if (type == "stiffened")
    {
        EoS = std::make_shared<StiffenedGasEoS>(gamma, p_ref);
    }
    else if (type == "JWL")
    {
        EoS = std::make_shared<JWLEoS>(gamma);
    }
    else
    {
        std::cout << "Gas not initialised" << std::endl;
        exit(1);
    }
}

/* -------------------------------------------------------------------*/
/*                                                                    */
/*                          Base class                                */
/*                                                                    */
/* -------------------------------------------------------------------*/
EquationOfState::EquationOfState(double gamma, double p_ref) : m_gamma(gamma), m_p_ref(p_ref) {}
EquationOfState::EquationOfState(double gamma) : m_gamma(gamma), m_p_ref(0.0) {}
EquationOfState::EquationOfState() : m_gamma(1.4), m_p_ref(0.0) {}

std::string EquationOfState::name() const
{
    return m_name;
}

double EquationOfState::p_ref() const
{
    return m_p_ref;
}

double EquationOfState::gamma() const
{
    return m_gamma;
}

double EquationOfState::enthalpy(const VectOfDouble &u, const VectOfDouble &w) const
{
    return (u[ENE] + w[P] + m_p_ref) / w[RHO];
};

/* -------------------------------------------------------------------*/
/*                                                                    */
/*                          Ideal Gas EoS                             */
/*                                                                    */
/* -------------------------------------------------------------------*/

IdealGasEoS::IdealGasEoS(double gamma) : EquationOfState(gamma) {}
VectOfDouble IdealGasEoS::conservativeToPrimitive(const VectOfDouble &u) const
{
    VectOfDouble w(nVar);
    w[RHO] = u[RHO];
    w[V] = u[MOM] / u[RHO];
    w[P] = (m_gamma - 1) * (u[ENE] - 0.5 * u[MOM] * u[MOM] / u[RHO]);
    return w;
};

VectOfDouble IdealGasEoS::primitiveToConservative(const VectOfDouble &w) const
{
    VectOfDouble u(nVar);
    u[RHO] = w[RHO];
    u[MOM] = w[RHO] * w[V];
    u[ENE] = w[P] / (m_gamma - 1) + 0.5 * w[RHO] * w[V] * w[V];
    return u;
};

double IdealGasEoS::c_s(const VectOfDouble &w) const
{
    return sqrt(m_gamma * w[P] / w[RHO]);
};

std::string IdealGasEoS::info() const
{
    std::ostringstream oss;
    oss << m_name << "_" << m_gamma;
    return oss.str();
};

/* -------------------------------------------------------------------*/
/*                                                                    */
/*                        Stiffened Gas EoS                           */
/*                                                                    */
/* -------------------------------------------------------------------*/

StiffenedGasEoS::StiffenedGasEoS(double gamma, double p_ref) : EquationOfState(gamma, p_ref) {};

VectOfDouble StiffenedGasEoS::conservativeToPrimitive(const VectOfDouble &u) const
{
    VectOfDouble w(nVar);
    w[RHO] = u[RHO];
    w[V] = u[MOM] / u[RHO];
    w[P] = (m_gamma - 1) * (u[ENE] - 0.5 * u[MOM] * u[MOM] / u[RHO]) - m_gamma * m_p_ref;
    return w;
};

VectOfDouble StiffenedGasEoS::primitiveToConservative(const VectOfDouble &w) const
{
    VectOfDouble u(nVar);
    u[RHO] = w[RHO];
    u[MOM] = w[RHO] * w[V];
    u[ENE] = (w[P] + m_gamma * m_p_ref) / (m_gamma - 1) + 0.5 * w[RHO] * w[V] * w[V];
    return u;
};

double StiffenedGasEoS::c_s(const VectOfDouble &w) const
{
    return sqrt(m_gamma * (w[P] + m_p_ref) / w[RHO]);
};

std::string StiffenedGasEoS::info() const
{
    std::ostringstream oss;
    oss << m_name << "_" << m_gamma << "_" << m_p_ref;
    return oss.str();
};

double StiffenedGasEoS::enthalpy(const VectOfDouble &u, const VectOfDouble &w) const
{
    return (u[ENE] + w[P] + m_p_ref) / w[RHO];
};

/* -------------------------------------------------------------------*/
/*                                                                    */
/*                               JWL EoS                              */
/*                                                                    */
/* -------------------------------------------------------------------*/

JWLEoS::JWLEoS(double gamma) : EquationOfState(gamma) {};

VectOfDouble JWLEoS::conservativeToPrimitive(const VectOfDouble &u) const
{
    VectOfDouble w(nVar);
    w[RHO] = u[RHO];
    w[V] = u[MOM] / u[RHO];
    w[P] = JWL_pref(w[RHO]) + m_gamma * (u[ENE] - 0.5 * u[MOM] * u[MOM] / u[RHO] - w[RHO] * JWL_eref(w[RHO]));
    return w;
};

VectOfDouble JWLEoS::primitiveToConservative(const VectOfDouble &w) const
{
    VectOfDouble u(nVar);
    u[RHO] = w[RHO];
    u[MOM] = w[RHO] * w[V];
    u[ENE] = (w[P] - JWL_pref(w[RHO])) / m_gamma + 0.5 * w[RHO] * w[V] * w[V] + w[RHO] * JWL_eref(w[RHO]);
    return u;
};

double JWLEoS::c_s(const VectOfDouble &w) const
{
    return sqrt((m_gamma + 1) * (w[P] - JWL_pref(w[RHO])) / w[0] + m_gamma * JWL_pref(w[RHO]) / w[RHO] + dJWL_eref(w[RHO]) - m_gamma * dJWL_eref(w[RHO]));
};

std::string JWLEoS::info() const
{
    std::ostringstream oss;
    oss << m_name << "_" << m_gamma;
    return oss.str();
};

double JWLEoS::JWL_pref(const double &rho) const
{
    double nu = 1 / rho;
    return JWL_A * exp(-JWL_R1 * nu / JWL_nu0) + JWL_B * exp(-JWL_R2 * nu / JWL_nu0);
};

double JWLEoS::JWL_eref(const double &rho) const
{
    double nu = 1 / rho;
    return JWL_nu0 * (JWL_A / JWL_R1 * exp(-JWL_R1 * nu / JWL_nu0) + JWL_B / JWL_R2 * exp(-JWL_R2 * nu / JWL_nu0));
}

double JWLEoS::dJWL_pref(const double &rho) const
{
    double nu = 1 / rho;
    return JWL_A * JWL_R1 * nu * nu / JWL_nu0 * exp(-JWL_R1 * nu / JWL_nu0) + JWL_B * nu * nu * JWL_R2 / JWL_nu0 * exp(-JWL_R2 * nu / JWL_nu0);
}

double JWLEoS::dJWL_eref(const double &rho) const
{
    double nu = 1 / rho;
    return JWL_A * nu * nu / JWL_nu0 * exp(-JWL_R1 * nu / JWL_nu0) + JWL_B * nu * nu / JWL_nu0 * exp(-JWL_R2 * nu / JWL_nu0);
}
