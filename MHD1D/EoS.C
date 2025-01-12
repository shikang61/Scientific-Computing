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
    else if (type == "tabulated")
    {
        EoS = std::make_shared<TabulatedEoS>();
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
double EquationOfState::c_a(const VectOfDouble &w) const
{
    return abs(w[BX]) / sqrt(w[RHO]);
};

double EquationOfState::c_f(const VectOfDouble &w) const
{
    double B = w[BX] * w[BX] + w[BY] * w[BY] + w[BZ] * w[BZ];
    double c_s2 = c_s(w) * c_s(w);
    return sqrt(0.5 * (c_s2 + B * B / w[RHO] + sqrt(pow(c_s2 + B * B / w[RHO], 2) - 4 * c_s2 * w[BX] * w[BX] / w[RHO])));
};

double EquationOfState::c_sl(const VectOfDouble &w) const
{
    double B = w[BX] * w[BX] + w[BY] * w[BY] + w[BZ] * w[BZ];
    double c_s2 = c_s(w) * c_s(w);
    return sqrt(0.5 * (c_s2 + B * B / w[RHO] - sqrt(pow(c_s2 + B * B / w[RHO], 2) - 4 * c_s2 * w[BX] * w[BX] / w[RHO])));
};

double EquationOfState::enthalpy(const VectOfDouble &u, const VectOfDouble &w) const
{
    return (u[ENE] + w[P]) / w[RHO];
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
    w[VX] = u[MOMX] / u[RHO];
    w[VY] = u[MOMY] / u[RHO];
    w[VZ] = u[MOMZ] / u[RHO];
    w[P] = (m_gamma - 1) * (u[ENE] - 0.5 * (u[MOMX] * u[MOMX] + u[MOMY] * u[MOMY] + u[MOMZ] * u[MOMZ]) / u[RHO] - 0.5 * (u[BY] * u[BY] + u[BZ] * u[BZ] + u[BX] * u[BX]));
    w[BX] = u[BX];
    w[BY] = u[BY];
    w[BZ] = u[BZ];
    return w;
};

VectOfDouble IdealGasEoS::primitiveToConservative(const VectOfDouble &w) const
{
    VectOfDouble u(nVar);
    u[RHO] = w[RHO];
    u[MOMX] = w[RHO] * w[VX];
    u[MOMY] = w[RHO] * w[VY];
    u[MOMZ] = w[RHO] * w[VZ];
    u[ENE] = w[P] / (m_gamma - 1) + 0.5 * w[RHO] * (w[VX] * w[VX] + w[VY] * w[VY] + w[VZ] * w[VZ]) + 0.5 * (w[BY] * w[BY] + w[BZ] * w[BZ] + w[BX] * w[BX]);
    u[BX] = w[BX];
    u[BY] = w[BY];
    u[BZ] = w[BZ];
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
    w[VX] = u[MOMX] / u[RHO];
    w[VY] = u[MOMY] / u[RHO];
    w[VZ] = u[MOMZ] / u[RHO];
    w[P] = (m_gamma - 1) * (u[ENE] - 0.5 * (u[MOMX] * u[MOMX] + u[MOMY] * u[MOMY] + u[MOMZ] * u[MOMZ]) / u[RHO] - 0.5 * (u[BY] * u[BY] + u[BZ] * u[BZ] + u[BX] * u[BX])) - m_gamma * m_p_ref;
    w[BX] = u[BX];
    w[BY] = u[BY];
    w[BZ] = u[BZ];
    return w;
};

VectOfDouble StiffenedGasEoS::primitiveToConservative(const VectOfDouble &w) const
{
    VectOfDouble u(nVar);
    u[RHO] = w[RHO];
    u[MOMX] = w[RHO] * w[VX];
    u[MOMY] = w[RHO] * w[VY];
    u[MOMZ] = w[RHO] * w[VZ];
    u[ENE] = (w[P] + m_gamma * m_p_ref) / (m_gamma - 1) + 0.5 * w[RHO] * (w[VX] * w[VX] + w[VY] * w[VY] + w[VZ] * w[VZ]) + 0.5 * (w[BY] * w[BY] + w[BZ] * w[BZ] + w[BX] * w[BX]);
    u[BX] = w[BX];
    u[BY] = w[BY];
    u[BZ] = w[BZ];
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
/*                          Tabulated EoS                             */
/*                                                                    */
/* -------------------------------------------------------------------*/
TabulatedEoS::TabulatedEoS()
{
    readData("Plasma19 EoS.txt");
    std::cout << "Using tabulated EoS. Data read: Plasma19 EoS.txt" << std::endl;
}

VectOfDouble TabulatedEoS::conservativeToPrimitive(const VectOfDouble &u) const
{
    double frac_rho_index, frac_p_index;
    VectOfDouble w(nVar);
    w[RHO] = u[RHO];
    w[VX] = u[MOMX] / u[RHO];
    w[VY] = u[MOMY] / u[RHO];
    w[VZ] = u[MOMZ] / u[RHO];
    frac_rho_index = search_density_index(u[RHO]);
    double epsilon = (u[ENE] - 0.5 * (u[MOMX] * u[MOMX]) / u[RHO]) / u[RHO];
    std::tie(frac_p_index, w[P]) = guess_pressure_from_rho_and_E(frac_rho_index, epsilon, specificEnergies);
    w[BX] = u[BX];
    w[BY] = u[BY];
    w[BZ] = u[BZ];
    return w;
};

VectOfDouble TabulatedEoS::primitiveToConservative(const VectOfDouble &w) const
{
    VectOfDouble u(nVar);
    double frac_p_index, frac_rho_index;
    u[RHO] = w[RHO];
    u[MOMX] = w[RHO] * w[VX];
    u[MOMY] = w[RHO] * w[VY];
    u[MOMZ] = w[RHO] * w[VZ];
    frac_p_index = search_pressure_index(w[P]);
    frac_rho_index = search_density_index(w[RHO]);
    double epsilon = bilinear_interp(frac_p_index, frac_rho_index, specificEnergies);
    u[ENE] = w[RHO] * epsilon + 0.5 * w[RHO] * w[VX] * w[VX];
    u[BX] = w[BX];
    u[BY] = w[BY];
    u[BZ] = w[BZ];
    return u;
};

double TabulatedEoS::c_s(const VectOfDouble &w) const
{
    VectOfDouble u(nVar);
    double pressure, frac_p_index, frac_rho_index;
    u = primitiveToConservative(w);
    frac_rho_index = search_density_index(w[RHO]);
    double epsilon = (u[ENE] - 0.5 * (u[MOMX] * u[MOMX]) / u[RHO]) / u[RHO];
    std::tie(frac_p_index, pressure) = guess_pressure_from_rho_and_E(frac_rho_index, epsilon, specificEnergies);
    return bilinear_interp(frac_p_index, frac_rho_index, soundSpeeds);
};

std::string TabulatedEoS::info() const
{
    std::ostringstream oss;
    oss << m_name;
    return oss.str();
};

void TabulatedEoS::readData(const std::string &filename)
{
    int ND;
    int NP;
    std::ifstream inputFile(filename);
    inputFile >> ND >> NP;
    int n = 8;
    std::string line;

    for (int i = 0; i < n && std::getline(inputFile, line); ++i)
    {
        // Do nothing, just skip the lines
    }

    // Resize vectors based on ND and NP
    densities.resize(ND);
    pressures.resize(NP);
    soundSpeeds.resize(NP, VectOfDouble(ND));      // Store sound speeds
    specificEnergies.resize(NP, VectOfDouble(ND)); // Store specific internal energies

    // Read densities
    for (int i = 0; i < ND; ++i)
    {
        inputFile >> densities[i];
    }
    // Read pressures
    for (int i = 0; i < NP; ++i)
    {
        inputFile >> pressures[i];
    }

    // Read tables: sound speeds and specific internal energies
    for (int row = 0; row < ND; ++row)
    {
        for (int col = 0; col < NP; ++col)
        {
            inputFile >> soundSpeeds[row][col];
        }
    }
    for (int row = 0; row < ND; ++row)
    {
        for (int col = 0; col < NP; ++col)
        {
            inputFile >> specificEnergies[row][col];
        }
    }
    inputFile.close();
}
double TabulatedEoS::search_pressure_index(const double &p) const
{
    int Rptr;
    for (int i = 1; i < pressures.size(); ++i)
    {
        Rptr = i;
        if (p <= pressures[Rptr])
        {
            return (pressures[Rptr] - p) / (pressures[Rptr] - pressures[Rptr - 1]) * (Rptr - 1) + (p - pressures[Rptr - 1]) / (pressures[Rptr] - pressures[Rptr - 1]) * Rptr;
        };
    };
    return 0.0;
}
double TabulatedEoS::search_density_index(const double &rho) const
{
    int Rptr;
    for (int i = 1; i < densities.size(); ++i)
    {
        Rptr = i;
        if (rho <= densities[Rptr])
        {
            return (densities[Rptr] - rho) / (densities[Rptr] - densities[Rptr - 1]) * (Rptr - 1) + (rho - densities[Rptr - 1]) / (densities[Rptr] - densities[Rptr - 1]) * Rptr;
        };
    };
    return 0.0;
}

double TabulatedEoS::bilinear_interp(const double &frac_p_index, const double &frac_rho_index, const VectofVectDouble &table) const
{
    int p_L = std::floor(frac_p_index);
    int p_R = std::ceil(frac_p_index);
    int rho_L = std::floor(frac_rho_index);
    int rho_R = std::ceil(frac_rho_index);
    double Q11 = table[p_L][rho_L];
    double Q21 = table[p_L][rho_R];
    double Q12 = table[p_R][rho_L];
    double Q22 = table[p_R][rho_R];
    // return Q11*(p_R-frac_p_index)*(rho_R-frac_rho_index) + Q21*(frac_p_index-p_L)*(rho_R-frac_rho_index) + Q12*(p_R-frac_p_index)*(frac_rho_index-rho_L) + Q22*(frac_p_index-p_L)*(frac_rho_index-rho_L);
    double N = (Q21 - Q11) * (frac_rho_index - rho_L) + Q11;
    double S = (Q22 - Q12) * (frac_rho_index - rho_L) + Q12;
    return (S - N) * (frac_p_index - p_L) + N;
}

std::tuple<double, double> TabulatedEoS::guess_pressure_from_rho_and_E(const double &frac_rho_index, const double &target, const VectofVectDouble &table) const
{
    int rho_L = std::floor(frac_rho_index);
    // std::cout << rho_L << std::endl;
    int p_R = 1;
    int p_L = 0;
    for (int i = 1; i < table[0].size(); i++)
    {
        double temp_i = i + 0.0;
        double temp_im1 = i - 1.0;
        // std::cout << frac_rho_index << std::endl;
        double bigger = bilinear_interp(temp_i, frac_rho_index, table);
        double smaller = bilinear_interp(temp_im1, frac_rho_index, table);
        // std::cout << "i: " << i << " " << smaller << " " << bigger << " " << target << std::endl;
        if ((bigger >= target) && (smaller <= target))
        {
            p_R = i;
            break;
        }
        // else
        // {
        //     // std::cout << "No pressure found." << std::endl;
        //     continue;
        // }
    };
    // std::cout << p_R << std::endl;
    p_L = p_R - 1;
    double temp_p_L = p_L + 0.0;
    double temp_p_R = p_R + 0.0;
    double tolerance = 1e-6;
    double frac_p_index = 0.5 * (temp_p_L + temp_p_R);
    int iteration = 0;
    int max_iterations = 1e7;
    while (abs(target - bilinear_interp(frac_p_index, frac_rho_index, table)) / target > tolerance)
    {
        // std::cout << "target: " << target << " interp: " << bilinear_interp(frac_p_index, frac_rho_index, table) << std::endl;
        // std::cout << temp_p_R << std::endl;
        // std::cout << frac_p_index << std::endl;
        if (table[p_L][rho_L] > table[p_R][rho_L])
        {
            if (target > bilinear_interp(frac_p_index, frac_rho_index, table))
            {
                // std::cout << "left" << std::endl;
                temp_p_R = frac_p_index;
                // std::cout << frac_p_index << std::endl;
                frac_p_index = 0.5 * (temp_p_L + temp_p_R);
            }
            else
            {
                // std::cout << "right" << std::endl;
                temp_p_L = frac_p_index;
                // std::cout << frac_p_index << std::endl;
                frac_p_index = 0.5 * (temp_p_R + temp_p_L);
            };
        }
        else
        {
            if (target > bilinear_interp(frac_p_index, frac_rho_index, table))
            {
                temp_p_L = frac_p_index;
                frac_p_index = 0.5 * (temp_p_R + temp_p_L);
            }
            else
            {
                temp_p_R = frac_p_index;
                frac_p_index = 0.5 * (temp_p_R + temp_p_L);
            };
        };
        if (iteration > max_iterations)
        {
            std::cout << "Max iterations reached without convergence." << std::endl;
            break;
        };
        iteration++;
    };
    double p = pressures[p_L] * (p_R - frac_p_index) + pressures[p_R] * (frac_p_index - p_L);
    return std::tie(frac_p_index, p);
}