#include "MHD1DSolver.H"
#define DEBUG(message) std::cout << "You are at line " << __LINE__ << " doing: " << #message << std::endl

// Simulation setup
int nCells = 200;
double x0 = 0.0;
double x1 = 1.0; // others
// double x1 = 800.0; // test_BrioWu
double C = 0.8;
int ng;
double tStop;
std::array<std::string, nVar> boundary = {"T", "T", "T", "T", "T", "T", "T", "T"}; // T = transmissive, R = reflective, P = periodic

// Plotting parameters
bool animate = false;

// Specify gas parameters
std::string type = "tabulated"; // stiffened, ideal, tabulated
double gamma = 1.4;             // MHD_test1 and toro
// double gamma = 2; // test_BrioWu
// double gamma = 5.0 / 3.0; // test_RyuJones
double p_ref = 3e8;

int main()
{
    // EoS_ptr EoS = nullptr;
    // createEoS(EoS, type, gamma, p_ref);
    // // VectOfDouble w = {0.125, 0.0, 0.0, 0.0, 0.1*101325, 0.0, 0.0, 0.0};
    // VectOfDouble w = {1.25, 0.0, 0.0, 0.0, 101325, 0.0, 0.0, 0.0};
    // // VectOfDouble w = {1.0, 0.0, 0.0, 0.0, 100.0*atm, 0.0, 0.0, 0.0};
    // VectOfDouble u = EoS->primitiveToConservative(w);
    // std::cout << "u: " << u[0] << " " << u[1] << " " << u[2] << " " << u[3] << " " << u[4] << " " << u[5] << " " << u[6] << " " << u[7] << std::endl;
    // VectOfDouble wprime = EoS->conservativeToPrimitive(u);
    // std::cout << "w_prime: " << wprime[0] << " " << wprime[1] << " " << wprime[2] << " " << wprime[3] << " " << wprime[4] << " " << wprime[5] << " " << wprime[6] << " " << wprime[7] << std::endl;
    // std::cout << EoS->c_s(w) << std::endl;

    solveEuler("SLIC", "Toro_1", "1Dcartesian");
    // solveEuler("SLIC", "Toro_2", "1Dcartesian");
    // solveEuler("SLIC", "Toro_3", "1Dcartesian");
    solveEuler("SLIC", "Toro_4", "1Dcartesian");
    // solveEuler("SLIC", "Toro_5", "1Dcartesian");
    // solveEuler("SLIC", "MHD_test2", "1Dcartesian");

    // solveEuler("SLIC", "MHD_test1", "1Dcartesian");
    // solveEuler("SLIC", "test_BrioWu1", "1Dcartesian");
    // solveEuler("SLIC", "test_BrioWu2", "1Dcartesian");
    // solveEuler("SLIC", "test_RyuJones", "1Dcartesian");
};

/* -------------------------------------------------------------------*/
/*                                                                    */
/*                              Solver                                */
/*                                                                    */
/* -------------------------------------------------------------------*/

void solveEuler(std::string scheme, std::string initial, std::string geometry)
{
    // Define the variables
    EoS_ptr EoS = nullptr;
    method_ptr numMethod = nullptr;
    createEoS(EoS, type, gamma, p_ref);
    createNumericalMethod(numMethod, scheme);
    ng = numMethod->ng();
    ODE_ptr source = std::make_shared<ODE_solver>(geometry, ng);

    VectofVectDouble u(nCells + 2 * ng, VectOfDouble(nVar));
    VectofVectDouble uPlus1(nCells + 2 * ng, VectOfDouble(nVar));
    VectofVectDouble uBar(nCells + 2 * ng, VectOfDouble(nVar));

    // Construct the FV domain
    double dx;
    VectOfDouble xCells(nCells + 2 * ng);
    ConstructFVDomain(dx, xCells);

    // Set the initial conditions
    setInitialConditions(xCells, u, tStop, initial, numMethod, EoS);

    std::ostringstream filename = generate_filename(numMethod, initial, geometry, EoS);
    std::ofstream output("./" + filename.str() + ".dat");
    std::cout << "Solving " << filename.str() << std::endl;
    int step = 0;
    double dt;
    double t = 0.0;
    // Time stepping loop
    do
    {
        // Store the data
        writeToFile(output, t, xCells, u, EoS);

        // Compute time step
        computeDt(u, dx, dt, step, EoS);
        t += dt;
        step++;

        // Update
        source->sourceUpdate(u, uBar, xCells, dt, EoS);
        applyBoundaryConditions(uBar, boundary, nCells);
        numMethod->UpdatewithFluxes(uBar, uPlus1, dx, dt, EoS);
        applyBoundaryConditions(uPlus1, boundary, nCells);

        u = uPlus1;
    } while (t < tStop);

    output.close();
    plot(filename, initial);
};

/* -------------------------------------------------------------------*/
/*                                                                    */
/*                      PrintUtility Routine                          */
/*                                                                    */
/* -------------------------------------------------------------------*/
std::ostringstream generate_filename(method_ptr &numMethod, std::string initial, std::string &geometry, EoS_ptr &EoS)
{
    std::ostringstream filename;
    filename << "Euler_" << numMethod->name() << "_" << initial
             << "_n=" << nCells
             << std::fixed << std::setprecision(3) << "_t=" << tStop
             << "_C=" << C
             << "_boundary_" << boundary[0]
             << "_" << geometry
             << "_" << EoS->info();
    return filename;
};

void writeToFile(std::ofstream &output, double t, VectOfDouble xCells, VectofVectDouble &u, EoS_ptr &EoS)
{
    for (int i = ng; i < nCells + ng; i++)
    {
        double x = xCells[i];
        VectOfDouble w = EoS->conservativeToPrimitive(u[i]);
        double B_mag = sqrt(w[BX] * w[BX] + w[BY] * w[BY] + w[BZ] * w[BZ]);
        output << t << " " << x << " " << w[RHO] << " " << w[VX] << " " << w[VY] << " " << w[VZ] << " "
               << w[P] << " " << w[BX] << " " << w[BY] << " " << w[BZ] << " " << B_mag << std::endl;
    }
    output << std::endl
           << std::endl;
}

void plot(std::ostringstream &filename, std::string initial)
{
    if (animate)
    {
        system(("gnuplot -c gif.gp " + filename.str()).c_str());
    };

    if (initial == "isolated_1" || initial == "isolated_2" || initial == "woodward")
    {
        system(("gnuplot -c png_isolated.gp " + filename.str()).c_str());
    }
    else
    {
        system(("gnuplot -c png.gp " + filename.str()).c_str());
    }

    remove(("./" + filename.str() + ".dat").c_str());
};
;

/* -------------------------------------------------------------------*/
/*                                                                    */
/*                  Domain Construction Routine                       */ // OK
/*                                                                    */
/* -------------------------------------------------------------------*/

void ConstructFVDomain(double &dx, VectOfDouble &xCells)
{
    dx = (x1 - x0) / nCells;
    for (int i = 0; i < xCells.size(); i++)
    {
        xCells[i] = x0 + dx * (i - ng + 0.5);
    };
};

/* -------------------------------------------------------------------*/
/*                                                                    */
/*                  Initial Conditions Routine                        */ // OK
/*                                                                    */
/* -------------------------------------------------------------------*/
void setInitialConditions(const VectOfDouble &xCells, VectofVectDouble &u, double &tStop, std::string &initial, method_ptr &numMethod, EoS_ptr &EoS)
{
    VectOfDouble w_L(nVar), w_R(nVar), u_L(nVar), u_R(nVar);
    double xDiscontinuity;
    double dim_factor = sqrt(atm);
    if (initial == "MHD_test1")
    {
        xDiscontinuity = 0.5;
        tStop = 0.25;
        w_L = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0};
        w_R = {0.125, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0};
    }
    else if (initial == "Toro_1")
    {
        xDiscontinuity = 0.5;
        tStop = 0.25;
        tStop = tStop / dim_factor;
        w_L = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0};
        w_R = {0.125, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0};
        w_L[VX] *= dim_factor;
        w_R[VX] *= dim_factor;
        w_L[P] = w_L[P] * atm;
        w_R[P] = w_R[P] * atm;
    }
    else if (initial == "Toro_2")
    {
        // May leave the range of validity as it approaches vacuum
        xDiscontinuity = 0.5;
        tStop = 0.15;
        tStop = tStop / dim_factor;
        w_L = {1.0, -2.0, 0.0, 0.0, 0.4, 0.0, 0.0, 0.0};
        w_R = {1.0, 2.0, 0.0, 0.0, 0.4, 0.0, 0.0, 0.0};
        w_L[VX] = w_L[VX] * dim_factor;
        w_L[P] = w_L[P] * atm;
        w_R[P] = w_R[P] * atm;
        w_R[VX] = w_R[VX] * dim_factor;
    }
    else if (initial == "Toro_3")
    {
        // Fail as not enough range of pressure (lower) in tabulated EoS
        xDiscontinuity = 0.5;
        tStop = 0.012;
        tStop = tStop / dim_factor;
        w_L = {1.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0};
        w_R = {1.0, 0.0, 0.0, 0.0, 0.001, 0.0, 0.0, 0.0};
        w_L[VX] *= dim_factor;
        w_L[P] = w_L[P] * atm;
        w_R[P] = w_R[P] * atm;
        w_R[VX] *= dim_factor;
    }
    else if (initial == "Toro_4")
    {
        xDiscontinuity = 0.5;
        tStop = 0.035;
        tStop = tStop / dim_factor;
        w_L = {1.0, 0.0, 0.0, 0.0, 0.01, 0.0, 0.0, 0.0};
        w_R = {1.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0};
        w_L[VX] *= dim_factor;
        w_L[P] = w_L[P] * atm;
        w_R[P] = w_R[P] * atm;
        w_R[VX] *= dim_factor;
    }
    else if (initial == "Toro_5")
    {
        // Fail as not enough range of pressure (upper) in tabulated EoS
        xDiscontinuity = 0.5;
        tStop = 0.035;
        tStop = tStop / dim_factor;
        w_L = {5.99924, 19.5875, 0.0, 0.0, 460.894, 0.0, 0.0, 0.0};
        w_R = {5.99924, -6.19633, 0.0, 0.0, 46.0950, 0.0, 0.0, 0.0};
        w_L[VX] *= dim_factor;
        w_L[P] = w_L[P] * atm;
        w_R[P] = w_R[P] * atm;
        w_R[VX] *= dim_factor;
    }
    else if (initial == "test_BrioWu1")
    {
        xDiscontinuity = 400;
        tStop = 80;
        w_L = {1.0, 0.0, 0.0, 0.0, 1.0, 0.75, 1.0, 0.0};
        w_R = {0.125, 0.0, 0.0, 0.0, 0.1, 0.75, -1.0, 0.0};
    }
    else if (initial == "test_BrioWu2")
    {
        xDiscontinuity = 400;
        tStop = 80;
        w_L = {1.0, 0.0, 0.0, 0.0, 1.0, 0.75, 0.0, 1.0};
        w_R = {0.125, 0.0, 0.0, 0.0, 0.1, 0.75, 0.0, -1.0};
    }
    else if (initial == "test_RyuJones")
    {
        xDiscontinuity = 0.5;
        tStop = 0.2;
        w_L = {1.08, 1.2, 0.01, 0.5, 3.6 / sqrt(4 * M_PI), 2 / sqrt(4 * M_PI), 2 / sqrt(4 * M_PI), 0.95};
        w_R = {1.0, 0.0, 0.0, 0.0, 4 / sqrt(4 * M_PI), 2 / sqrt(4 * M_PI), 2 / sqrt(4 * M_PI), 1.0};
    }
    u_L = EoS->primitiveToConservative(w_L);
    u_R = EoS->primitiveToConservative(w_R);

    for (int i = 0; i < u.size(); i++)
    {
        double x = xCells[i];
        if (x <= xDiscontinuity)
        {
            u[i] = u_L;
        }
        else
        {
            u[i] = u_R;
        }
    }
    return;
}

/* -------------------------------------------------------------------*/
/*                                                                    */
/*                  Compute Time Step Routine                         */ // OK
/*                                                                    */
/* -------------------------------------------------------------------*/

void computeDt(VectofVectDouble &u, double &dx, double &dt, int &step, EoS_ptr &EoS)
{
    VectOfDouble a(nCells);
    double a_max;

    for (int i = 0; i < nCells; i++)
    {
        int pos = i + ng;
        VectOfDouble w = EoS->conservativeToPrimitive(u[pos]);
        double v = sqrt(w[VX] * w[VX] + w[VY] * w[VY] + w[VZ] * w[VZ]);
        // Can write a function to return all wavespeeds in the system and take the largest,
        // in this case, we know c_f is the fastest (but it is not general)
        a[i] = v + EoS->c_f(w);
    };
    a_max = *std::max_element(a.begin(), a.end());

    if (step < 5)
    {
        dt = 0.4 * C * dx / abs(a_max);
    }
    else
    {
        dt = C * dx / abs(a_max);
    }
}

/* -------------------------------------------------------------------*/
/*                                                                    */
/*                 Boundary Conditions Routine                        */ // OK
/*                                                                    */
/* -------------------------------------------------------------------*/

void applyBoundaryConditions(VectofVectDouble &u, std::array<std::string, nVar> boundary, int nCells)
{
    for (int i = 0; i < nVar; i++)
    {
        for (int ig = 0; ig < ng; ig++)
        {
            if (boundary[i] == "T")
            {
                u[ig][i] = u[ng + ig][i];
                u[nCells + ng + ig][i] = u[nCells + ig][i];
            }
            else if (boundary[i] == "R")
            {
                u[ig][i] = -u[ng + ig][i];
                u[nCells + ng + ig][i] = -u[nCells + ig][i];
            }
            else if (boundary[i] == "P")
            {
                u[ig][i] = u[nCells + ig][i];
                u[nCells + ng + ig][i] = u[ng + ig][i];
            };
        }
    }
};
