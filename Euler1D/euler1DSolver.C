#include "euler1DSolver.H"
#define DEBUG(message) std::cout << "You are at line " << __LINE__ << " doing: " << #message << std::endl
#define SHOW(variable) std::cout << "You are at line " << __LINE__ << " doing: " << variable << std::endl

// Simulation setup
int nCells = 310;
double x0 = 0.0;
double x1 = 1.0;
double C = 0.8;
int ng;
double tStop;
std::array<std::string, nVar> boundary = {"T", "T", "T"}; // T = transmissive, R = reflective, P = periodic

// Plotting parameters
bool animate = true;
bool exactSolution = false;

// Specify gas parameters - can use a function to set these?
Mat_ptr material_1 = std::make_shared<Material>("ideal", 1.4, 0.0, false);

int main()
{
    solveEuler("MUSCL_Godunov", "toro1", "1Dcartesian", material_1);
    // solveEuler("Godunov", "toro2", "1Dcartesian", material_1);
    // solveEuler("Godunov", "toro3", "1Dcartesian", material_1);
    // solveEuler("SLIC", "toro4", "1Dcartesian", material_1);
    // solveEuler("Godunov", "toro5", "1Dcartesian", material_1);
    // solveEuler("FORCE", "isolated_1", "1Dcartesian");
    // solveEuler("FVTV", "isolated_1", "1Dcartesian");
    // solveEuler("SLIC", "toro2", "1Dcartesian");
    // solveEuler("SLIC", "toro3", "1Dcartesian");
    // solveEuler("SLIC", "toro4", "1Dcartesian");
    // solveEuler("SLIC", "toro5", "1Dcartesian");

    // solveEuler("No", "testA", "1Dcartesian");
};

/* -------------------------------------------------------------------*/
/*                                                                    */
/*                              Solver                                */
/*                                                                    */
/* -------------------------------------------------------------------*/

void solveEuler(std::string scheme, std::string initial, std::string geometry, Mat_ptr &material)
{
    EoS_ptr EoS = nullptr;
    Method_ptr numMethod = nullptr;
    createEoS(EoS, material->type(), material->gamma(), material->p_ref());
    createNumericalMethod(numMethod, scheme);
    ng = numMethod->ng();
    ODE_ptr source = std::make_shared<ODE_solver>(geometry, ng);

    // Define the variables
    VectofVectDouble u(nCells + 2 * ng, VectOfDouble(nVar));
    VectofVectDouble uPlus1(nCells + 2 * ng, VectOfDouble(nVar));
    VectofVectDouble uBar(nCells + 2 * ng, VectOfDouble(nVar));
    VectofVectDouble w_ExactRiemann(nCells + 2 * ng, VectOfDouble(nVar));

    // Construct the FV domain
    double dx;
    VectOfDouble xCells(nCells + 2 * ng);
    ConstructFVDomain(dx, xCells);

    // Set the initial conditions
    setInitialConditions(xCells, u, tStop, w_ExactRiemann, initial, numMethod, exactSolution, EoS);

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
        writeToFile(output, t, xCells, u, w_ExactRiemann, EoS);

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
/*                  Initial Conditions Routine                        */ // OK
/*                                                                    */
/* -------------------------------------------------------------------*/

void setInitialConditions(const VectOfDouble &xCells, VectofVectDouble &u, double &tStop, VectofVectDouble &w_ExactRiemann, std::string &initial, Method_ptr &numMethod, bool exactSolution, EoS_ptr &EoS)
{
    VectOfDouble w_L(nVar), w_R(nVar), u_L(nVar), u_R(nVar);
    double xDiscontinuity;
    if (initial == "woodward")
    {
        VectOfDouble w_C(nVar), u_C(nVar);
        double lDiscontinuity = 0.1;
        double rDiscontinuity = 0.9;
        xDiscontinuity = 0.5;
        tStop = 0.038;
        w_L = {1.0, 0.0, 1000.0};
        w_C = {1.0, 0.0, 0.01};
        w_R = {1.0, 0.0, 100.0};
        u_L = EoS->primitiveToConservative(w_L);
        u_C = EoS->primitiveToConservative(w_C);
        u_R = EoS->primitiveToConservative(w_R);

        for (int i = 0; i < u.size(); i++)
        {
            double x = xCells[i];
            if (x <= lDiscontinuity)
            {
                u[i] = u_L;
            }
            else if (x > lDiscontinuity && x < rDiscontinuity)
            {
                u[i] = u_C;
            }
            else
            {
                u[i] = u_R;
            };
        };
        return;
    }

    if (initial == "toro1")
    {
        if (EoS->name() == "stiffened")
        {
            xDiscontinuity = 0.5;
            tStop = 2e-4;
            w_L = {1500.0, 0.0, 3000.0 * 101325.0};
            w_R = {1000.0, 0.0, 101325.0};
        }
        else
        {
            xDiscontinuity = 0.5;
            tStop = 0.25;
            w_L = {1.0, 0.0, 1.0};
            w_R = {0.125, 0.0, 0.1};
        }
    }
    else if (initial == "toro2")
    {
        if (EoS->name() == "stiffened")
        {
            xDiscontinuity = 0.5;
            tStop = 2e-4;
            w_L = {1000.0, -350.0, 2 * 101325.0};
            w_R = {1000.0, 350.0, 2 * 101325.0};
        }
        else
        {
            xDiscontinuity = 0.5;
            tStop = 0.15;
            w_L = {1.0, -2.0, 0.4};
            w_R = {1.0, 2.0, 0.4};
        }
    }
    else if (initial == "toro3")
    {
        xDiscontinuity = 0.5;
        tStop = 0.012;
        w_L = {1.0, 0.0, 1000.0};
        w_R = {1.0, 0.0, 0.01};
    }
    else if (initial == "toro4")
    {
        xDiscontinuity = 0.5;
        tStop = 0.035;
        w_L = {1.0, 0.0, 0.01};
        w_R = {1.0, 0.0, 100.0};
    }
    else if (initial == "toro5")
    {
        xDiscontinuity = 0.5;
        tStop = 0.035;
        w_L = {5.99924, 19.5975, 460.894};
        w_R = {5.99242, -6.19633, 46.0950};
    }
    else if (initial == "test1")
    {
        xDiscontinuity = 0.5;
        tStop = 0.5;
        w_L = {1.0, 0.0, 1.0};
        w_R = {1.0, 0.0, 1.0};
    }
    else if (initial == "test_2")
    {
        xDiscontinuity = 0.5;
        tStop = 0.5;
        w_L = {2.0, 1.0, 3.0};
        w_R = {2.0, 1.0, 3.0};
    }
    else if (initial == "test3")
    {
        xDiscontinuity = 0.5;
        tStop = 0.5;
        w_L = {1.0, 0.0, 1.0};
        w_R = {0.1, 0.0, 1.0};
    }
    else if (initial == "test4")
    {
        xDiscontinuity = 0.25;
        tStop = 0.5;
        w_L = {1.0, 1.0, 1.0};
        w_R = {0.1, 1.0, 1.0};
    }
    else if (initial == "testA")
    {
        xDiscontinuity = 0.505;
        tStop = 0.0007;
        w_L = {1.0, 0.0, 1e5};
        w_R = {0.125, 0.0, 1e4};
    }
    else if (initial == "testB")
    {
        xDiscontinuity = 0.505;
        tStop = 0.0007;
        w_L = {1.0, 0.0, 1e5};
        w_R = {0.125, 0.0, 1e4};
    }
    else if (initial == "isolated_1")
    {
        xDiscontinuity = 0.5;
        tStop = 2.0;
        w_L = {1.4, 0.0, 1.0};
        w_R = {1.0, 0.0, 1.0};
    }
    else if (initial == "isolated_2")
    {
        xDiscontinuity = 0.4;
        tStop = 2.0;
        w_L = {1.4, 0.1, 1.0};
        w_R = {1.0, 0.1, 1.0};
    }
    else if (initial == "explosion")
    {
        xDiscontinuity = 0.4;
        tStop = 0.25;
        w_L = {1.0, 0.0, 1.0};
        w_R = {0.125, 0.0, 0.1};
    }
    else
    {
        std::cout << "Initial condition not set properly" << std::endl;
        exit(1);
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
        if (exactSolution)
        {
            w_ExactRiemann[i] = exactRiemannSolver(x, xDiscontinuity, tStop, w_L, w_R, EoS);
        }
    }
    return;
}

/* -------------------------------------------------------------------*/
/*                                                                    */
/*                      PrintUtility Routine                          */
/*                                                                    */
/* -------------------------------------------------------------------*/
std::ostringstream generate_filename(Method_ptr &numMethod, std::string initial, std::string &geometry, EoS_ptr &EoS)
{
    std::ostringstream filename;
    filename << "Euler_" << numMethod->name() << "_" << initial
             << "_n=" << nCells
             << std::fixed << std::setprecision(5) << "_t=" << tStop
             << "_C=" << C
             << "_boundary_" << boundary[0] << boundary[1] << boundary[2]
             << "_" << geometry
             << "_EoS:" << EoS->info();
    return filename;
};

void writeToFile(std::ofstream &output, double t, VectOfDouble xCells, VectofVectDouble &u, VectofVectDouble &w_ExactRiemann, EoS_ptr &EoS)
{
    for (int i = ng; i < nCells + ng; i++)
    {
        double x = xCells[i];
        VectOfDouble w = EoS->conservativeToPrimitive(u[i]);
        output << t << " " << x << " " << w[RHO] << " " << w[V] << " " << w[P] << " "
               << w_ExactRiemann[i][RHO] << " " << w_ExactRiemann[i][V] << " " << w_ExactRiemann[i][P] << std::endl;
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
    }
};

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
        a[i] = abs(w[V]) + EoS->c_s(w);
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
