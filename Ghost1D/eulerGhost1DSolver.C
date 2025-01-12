#include "eulerGhost1DSolver.H"
#define DEBUG(message) std::cout << "You are at line " << __LINE__ << " doing: " << #message << std::endl

// Simulation setup
int nCells = 200;
double x0 = 0.0;
double x1 = 1.0;
double C = 0.8;
int ng;
double tStop;
std::array<std::string, nVar> boundary = {"T", "T", "T"}; // T = transmissive, R = reflective, P = periodic

// Plotting parameters
bool animate = true;
bool exactSolution = false;

double gamma;
// Specify gas parameters

// CDS test
// std::string type = "stiffened";
// double gamma1 = 4.4;
// double p_ref1 = 6e8;
// double gamma2 = 1.4;
// double p_ref2 = 0.0;

// heliumSlab test
std::string type = "ideal";
double gamma1 = 1.4;
double p_ref1 = 0.0;
double gamma2 = 1.68;
double p_ref2 = 0.0;

// Water variables
// double gamma = 2.35;
// double p_ref = 1e9;

// JWL variables
// double gamma = 0.25;
// double p_ref = 0.0;

int main()
{

    // solveEuler("SLIC", "stationary1", "1Dcartesian");
    // solveEuler("SLIC", "stationary2", "1Dcartesian");
    // solveEuler("SLIC", "moving1", "1Dcartesian");
    // solveEuler("SLIC", "moving2", "1Dcartesian");
    // solveEuler("SLIC", "moving3", "1Dcartesian");
    // solveEuler("no", "CDS", "1Dcartesian");
    // solveEuler("SLIC", "toro_2", "1Dcartesian");
    // solveEuler("SLIC", "toro_3", "1Dcartesian");
    // solveEuler("SLIC", "toro_4", "1Dcartesian");
    // solveEuler("SLIC", "toro_5", "1Dcartesian");
    // solveEuler("SLIC", "FedkiwA", "1Dcartesian");
    // solveEuler("Godunov", "FedkiwB", "1Dcartesian");
    solveEuler("SLIC", "heliumSlab", "1Dcartesian");

    // solveEuler("SLIC", "toro_2", "1Dcartesian");
    // solveEuler("SLIC", "toro_3", "1Dcartesian");
    // solveEuler("SLIC", "toro_4", "1Dcartesian");
    // solveEuler("SLIC", "toro_5", "1Dcartesian");
    // solveEuler("No", "test_4", "1Dcartesian");
    // solveEuler("SLIC", "test_4", "1Dcartesian");
    // solveEuler("No", "test_4", "1Dcartesian");
    // solveEuler("No", "TestA", "1Dcartesian");
    // solveEuler("No", "TestB", "1Dcartesian");
    // solveEuler("No", "TestC", "1Dcartesian");
    // solveEuler("No", "TestD1", "1Dcartesian");
    // solveEuler("No", "TestD2", "1Dcartesian");
    // solveEuler("No", "CDS", "1Dcartesian");
};

/* -------------------------------------------------------------------*/
/*                                                                    */
/*                              Solver                                */
/*                                                                    */
/* -------------------------------------------------------------------*/

void solveEuler(std::string scheme, std::string initial, std::string geometry)
{
    EoS_ptr EoS_1 = nullptr;
    EoS_ptr EoS_2 = nullptr;
    method_ptr numMethod = nullptr;
    createEoS(EoS_1, type, gamma1, p_ref1);
    createEoS(EoS_2, type, gamma2, p_ref2);
    createNumericalMethod(numMethod, scheme);
    ng = numMethod->ng();
    ODE_ptr source = std::make_shared<ODE_solver>(geometry, ng);
    LevelSet_ptr levelSet = std::make_shared<LevelSet>(ng);

    // Define the variables
    VectofVectDouble u(nCells + 2 * ng, VectOfDouble(nVar));
    VectofVectDouble u2(nCells + 2 * ng, VectOfDouble(nVar));

    VectofVectDouble uPlus1(nCells + 2 * ng, VectOfDouble(nVar));
    VectofVectDouble u2Plus1(nCells + 2 * ng, VectOfDouble(nVar));

    VectofVectDouble uBar(nCells + 2 * ng, VectOfDouble(nVar));
    VectofVectDouble uBar2(nCells + 2 * ng, VectOfDouble(nVar));

    VectOfDouble phi(nCells + 2 * ng);
    VectOfDouble phiPlus1(nCells + 2 * ng);
    VectofVectDouble w_ExactRiemann(nCells + 2 * ng, VectOfDouble(nVar));

    // Construct the FV domain
    double dx;
    VectOfDouble xCells(nCells + 2 * ng);
    ConstructFVDomain(dx, xCells);

    // Set the initial conditions
    setInitialConditions(xCells, u, u2, phi, tStop, w_ExactRiemann, initial, scheme, exactSolution, EoS_1, EoS_2);

    std::ostringstream filename = generate_filename(scheme, initial, geometry, EoS_1, EoS_2);
    std::ofstream output("./" + filename.str() + ".dat");
    std::cout << "Solving " << filename.str() << std::endl;
    int step = 0;
    double dt;
    double dt1;
    double dt2;
    double t = 0.0;
    // Time stepping loop
    do
    {
        // Store the data
        writeToFile(output, t, xCells, u, u2, phi, w_ExactRiemann, EoS_1, EoS_2);

        applyGhostBoundaryConditions(u, u2, phi, EoS_1, EoS_2);

        // Compute time step
        computeDt(u, dx, dt1, step, EoS_1);
        computeDt(u2, dx, dt2, step, EoS_2);
        dt = std::min(dt1, dt2);
        t += dt;
        step++;

        // Update
        levelSet->reinitialisation(phi, dx, x1);
        levelSet->phi_evolution(phi, phiPlus1, dx, dt, u, u2, EoS_1, EoS_2);
        source->sourceUpdate(u, uBar, xCells, dt, EoS_1);
        source->sourceUpdate(u2, uBar2, xCells, dt, EoS_2);
        applyBoundaryConditions(uBar, boundary, EoS_1);
        applyBoundaryConditions(uBar2, boundary, EoS_2);
        numMethod->UpdatewithFluxes(uBar, uPlus1, dx, dt, EoS_1);
        numMethod->UpdatewithFluxes(uBar2, u2Plus1, dx, dt, EoS_2);
        applyBoundaryConditions(uPlus1, boundary, EoS_1);
        applyBoundaryConditions(u2Plus1, boundary, EoS_2);

        u = uPlus1;
        u2 = u2Plus1;
        phi = phiPlus1;

    } while (t < tStop);
    output.close();
    plot(filename, initial);
};

/* -------------------------------------------------------------------*/
/*                                                                    */
/*                  Initial Conditions Routine                        */ // OK
/*                                                                    */
/* -------------------------------------------------------------------*/
void setInitialConditions(const VectOfDouble &xCells, VectofVectDouble &u, VectofVectDouble &u2, VectOfDouble &phi, double &tStop,
                          VectofVectDouble &w_ExactRiemann, std::string &initial, std::string &scheme, bool exactSolution,
                          EoS_ptr &EoS_1, EoS_ptr &EoS_2)
{
    VectOfDouble w_L(nVar), w_R(nVar), u_L(nVar), u_R(nVar);

    double xDiscontinuity;

    if (initial == "toro_1")
    {
        if (((scheme == "Godunov") || (scheme == "HLLC") || (scheme == "HLL") || (scheme == "FVTV")) && (type == "ideal"))
        {
            xDiscontinuity = 0.3;
            tStop = 0.2;
            w_L = {1.0, 0.75, 1.0, 1.4, 0.0};
            w_R = {0.125, 0.0, 0.1, 1.4, 0.0};
        }
        else if (type == "stiffened")
        {
            xDiscontinuity = 0.5;
            tStop = 2e-4;
            gamma = 1.4;
            gamma2 = 1.4;
            w_L = {1500.0, 0.0, 3000.0 * 101325.0};
            w_R = {1000.0, 0.0, 101325.0};
        }
        else
        {
            xDiscontinuity = 0.5;
            tStop = 0.25;
            w_L = {1.0, 0.0, 1.0, 1.4, 0.0};
            w_R = {0.125, 0.0, 0.1, 1.4, 0.0};
        }
    }
    else if (initial == "toro_1flipped")
    {
        xDiscontinuity = 0.5;
        tStop = 0.25;
        w_L = {0.125, 0.0, 0.1, 1.4, 0.0};
        w_R = {1.0, 0.0, 1.0, 1.4, 0.0};
    }
    else if (initial == "toro_2")
    {
        if (type == "stiffened")
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
            w_L = {1.0, -2.0, 0.4, 1.4, 0.0};
            w_R = {1.0, 2.0, 0.4, 1.4, 0.0};
        }
    }
    else if (initial == "toro_3")
    {
        xDiscontinuity = 0.5;
        tStop = 0.012;
        w_L = {1.0, 0.0, 1000.0, 1.4, 0.0};
        w_R = {1.0, 0.0, 0.01, 1.4, 0.0};
    }
    else if (initial == "toro_4")
    {
        if (scheme == "Godunov" || scheme == "HLLC" || (scheme == "HLL") || (scheme == "FVTV"))
        {
            xDiscontinuity = 0.4;
            tStop = 0.035;
            w_L = {5.99924, 19.5975, 460.894, 1.4, 0.0};
            w_R = {5.99242, -6.19633, 46.0950, 1.4, 0.0};
        }
        else
        {
            xDiscontinuity = 0.5;
            tStop = 0.035;
            w_L = {1.0, 0.0, 0.01, 1.4, 0.0};
            w_R = {1.0, 0.0, 100.0, 1.4, 0.0};
        }
    }
    else if (initial == "toro_5")
    {
        if (scheme == "Godunov" || scheme == "HLLC" || (scheme == "HLL") || (scheme == "FVTV"))
        {
            xDiscontinuity = 0.8;
            tStop = 0.012;
            w_L = {1.0, -19.5975, 1000.0, 1.4, 0.0};
            w_R = {1.0, -19.5975, 0.01, 1.4, 0.0};
        }
        else
        {
            xDiscontinuity = 0.5;
            tStop = 0.035;
            w_L = {5.99924, 19.5975, 460.894, 1.4, 0.0};
            w_R = {5.99242, -6.19633, 46.0950, 1.4, 0.0};
        }
    }
    else if (initial == "test_1")
    {
        xDiscontinuity = 0.5;
        tStop = 1.0;
        w_L = {1.0, 0.0, 1.0};
        w_R = {1.0, 0.0, 1.0};
    }
    else if (initial == "test_2")
    {
        xDiscontinuity = 0.5;
        tStop = 1.0;
        w_L = {2.0, 1.0, 3.0};
        w_R = {2.0, 1.0, 3.0};
    }
    else if (initial == "test_3")
    {
        xDiscontinuity = 0.5;
        tStop = 1.0;
        w_L = {1.0, 0.0, 1.0};
        w_R = {0.1, 0.0, 1.0};
    }
    else if (initial == "test_4")
    {
        xDiscontinuity = 0.25;
        tStop = 0.5;
        w_L = {1.0, 1.0, 1.0, gamma, 0.0};
        w_R = {0.1, 1.0, 1.0, gamma2, 0.0};
    }
    else if (initial == "explosion")
    {
        xDiscontinuity = 0.4;
        tStop = 0.25;
        w_L = {1.0, 0.0, 1.0};
        w_R = {0.125, 0.0, 0.1};
    }
    else if (initial == "JWL")
    {
        xDiscontinuity = 0.5;
        tStop = 12e-6;
        w_L = {1700.0, 0.0, 1e12};
        w_R = {1000.0, 0.0, 5e10};
    }
    else if (initial == "TestA")
    {
        xDiscontinuity = 0.505;
        tStop = 0.0007;
        gamma = 1.4;
        gamma2 = 1.2;
        w_L = {1.0, 0.0, 1e5, gamma, 0.0};
        w_R = {0.125, 0.0, 1e4, gamma2, 0.0};
    }
    else if (initial == "TestB")
    {
        tStop = 0.0012;
        VectOfDouble w_M(nVar);
        xDiscontinuity = 0.505;
        double x_left = 0.055;
        gamma = 1.4;
        w_L = {1.3333, 0.3535 * sqrt(1e5), 1.5e5, gamma, 0.0};
        w_M = {1.0, 0.0, 1e5, gamma, 0.0};
        gamma2 = 1.2;
        w_R = {0.1379, 0.0, 1e5, gamma2, 0.0};
        double S_shock = (w_L[0] * w_L[1] - w_M[0] * w_M[1]) / (w_L[0] - w_M[0]);
        double time_to_progress_to_RP = (xDiscontinuity - x_left) / S_shock;
        tStop -= time_to_progress_to_RP;
    }
    else if (initial == "TestC")
    {
        tStop = 0.0017;
        VectOfDouble w_M(nVar);
        xDiscontinuity = 0.505;
        double x_left = 0.055;
        gamma = 1.4;
        w_L = {1.3333, 0.3535 * sqrt(1e5), 1.5e5, gamma, 0.0};
        w_M = {1.0, 0.0, 1e5, gamma, 0.0};
        gamma2 = 1.2;
        w_R = {3.1538, 0.0, 1e5, gamma2, 0.0};
        double S_shock = (w_L[0] * w_L[1] - w_M[0] * w_M[1]) / (w_L[0] - w_M[0]);
        double time_to_progress_to_RP = (xDiscontinuity - x_left) / S_shock;
        tStop -= time_to_progress_to_RP;
    }
    else if (initial == "TestD1")
    {
        tStop = 0.0005;
        VectOfDouble w_M(nVar);
        xDiscontinuity = 0.505;
        double x_left = 0.055;
        gamma = 1.4;
        w_L = {4.3333, 3.2817 * sqrt(1e5), 1.5e6, gamma, 0.0};
        w_M = {1.0, 0.0, 1e5, gamma, 0.0};
        gamma2 = 1.2;
        w_R = {0.1379, 0.0, 1e5, gamma2, 0.0};
        double S_shock = (w_L[0] * w_L[1] - w_M[0] * w_M[1]) / (w_L[0] - w_M[0]);
        double time_to_progress_to_RP = (xDiscontinuity - x_left) / S_shock;
        tStop -= time_to_progress_to_RP;
    }
    else if (initial == "TestD2")
    {
        tStop = 0.0007;
        VectOfDouble w_M(nVar);
        xDiscontinuity = 0.505;
        double x_left = 0.055;
        gamma = 1.4;
        w_L = {4.3333, 3.2817 * sqrt(1e5), 1.5e6, gamma, 0.0};
        w_M = {1.0, 0.0, 1e5, gamma, 0.0};
        gamma2 = 1.2;
        w_R = {3.1538, 0.0, 1e5, gamma2, 0.0};
        double S_shock = (w_L[0] * w_L[1] - w_M[0] * w_M[1]) / (w_L[0] - w_M[0]);
        double time_to_progress_to_RP = (xDiscontinuity - x_left) / S_shock;
        tStop -= time_to_progress_to_RP;
    }
    else if (initial == "CDS")
    {
        tStop = 237.44e-6;
        xDiscontinuity = 0.7;
        w_L = {1000.0, 0.0, 1e9};
        w_R = {50.0, 0.0, 1e5};
    }
    else if (initial == "stationary1")
    {
        xDiscontinuity = 0.5;
        tStop = 0.2;
        w_L = {1.0, 0.0, 1.0, 1.4, 0.0};
        w_R = {0.5, 0.0, 1.0, 1.4, 0.0};
    }
    else if (initial == "stationary2")
    {
        xDiscontinuity = 0.5;
        tStop = 0.2;
        w_L = {1.0, 0.0, 1.0, 1.4, 0.0};
        w_R = {0.5, 0.0, 1.0, 1.67, 0.0};
    }
    else if (initial == "moving1")
    {
        xDiscontinuity = 0.5;
        tStop = 0.5;
        w_L = {1.0, 0.5, 1.0, 1.4, 0.0};
        w_R = {0.5, 0.5, 1.0, 1.4, 0.0};
    }
    else if (initial == "moving2")
    {
        xDiscontinuity = 0.5;
        tStop = 0.5;
        w_L = {1.0, 0.5, 1.0, 1.4, 0.0};
        w_R = {0.5, 0.5, 1.0, 1.67, 0.0};
    }
    else if (initial == "moving3")
    {
        xDiscontinuity = 0.5;
        tStop = 0.5;
        w_L = {1.0, -0.5, 1.0, 1.4, 0.0};
        w_R = {0.5, -0.5, 1.0, 1.67, 0.0};
    }
    else if (initial == "FedkiwA")
    {
        xDiscontinuity = 0.5;
        tStop = 0.0007;
        w_L = {1.0, 0.0, 1e5, 1.4, 0.0};
        w_R = {0.125, 0.0, 1e4, 1.2, 0.0};
    }
    else if (initial == "FedkiwB")
    {
        xDiscontinuity = 0.5;
        tStop = 0.0012;
        w_L = {1.333, 0.3535 * sqrt(1e5), 1.5e5, 1.4, 0.0};
        w_R = {0.1379, 0.0, 1e5, 1.67, 0.0};
    }
    else if (initial == "heliumSlab")
    {
        xDiscontinuity = 0.5;
        tStop = 0.3;
        w_L = {1.3765, 0.3948, 1.57};
        w_R = {0.138, 0.0, 1.0};
    }
    u_L = EoS_1->primitiveToConservative(w_L);
    u_R = EoS_2->primitiveToConservative(w_R);

    for (int i = 0; i < u.size(); i++)
    {
        double x = xCells[i];
        // if(x <= xDiscontinuity){
        // 	u[i] = u_L;
        // }
        // else{
        // 	u[i] = u_R;
        // }
        if (initial == "FedkiwB")
        {
            VectOfDouble w_L2(nVar), u_L2(nVar);
            w_L2 = {1.0, 0.0, 1e5, 1.4, 0.0};
            u_L2 = EoS_1->primitiveToConservative(w_L2);
            if (x < 0.05)
            {
                u[i] = u_L;
            }
            else
            {
                u[i] = u_L2;
            }
        }
        else if (initial == "heliumSlab")
        {
            VectOfDouble w_L2(nVar), u_L2(nVar);
            w_L2 = {1.0, 0.0, 1.0};
            u_L2 = EoS_1->primitiveToConservative(w_L2);
            if (x < 0.25)
            {
                u[i] = u_L;
            }
            else
            {
                u[i] = u_L2;
            }
        }
        else
        {
            u[i] = u_L;
        }

        u2[i] = u_R;

        if (initial == "heliumSlab")
        {
            phi[i] = 0.1 - abs(x - 0.5);
        }
        else
        {
            phi[i] = x - xDiscontinuity;
        }
        if (exactSolution)
        {
            if (initial == "FedkiwB")
            {
                VectOfDouble w_L2(nVar);
                w_L2 = {1.0, 0.0, 1e5, 1.4, 0.0};
                double S_shock = (w_L[0] * w_L[1] - w_L2[0] * w_L2[1]) / (w_L[0] - w_L2[0]);
                double time_to_progress_to_RP = (xDiscontinuity - 0.05) / S_shock;
                double time = tStop - time_to_progress_to_RP;
                w_ExactRiemann[i] = exactRiemannSolver(x, xDiscontinuity, time, w_L, w_R, EoS_1, EoS_1);
            }
            else
            {
                w_ExactRiemann[i] = exactRiemannSolver(x, xDiscontinuity, tStop, w_L, w_R, EoS_1, EoS_2);
            }
        }
    }
    return;
}

/* -------------------------------------------------------------------*/
/*                                                                    */
/*                      PrintUtility Routine                          */
/*                                                                    */
/* -------------------------------------------------------------------*/
std::ostringstream generate_filename(std::string scheme, std::string initial, std::string &geometry, EoS_ptr &EoS_1, EoS_ptr &EoS_2)
{
    std::ostringstream filename;
    filename << "Euler_" << scheme << "_" << initial
             << "_n=" << nCells
             << std::fixed << std::setprecision(3) << "_t=" << tStop
             << "_boundary_" << boundary[0] << boundary[1] << boundary[2]
             << "_" << geometry
             << "_" << EoS_1->info()
             << "_" << EoS_2->info();
    return filename;
};

void writeToFile(std::ofstream &output, double t, VectOfDouble xCells, VectofVectDouble &u, VectofVectDouble &u2, VectOfDouble &phi, VectofVectDouble &w_ExactRiemann, EoS_ptr &EoS_1, EoS_ptr &EoS_2)
{
    for (int i = ng; i < nCells + ng; i++)
    {
        double x = xCells[i];
        VectOfDouble w = EoS_1->conservativeToPrimitive(u[i]);
        VectOfDouble w2 = EoS_2->conservativeToPrimitive(u2[i]);
        output << t << " " << x << " " << w[0] << " " << w[1] << " " << w[2] << " "
               << w2[0] << " " << w2[1] << " " << w2[2] << " " << phi[i] << " " << w_ExactRiemann[i][0] << " " << w_ExactRiemann[i][1] << " " << w_ExactRiemann[i][2] << std::endl;
    };
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
        a[i] = abs(w[1]) + EoS->c_s(w);
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

void applyBoundaryConditions(VectofVectDouble &u, std::array<std::string, nVar> boundary, EoS_ptr &EoS)
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
        };
    };
};

void applyGhostBoundaryConditions(VectofVectDouble &u, VectofVectDouble &u2, VectOfDouble &phi, EoS_ptr &EoS_1, EoS_ptr &EoS_2)
{
    VectOfDouble interface_vec;
    VectOfDouble w(nVar), w2(nVar), new_u(nVar), wI(nVar), w2I(nVar), wLI(nVar), w2LI(nVar), wRI(nVar), w2RI(nVar);
    for (int i = ng; i < nCells + ng; i++)
    {
        if (phi[i] * phi[i + 1] < 0.0)
        {
            interface_vec.push_back(i);
        };
    };
    if (interface_vec.size() == 1)
    {
        int interface = interface_vec[0];
        wI = EoS_1->conservativeToPrimitive(u[interface]);
        w2I = EoS_2->conservativeToPrimitive(u2[interface + 1]);
        // Copy pressure and velocity into ghost cells
        for (int i = ng; i < nCells + ng; i++)
        {
            w = EoS_1->conservativeToPrimitive(u[i]);
            w2 = EoS_2->conservativeToPrimitive(u2[i]);
            if (i <= interface)
            {
                w2[1] = w[1];
                w2[2] = w[2];
                w2[0] = w2I[0] * pow((w2[2] / w2I[2]), 1.0 / (w2I[3]));
                new_u = EoS_2->primitiveToConservative(w2);
                u2[i][0] = new_u[0];
                u2[i][1] = new_u[1];
                u2[i][2] = new_u[2];
            }
            else
            {
                w[1] = w2[1];
                w[2] = w2[2];
                w[0] = wI[0] * pow((w[2] / wI[2]), 1.0 / (wI[3]));
                new_u = EoS_1->primitiveToConservative(w);
                u[i][0] = new_u[0];
                u[i][1] = new_u[1];
                u[i][2] = new_u[2];
            };
        };
    }
    else if (interface_vec.size() == 2)
    {
        int interface1 = interface_vec[0];
        int interface2 = interface_vec[1];
        wLI = EoS_1->conservativeToPrimitive(u[interface1]);
        w2LI = EoS_2->conservativeToPrimitive(u2[interface1 + 1]);
        w2RI = EoS_2->conservativeToPrimitive(u2[interface2]);
        wRI = EoS_1->conservativeToPrimitive(u[interface2 + 1]);
        double gamma_1 = EoS_1->gamma();
        double gamma_2 = EoS_2->gamma();
        for (int i = ng; i < nCells + ng; i++)
        {
            w = EoS_1->conservativeToPrimitive(u[i]);
            w2 = EoS_2->conservativeToPrimitive(u2[i]);
            if (((phi[i] > 0) && (phi[i] > phi[i - 1])) || ((phi[i] < 0) && (phi[i] < phi[i - 1])))
            {
                if (phi[i] > 0)
                { // set ghost zone of material 1 after interface 1
                    w[1] = w2[1];
                    w[2] = w2[2];
                    w[0] = wLI[0] * pow((w[2] / wLI[2]), 1.0 / gamma_1);
                    new_u = EoS_1->primitiveToConservative(w);
                    u[i][0] = new_u[0];
                    u[i][1] = new_u[1];
                    u[i][2] = new_u[2];
                }
                else
                { // set ghost zone of material 2 after interface 2
                    w2[1] = w[1];
                    w2[2] = w[2];
                    w2[0] = w2RI[0] * pow((w2[2] / w2RI[2]), 1.0 / gamma_2);
                    new_u = EoS_2->primitiveToConservative(w2);
                    u2[i][0] = new_u[0];
                    u2[i][1] = new_u[1];
                    u2[i][2] = new_u[2];
                }
            }
            else
            {
                if (phi[i] > 0)
                { // set ghost zone of material 1 after interface 2
                    w[1] = w2[1];
                    w[2] = w2[2];
                    w[0] = wRI[0] * pow((w[2] / wRI[2]), 1.0 / gamma_1);
                    new_u = EoS_1->primitiveToConservative(w);
                    u[i][0] = new_u[0];
                    u[i][1] = new_u[1];
                    u[i][2] = new_u[2];
                }
                else
                { // set ghost zone of material 2 after interface 1
                    w2[1] = w[1];
                    w2[2] = w[2];
                    w2[0] = w2LI[0] * pow((w2[2] / w2LI[2]), 1.0 / gamma_2);
                    new_u = EoS_2->primitiveToConservative(w2);
                    u2[i][0] = new_u[0];
                    u2[i][1] = new_u[1];
                    u2[i][2] = new_u[2];
                }
            }
        }
    }
};
