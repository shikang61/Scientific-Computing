#include "euler2DSolver.H"

#define DEBUG(message) std::cout << "You are at line " << __LINE__ << " doing: " << #message << std::endl

// Global parameters
int Nx = 100;
int Ny = 100;
double x0 = 0.0;
double x1 = 1.0;
double yBottom = -1.0;
double yTop = 1.0;
double C = 0.8;
int ng;
double tStop;
std::array<std::string, nVar> boundary = {"T", "R", "R", "T"}; // T = transmissive, R = reflective, P = periodic

// Plotting parameters
bool animate = true;
bool plotGraph = true;

// Specify gas parameters
std::string type = "ideal"; // stiffened, ideal
double gamma = 1.4;
double p_ref = 0.0;

std::pair<std::string, int> sliceX = {"X", 50};
std::pair<std::string, int> sliceY = {"Y", 50};
std::pair<std::string, int> plot2D = {"2D", 0};

int main()
{
    // solveEuler("SLIC", "toro_1d", plot2D, "cartesian");
    // solveEuler("SLIC", "explosion", plot2D, "cylindrical");
    // solveEuler("SLIC", "explosion", sliceY, "cylindrical");
    // solveEuler("SLIC", "2D_RP_1", plot2D, "cylindrical");
    // solveEuler("SLIC", "2D_RP_2", plot2D, "cylindrical");
    solveEuler("SLIC", "explosion", plot2D, "spherical");
    solveEuler("SLIC", "explosion", sliceY, "spherical");
    return 0;
};

/* -------------------------------------------------------------------*/
/*                                                                    */
/*                              Solver                                */
/*                                                                    */
/* -------------------------------------------------------------------*/

void solveEuler(std::string scheme, std::string initial, std::pair<std::string, int> &slice, std::string geometry)
{
    EoS_ptr EoS = nullptr;
    method_ptr numMethod_x = nullptr;
    method_ptr numMethod_y = nullptr;
    createEoS(EoS, type, gamma, p_ref);
    createNumericalMethod(numMethod_x, scheme, "x");
    createNumericalMethod(numMethod_y, scheme, "y");
    ng = numMethod_x->ng();
    ODE_ptr source = std::make_shared<ODE_solver>(geometry, ng);

    // Define the variables
    Vect2D u(Nx + 2 * ng, Vect1D(Ny + 2 * ng, VectOfDouble(nVar)));
    Vect2D uBar(Nx + 2 * ng, Vect1D(Ny + 2 * ng, VectOfDouble(nVar)));
    Vect2D uPlus1(Nx + 2 * ng, Vect1D(Ny + 2 * ng, VectOfDouble(nVar)));
    Vect2D w_ExactRiemann(Nx + 2 * ng, Vect1D(Ny + 2 * ng, VectOfDouble(nVar)));

    // Construct the FV domain
    double dx, dy;
    Vect2D xyCells(Nx + 2 * ng, Vect1D(Ny + 2 * ng, VectOfDouble(2)));
    ConstructFVDomain(dx, dy, xyCells);

    // Set the initial conditions
    setInitialConditions(xyCells, u, tStop, initial, numMethod_x, geometry, EoS);

    std::ostringstream filename = generate_filename(numMethod_x, initial, slice, geometry, EoS);
    std::ofstream output("./" + filename.str() + ".dat");

    std::cout << "Solving " << filename.str() << std::endl;

    int step = 0;
    double dt;
    double t = 0.0;
    // Time stepping loop
    do
    {
        // Store the data
        writeToFile(output, t, xyCells, u, slice, EoS);

        // Compute time step
        computeDt(u, dx, dy, dt, step, EoS);
        t += dt;
        step++;

        // Update
        source->sourceUpdate(u, uBar, xyCells, dt, EoS);
        applyBoundaryConditions(uBar, boundary, Nx, Ny, geometry);
        numMethod_x->UpdatewithFluxes(uBar, uBar, dx, dy, dt, EoS);
        applyBoundaryConditions(uBar, boundary, Nx, Ny, geometry);
        numMethod_y->UpdatewithFluxes(uBar, uPlus1, dx, dy, dt, EoS);
        applyBoundaryConditions(uPlus1, boundary, Nx, Ny, geometry);

        u = uPlus1;
    } while (t < tStop);
    output.close();

    if (plotGraph)
    {
        plot(filename, initial, slice);
    };

    std::cout << "Finished solving " << filename.str() << std::endl;
    std::cout << "-----------------------------------" << std::endl;
};

/* -------------------------------------------------------------------*/
/*                                                                    */
/*                      PrintUtility Routine                          */ // OK
/*                                                                    */
/* -------------------------------------------------------------------*/
std::ostringstream generate_filename(method_ptr &numMethod, std::string &initial, std::pair<std::string, int> &slice, std::string &geometry, EoS_ptr &EoS)
{
    std::ostringstream filename;
    filename << "Euler2D_" << numMethod->name() << "_" << initial
             << "_" << geometry
             << "_slice" << slice.first << slice.second
             << "_Nx=" << Nx << "_Ny=" << Ny
             << std::fixed << std::setprecision(3) << "_t=" << tStop
             << std::fixed << std::setprecision(1) << "_C=" << C
             << "_" << boundary[0] << boundary[1] << boundary[2] << boundary[3];
    return filename;
};

void writeToFile(std::ofstream &output, double &t, Vect2D &xyCells, Vect2D &u, std::pair<std::string, int> &slice, EoS_ptr &EoS)
{
    if (slice.first == "X")
    {
        int j = slice.second + ng - 1;
        for (int i = ng; i < Nx + ng; i++)
        {
            double x = xyCells[i][j][0];
            VectOfDouble w = EoS->conservativeToPrimitive(u[i][j]);
            output << t << " " << x << " " << w[RHO] << " " << w[VX] << " " << w[VY] << " " << w[P] << std::endl;
        };
    }
    else if (slice.first == "Y")
    {
        int i = slice.second + ng - 1;
        for (int j = ng; j < Ny + ng; j++)
        {
            double y = xyCells[i][j][1];
            VectOfDouble w = EoS->conservativeToPrimitive(u[i][j]);
            output << t << " " << y << " " << w[RHO] << " " << w[VX] << " " << w[VY] << " " << w[P] << std::endl;
        };
    }
    else if (slice.first == "2D")
    {
        for (int i = ng; i < Nx + ng; i++)
        {
            for (int j = ng; j < Ny + ng; j++)
            {
                double x = xyCells[i][j][0];
                double y = xyCells[i][j][1];
                VectOfDouble w = EoS->conservativeToPrimitive(u[i][j]);
                output << t << " " << x << " " << y << " " << w[RHO] << " " << w[VX] << " " << w[VY] << " " << w[P] << std::endl;
            };
            output << std::endl;
        };
    };
    output << std::endl
           << std::endl;
};

void plot(std::ostringstream &filename, std::string &initial, std::pair<std::string, int> &slice)
{
    if (slice.first == "2D")
    {
        if (animate)
        {
            system(("gnuplot -c gif2D.gp " + filename.str()).c_str());
            // system(("python gif3D.py " + filename.str()).c_str());
        };
        system(("gnuplot -c png2D.gp " + filename.str()).c_str());
        system(("gnuplot -c png3D.gp " + filename.str()).c_str());
        // system(("python png3D.py " + filename.str()).c_str());
    }
    else
    {
        if (animate)
        {
            system(("gnuplot -c gif1D.gp " + filename.str()).c_str());
        };
        if (initial == "isolated_1" || initial == "isolated_2" || initial == "woodward")
        {
            system(("gnuplot -c png_isolated.gp " + filename.str()).c_str());
        }
        else
        {
            system(("gnuplot -c png1D.gp " + filename.str()).c_str());
        };
    };
    remove(("./" + filename.str() + ".dat").c_str());
};

/* -------------------------------------------------------------------*/
/*                                                                    */
/*                  Domain Construction Routine                       */ // OK
/*                                                                    */
/* -------------------------------------------------------------------*/

void ConstructFVDomain(double &dx, double &dy, Vect2D &xyCells)
{
    dx = (x1 - x0) / Nx;
    dy = (yTop - yBottom) / Ny;
    for (int j = 0; j < Ny + 2 * ng; j++)
    {
        for (int i = 0; i < Nx + 2 * ng; i++)
        {
            xyCells[i][j][0] = x0 + dx * (i - ng + 0.5);
            xyCells[i][j][1] = yBottom + dy * (j - ng + 0.5);
        }
    }
};

/* -------------------------------------------------------------------*/
/*                                                                    */
/*                  Initial Conditions Routine                        */ // OK
/*                                                                    */
/* -------------------------------------------------------------------*/
void setInitialConditions(const Vect2D &xyCells, Vect2D &u, double &tStop, std::string &initial, method_ptr &numMethod, std::string &geometry, EoS_ptr &EoS)
{
    VectOfDouble w_Q1(nVar), w_Q2(nVar), w_Q3(nVar), w_Q4(nVar), u_Q1(nVar), u_Q2(nVar), u_Q3(nVar), u_Q4(nVar);
    ;
    double xDiscontinuity = 0.0;
    double yDiscontinuity = 0.0;

    /*  Q2 | Q1
       ----|----
        Q3 | Q4
    */
    if (initial == "toro_1x")
    {
        xDiscontinuity = 0.5;
        yDiscontinuity = 0.5;
        tStop = 0.25;
        w_Q2 = {1.0, 0.0, 0.0, 1.0};
        w_Q3 = {1.0, 0.0, 0.0, 1.0};
        w_Q1 = {0.125, 0.0, 0.0, 0.1};
        w_Q4 = {0.125, 0.0, 0.0, 0.1};
    }
    else if (initial == "toro_3x")
    {
        xDiscontinuity = 0.5;
        yDiscontinuity = 0.5;
        tStop = 0.012;
        w_Q2 = {1.0, 0.0, 0.0, 1000.0};
        w_Q3 = {1.0, 0.0, 0.0, 1000.0};
        w_Q1 = {0.125, 0.0, 0.0, 0.01};
        w_Q4 = {0.125, 0.0, 0.0, 0.01};
    }
    else if (initial == "toro_1y")
    {
        xDiscontinuity = 0.5;
        yDiscontinuity = 0.5;
        tStop = 0.25;
        w_Q1 = {1.0, 0.0, 0.0, 1.0};
        w_Q2 = {1.0, 0.0, 0.0, 1.0};
        w_Q3 = {0.125, 0.0, 0.0, 0.1};
        w_Q4 = {0.125, 0.0, 0.0, 0.1};
    }
    else if (initial == "toro_3y")
    {
        xDiscontinuity = 0.5;
        yDiscontinuity = 0.5;
        tStop = 0.012;
        w_Q1 = {1.0, 0.0, 0.0, 1000.0};
        w_Q2 = {1.0, 0.0, 0.0, 1000.0};
        w_Q3 = {0.125, 0.0, 0.0, 0.01};
        w_Q4 = {0.125, 0.0, 0.0, 0.01};
    }
    else if (initial == "toro_1d" || initial == "explosion")
    {
        tStop = 0.25;
        w_Q1 = {1.0, 0.0, 0.0, 1.0};
        w_Q2 = {0.125, 0.0, 0.0, 0.1};
        w_Q3 = {0.0, 0.0, 0.0, 0.0};
        w_Q4 = {0.0, 0.0, 0.0, 0.0};
    }
    else if (initial == "2D_RP_1")
    {
        xDiscontinuity = 0.5;
        yDiscontinuity = 0.5;
        tStop = 0.3;
        w_Q1 = {0.5313, 0.0, 0.0, 0.4};
        w_Q2 = {1, 0.7276, 0.0, 1};
        w_Q3 = {0.8, 0.0, 0.0, 1.0};
        w_Q4 = {1.0, 0.0, 0.7276, 1.0};
    }
    else if (initial == "2D_RP_2")
    {
        xDiscontinuity = 0.5;
        yDiscontinuity = 0.5;
        tStop = 0.3;
        w_Q1 = {1.0, 0.75, -0.5, 1.0};
        w_Q2 = {2.0, 0.75, 0.5, 1.0};
        w_Q3 = {1.0, -0.75, 0.5, 1.0};
        w_Q4 = {3.0, -0.75, -0.5, 1.0};
    }
    else if (initial == "inplosion")
    {
        tStop = 0.25;
        w_Q1 = {0.125, 0.0, 0.0, 0.1};
        w_Q2 = {1.0, 0.0, 0.0, 1.0};
        w_Q3 = {0.0, 0.0, 0.0, 0.0};
        w_Q4 = {0.0, 0.0, 0.0, 0.0};
    };

    u_Q1 = EoS->primitiveToConservative(w_Q1);
    u_Q2 = EoS->primitiveToConservative(w_Q2);
    u_Q3 = EoS->primitiveToConservative(w_Q3);
    u_Q4 = EoS->primitiveToConservative(w_Q4);
    double x, y;
    for (int j = 0; j < Ny + 2 * ng; j++)
    {
        for (int i = 0; i < Nx + 2 * ng; i++)
        {
            x = xyCells[i][j][0];
            y = xyCells[i][j][1];
            if (initial == "toro_1d")
            {
                if (x <= y)
                {
                    u[i][j] = u_Q1;
                }
                else
                {
                    u[i][j] = u_Q2;
                }
            }
            else if (initial == "explosion")
            {
                if (geometry == "cylindrical")
                {
                    if ((x - 1) * (x - 1) + (y - 1) * (y - 1) <= 0.4 * 0.4)
                    {
                        u[i][j] = u_Q1;
                    }
                    else
                    {
                        u[i][j] = u_Q2;
                    }
                }
                else if (geometry == "spherical")
                {
                    if (x * x + y * y <= 0.4 * 0.4)
                    {
                        u[i][j] = u_Q1;
                    }
                    else
                    {
                        u[i][j] = u_Q2;
                    }
                }
            }
            else
            {
                if (x >= xDiscontinuity && y >= yDiscontinuity)
                {
                    u[i][j] = u_Q1;
                }
                else if (x < xDiscontinuity && y >= yDiscontinuity)
                {
                    u[i][j] = u_Q2;
                }
                else if (x < xDiscontinuity && y < yDiscontinuity)
                {
                    u[i][j] = u_Q3;
                }
                else if (x >= xDiscontinuity && y < yDiscontinuity)
                {
                    u[i][j] = u_Q4;
                };
            };
        };
    };
};

/* -------------------------------------------------------------------*/
/*                                                                    */
/*                  Compute Time Step Routine                         */ // OK
/*                                                                    */
/* -------------------------------------------------------------------*/
double findMax(const Vect1D &a)
{
    double maxVal = 0.0;
    for (VectOfDouble row : a)
    {
        if (!row.empty())
        { // Check if row is non-empty
            double maxInRow = *std::max_element(row.begin(), row.end());
            maxVal = std::max(maxVal, maxInRow);
        }
    }
    return maxVal;
}

void computeDt(Vect2D &u, double &dx, double &dy, double &dt, int step, EoS_ptr &EoS)
{
    Vect1D a(Nx, VectOfDouble(Ny));
    double a_max;
    for (int i = 0; i < Nx; i++)
    {
        for (int j = 0; j < Ny; j++)
        {
            VectOfDouble w = EoS->conservativeToPrimitive(u[i + ng][j + ng]);
            a[i][j] = sqrt(w[VX] * w[VX] + w[VY] * w[VY]) + EoS->c_s(w);
        };
    };
    a_max = findMax(a);
    if (step < 10)
    {
        dt = 0.5 * C * std::min(dx, dy) / a_max;
    }
    else
    {
        dt = C * std::min(dx, dy) / a_max;
    }
}

/* -------------------------------------------------------------------*/
/*                                                                    */
/*                 Boundary Conditions Routine                        */ // OK
/*                                                                    */
/* -------------------------------------------------------------------*/

void applyBoundaryConditions(Vect2D &u, std::array<std::string, nVar> &boundary, int &Nx, int &Ny, const std::string &geometry)
{
    for (int n = 0; n < nVar; n++)
    {
        if (boundary[n] == "T")
        {
            for (int ig = 0; ig < ng; ig++)
            {
                // Left and right boundaries
                for (int j = ng; j < Ny + ng; j++)
                {
                    u[ig][j][n] = u[2 * ng - 1 - ig][j][n];
                    u[Nx + ng + ig][j][n] = u[Nx + 1 - ig][j][n];
                };
                // Top and bottom boundaries
                for (int i = ng; i < Nx + ng; i++)
                {
                    u[i][ig][n] = u[i][2 * ng - 1 - ig][n];
                    u[i][Ny + ng + ig][n] = u[i][Ny + 1 - ig][n];
                };
            };
        }
        else if (boundary[n] == "R")
        {
            for (int ig = 0; ig < ng; ig++)
            {
                // Left and right boundaries
                for (int j = 0; j < Ny + 2 * ng; j++)
                {
                    if (n == VX)
                    {
                        u[ig][j][n] = -u[2 * ng - 1 - ig][j][n];
                        u[Nx + ng + ig][j][n] = -u[Nx + 1 - ig][j][n];
                    }
                    else if (n == VY)
                    {
                        u[ig][j][n] = u[2 * ng - 1 - ig][j][n];
                        u[Nx + ng + ig][j][n] = u[Nx + 1 - ig][j][n];
                    };
                };
                // Top and bottom boundaries
                for (int i = 0; i < Nx + 2 * ng; i++)
                {
                    if (n == VY)
                    {
                        u[i][ig][n] = -u[i][2 * ng - 1 - ig][n];
                        u[i][Ny + ng + ig][n] = -u[i][Ny + 1 - ig][n];
                    }
                    else if (n == VX)
                    {
                        u[i][ig][n] = u[i][2 * ng - 1 - ig][n];
                        u[i][Ny + ng + ig][n] = u[i][Ny + 1 - ig][n];
                    };
                };
            };
        }
        else if (boundary[n] == "P")
        {
            for (int ig = 0; ig < ng; ig++)
            {
                for (int j = 0; j < Ny + 2 * ng; j++)
                {
                    u[ig][j][n] = u[Nx + ig][j][n];
                    u[Nx + ng + ig][j][n] = u[ng + ig][j][n];
                }
                for (int i = 0; i < Nx + 2 * ng; i++)
                {
                    u[i][ig][n] = u[i][Ny + ig][n];
                    u[i][Ny + ng + ig][n] = u[i][ng + ig][n];
                }
            }
        };
    };
};
