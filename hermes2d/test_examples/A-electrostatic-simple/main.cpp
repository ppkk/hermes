#define HERMES_REPORT_ALL
#include "definitions.h"
#include "exceptions.h"

// This example shows how to solve a simple PDE that describes stationary
// heat transfer in an object consisting of two materials (aluminum and
// copper). The object is heated by constant volumetric heat sources
// (generated, for example, by a DC electric current). The temperature
// on the boundary is fixed. We will learn how to:
//
//   - load the mesh,
//   - perform initial refinements,
//   - create a H1 space over the mesh,
//   - define weak formulation,
//   - initialize matrix solver,
//   - assemble and solve the matrix system,
//   - output the solution and element orders in VTK format
//     (to be visualized, e.g., using Paraview),
//   - visualize the solution using Hermes' native OpenGL-based functionality.
//
// PDE: Poisson equation -div(LAMBDA grad u) - VOLUME_HEAT_SRC = 0.
//
// Boundary conditions: Dirichlet u(x, y) = FIXED_BDY_TEMP on the boundary.
//
// Geometry: L-Shape domain (see file domain.mesh).
//
// The following parameters can be changed:

const bool HERMES_VISUALIZATION = true;           // Set to "false" to suppress Hermes OpenGL visualization.
const bool VTK_VISUALIZATION = true;              // Set to "true" to enable VTK output.
const int P_INIT = 1;                              // Uniform polynomial degree of mesh elements.
const int INIT_REF_NUM = 0;                       // Number of initial uniform mesh refinements.

// Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
Hermes::MatrixSolverType matrix_solver_type = Hermes::SOLVER_UMFPACK;

// Problem parameters.
const double PERM_AIR = 8.854e-12;
const double POTENTIAL0 = 0.0;
const double POTENTIAL1 = 10.0;

int main(int argc, char* argv[])
{
  // Time measurement.
  Hermes::TimePeriod cpu_time;

  // Load the mesh.
  Hermes::Hermes2D::Mesh mesh;
  Hermes::Hermes2D::MeshReaderH2DXML mloader;
  mloader.load("aa-electrostatic-ctverec-triangle-linear.xml", &mesh);

  // Perform initial mesh refinements (optional).
  for (int i = 0; i < INIT_REF_NUM; i++)
    mesh.refine_all_elements();

  // Initialize the weak formulation.
  CustomWeakFormPoisson wf("0", new Hermes::Hermes1DFunction<double>(PERM_AIR));

  // Initialize essential boundary conditions.
  Hermes::Hermes2D::DefaultEssentialBCConst<double> bc_essential_outer(Hermes::vector<std::string>("1", "2", "3", "4"),
    POTENTIAL0);
  Hermes::Hermes2D::DefaultEssentialBCConst<double> bc_essential_inner(Hermes::vector<std::string>("5", "6", "7", "8"),
    POTENTIAL1);
  Hermes::Hermes2D::EssentialBCs<double> bcs(Hermes::vector<Hermes::Hermes2D::EssentialBoundaryCondition<double>* >(&bc_essential_inner, &bc_essential_outer));

  // Create an H1 space with default shapeset.
  Hermes::Hermes2D::H1Space<double> space(&mesh, &bcs, P_INIT);
  int ndof = space.get_num_dofs();
  info("ndof = %d", ndof);

  // Initialize the FE problem.
  Hermes::Hermes2D::DiscreteProblem<double> dp(&wf, &space);

  // Initial coefficient vector for the Newton's method.
  double* coeff_vec = new double[ndof];
  memset(coeff_vec, 0, ndof*sizeof(double));

  // Perform Newton's iteration and translate the resulting coefficient vector into a Solution.
  Hermes::Hermes2D::Solution<double> sln;
  Hermes::Hermes2D::NewtonSolver<double> newton(&dp, matrix_solver_type);

  newton.set_verbose_output(true);

  try
  {
    newton.solve(coeff_vec);
  }
  catch(Hermes::Exceptions::Exception e)
  {
    e.printMsg();
    error("Newton's iteration failed.");
  }
  
  Hermes::Hermes2D::Solution<double>::vector_to_solution(newton.get_sln_vector(), &space, &sln);

  cpu_time.tick();
  printf("Duration %lf\n", cpu_time.accumulated());

  // VTK output.
  if (VTK_VISUALIZATION)
  {
    // Output solution in VTK format.
    Hermes::Hermes2D::Views::Linearizer lin;
    bool mode_3D = true;
    lin.save_solution_vtk(&sln, "sln.vtk", "Temperature", mode_3D);
    info("Solution in VTK format saved to file %s.", "sln.vtk");

    // Output mesh and element orders in VTK format.
    Hermes::Hermes2D::Views::Orderizer ord;
    ord.save_orders_vtk(&space, "ord.vtk");
    info("Element orders in VTK format saved to file %s.", "ord.vtk");
  }

  // Visualize the solution.
  if (HERMES_VISUALIZATION)
  {
    Hermes::Hermes2D::Views::ScalarView view("Solution", new Hermes::Hermes2D::Views::WinGeom(0, 0, 440, 350));
    // Hermes uses adaptive FEM to approximate higher-order FE solutions with linear
    // triangles for OpenGL. The second parameter of View::show() sets the error
    // tolerance for that. Options are HERMES_EPS_LOW, HERMES_EPS_NORMAL (default),
    // HERMES_EPS_HIGH and HERMES_EPS_VERYHIGH. The size of the graphics file grows
    // considerably with more accurate representation, so use it wisely.
    view.show(&sln, Hermes::Hermes2D::Views::HERMES_EPS_HIGH);
    Hermes::Hermes2D::Views::View::wait();
  }

  // Clean up.
  delete [] coeff_vec;

  return 0;
}
