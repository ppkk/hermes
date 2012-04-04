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

using namespace Hermes::Hermes2D::RefinementSelectors;

const bool HERMES_VISUALIZATION = true;           // Set to "false" to suppress Hermes OpenGL visualization.
const bool VTK_VISUALIZATION = true;              // Set to "true" to enable VTK output.
const int P_INIT = 2;                              // Uniform polynomial degree of mesh elements.
const int INIT_REF_NUM = 0;                       // Number of initial uniform mesh refinements.

const double THRESHOLD = 0.3;                     // This is a quantitative parameter of the adapt(...) function and
const int STRATEGY = 0;                           // Adaptive strategy:
//const CandList CAND_LIST = H2D_HP_ANISO;          // Predefined list of element refinement candidates. Possible values are
//const CandList CAND_LIST = H2D_H_ANISO;
const int MESH_REGULARITY = -1;                   // Maximum allowed level of hanging nodes:
const double CONV_EXP = 1.0;                      // Default value is 1.0. This parameter influences the selection of


// Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
Hermes::MatrixSolverType matrix_solver_type = Hermes::SOLVER_UMFPACK;

// Problem parameters.
const double PERM_AIR = 8.854e-12;
const double POTENTIAL0 = 0.0;
const double POTENTIAL1 = 10.0;

const int STEPS = 3;

using namespace Hermes::Hermes2D;

Space<double>* spaces[STEPS + 1];
Solution<double>* solutions[STEPS + 1];
Space<double>* refSpaces[STEPS + 1];
Solution<double>* refSolutions[STEPS + 1];

Views::OrderView order_view("Space");

void adaptiveStep(int step, WeakForm<double>* wf, Selector<double>* selector)
{
    refSpaces[step - 1] = Space<double>::construct_refined_space(spaces[step - 1]);
    refSolutions[step - 1] = new Solution<double>(refSpaces[step - 1]->get_mesh());
    int ndof = refSpaces[step - 1]->get_num_dofs();

    order_view.show(refSpaces[step - 1]);
    order_view.wait_for_keypress();

    // Initialize the FE problem.
    DiscreteProblem<double> dp(wf, refSpaces[step - 1]);

    // Initial coefficient vector for the Newton's method.
    double* coeff_vec = new double[ndof];
    memset(coeff_vec, 0, ndof*sizeof(double));

    // Perform Newton's iteration and translate the resulting coefficient vector into a Solution.
    NewtonSolver<double> newton(&dp, matrix_solver_type);

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

    Solution<double>::vector_to_solution(newton.get_sln_vector(), refSpaces[step - 1], refSolutions[step - 1]);

    OGProjection<double>::project_global(spaces[step - 1], refSolutions[step - 1], solutions[step - 1], matrix_solver_type);

    Mesh* newMesh = new Mesh();
    newMesh->copy(spaces[step - 1]->get_mesh());
    spaces[step] = spaces[step - 1]->dup(newMesh);
    solutions[step] = new Solution<double>(newMesh);

    Adapt<double>* adaptivity = new Adapt<double>(spaces[step]);
    double err_est_rel = adaptivity->calc_err_est(solutions[step - 1], refSolutions[step - 1]) * 100;

    // Report results.
    info("ndof_coarse: %d, ndof_fine: %d, err_est_rel: %g%%",
      spaces[step - 1]->get_num_dofs(), refSpaces[step - 1]->get_num_dofs(), err_est_rel);

    adaptivity->adapt(selector, THRESHOLD, STRATEGY, MESH_REGULARITY);
    info("new space: %d", spaces[step]->get_num_dofs());

    // Clean up.
    delete [] coeff_vec;
    delete adaptivity;
}

int main(int argc, char* argv[])
{

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


  Selector<double>* selector = new H1ProjBasedSelector<double>(H2D_HP_ANISO, CONV_EXP, H2DRS_DEFAULT_ORDER);
  //Selector<double>* selector = new H1ProjBasedSelector<double>(H2D_H_ANISO, CONV_EXP, H2DRS_DEFAULT_ORDER);
  //Selector<double>* selector = new HOnlySelector<double>();

  // Create an H1 space with default shapeset.
  spaces[0] = new H1Space<double>(&mesh, &bcs, P_INIT);
  solutions[0] = new Solution<double>(&mesh);
  int ndof = spaces[0]->get_num_dofs();
  info("ndof = %d", ndof);

  for(int step = 1; step <= STEPS; step++)
  {
      adaptiveStep(step, &wf, selector);
  }

//  Hermes::Hermes2D::Views::ScalarView view("Solution", new Hermes::Hermes2D::Views::WinGeom(0, 0, 440, 350));
//  view.show(&sln, Hermes::Hermes2D::Views::HERMES_EPS_HIGH);
//  Hermes::Hermes2D::Views::View::wait();


  return 0;
}

