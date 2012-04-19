#define HERMES_REPORT_INFO
#include "hermes2d.h"

#include "definitions.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;

// This test makes sure that subdomains work correctly.

// Uniform polynomial degree of mesh elements.
const int P_INIT = 2;
// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 0;
// Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
Hermes::MatrixSolverType matrix_solver_type = Hermes::SOLVER_UMFPACK;
const double NEWTON_TOL = 1e-4;
// Maximum allowed number of Newton iterations.
const int NEWTON_MAX_ITER = 100;

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh_current, mesh_heat;
  Hermes::vector<Mesh*> meshes(&mesh_current, &mesh_heat);
  MeshReaderH2DXML mloader;
  mloader.load("problem-traverse.xml", meshes);

  // Perform initial mesh refinements (optional).
  for(int i = 0; i < INIT_REF_NUM; i++)
    for(unsigned int meshes_i = 0; meshes_i < meshes.size(); meshes_i++)
      meshes[meshes_i]->refine_all_elements();

  //  mloader.save("problem2.xml", meshes);
  //  mloader.load("problem2.xml", meshes);

  // Initialize essential boundary conditions.
  DefaultEssentialBCConst<double> bc_essential_current(Hermes::vector<std::string>("8", "9", "6", "3"), 0.0);
  EssentialBCs<double> bcs_current(&bc_essential_current);

  DefaultEssentialBCConst<double> bc_essential_heat(Hermes::vector<std::string>("2", "5", "1", "6", "4", "0", "7", "3"), 0.0);
  EssentialBCs<double> bcs_heat(&bc_essential_heat);

  //  DefaultEssentialBCConst<double> bc_essential_complement(Hermes::vector<std::string>("Bottom Right", "Top Right", "Top Left", "Horizontal Left", "Vertical Bottom"), 0.0);
  //  EssentialBCs<double> bcs_complement(&bc_essential_complement);

  // Create H1 spaces with default shapeset.
  H1Space<double> space_current(&mesh_current, &bcs_current, P_INIT);
  int ndof_current = space_current.get_num_dofs();

  H1Space<double> space_heat(&mesh_heat, &bcs_heat, P_INIT);
  int ndof_heat = space_heat.get_num_dofs();

  //  H1Space<double> space_complement(&mesh_complement, &bcs_complement, P_INIT);
  //  int ndof_complement = space_complement.get_num_dofs();

  info("ndofs current: %d\n", ndof_current);
  info("ndofs heat: %d\n", ndof_heat);

  //  Hermes::vector<Space<double> *> all_spaces(&space_current, &space_heat);
  //  Hermes::vector<const Space<double> *> all_spaces_const(&space_current, &space_heat);

  //  // Calculate and report the number of degrees of freedom.
  //  int ndof = Space<double>::get_num_dofs(all_spaces_const);
  //  info("ndof = %d.", ndof);

  //  double* coeff_vec = new double[ndof];

  // Initialize weak formulation.
  CustomWeakFormCurrent wf_current;

  // Initialize the FE problem.
  DiscreteProblem<double> dp_current(&wf_current, &space_current);

  // Initialize the Newton solver.
  NewtonSolver<double> newton_current(&dp_current, matrix_solver_type);

  // Perform Newton's iteration.
  info("Solving nonlinear problem:");
  bool verbose = true;
  // Perform Newton's iteration and translate the resulting coefficient vector into previous time level solutions.
  newton_current.set_verbose_output(verbose);
  try
  {
    newton_current.solve(NULL, NEWTON_TOL, NEWTON_MAX_ITER);
  }
  catch(Hermes::Exceptions::Exception e)
  {
    e.printMsg();
    error("Newton's iteration failed for current.");
  };

  Solution<double> solution_current;
  Solution<double>::vector_to_solution(newton_current.get_sln_vector(), &space_current, &solution_current);


  // Initialize weak formulation.
  CustomWeakFormHeat wf_heat(&solution_current);

  // Initialize the FE problem.
  DiscreteProblem<double> dp_heat(&wf_heat, &space_heat);

  // Initialize the Newton solver.
  NewtonSolver<double> newton_heat(&dp_heat, matrix_solver_type);

  // Perform Newton's iteration.
  info("Solving nonlinear problem:");
  // Perform Newton's iteration and translate the resulting coefficient vector into previous time level solutions.
  newton_heat.set_verbose_output(verbose);
  try
  {
    newton_heat.solve(NULL, NEWTON_TOL, NEWTON_MAX_ITER);
  }
  catch(Hermes::Exceptions::Exception e)
  {
    e.printMsg();
    error("Newton's iteration failed for heat.");
  };

  //      Solution<double> solution_heat;

  //    {
  //      Hermes::vector<Solution<double> *> tmp(&xvel_prev_time, &yvel_prev_time, &p_prev_time, &temperature_prev_time);
  //      Solution<double>::vector_to_solutions(newton.get_sln_vector(), Hermes::vector<const Space<double> *>(&xvel_space,
  //          &yvel_space, &p_space, &temperature_space), tmp);
  //    }

  //    // Show the solution at the end of time step.
  //    sprintf(title, "Velocity [m/s], time %g s", current_time);
  //    vview.set_title(title);
  //    vview.show(&xvel_prev_time, &yvel_prev_time);
  //    sprintf(title, "Pressure [Pa], time %g s", current_time);
  //    pview.set_title(title);
  //    pview.show(&p_prev_time);
  //    sprintf(title, "Temperature [C], time %g s", current_time);
  //    tempview.set_title(title);
  //    tempview.show(&temperature_prev_time, Views::HERMES_EPS_HIGH);
  //  }

  // Wait for all views to be closed.
  Views::View::wait();
  return 0;


  //  if(ndof_whole_domain == 9)
  //  {
  //    info("Success!");
  //    return TEST_SUCCESS;
  //  }
  //  else
  //  {
  //    info("Failure!");
  //    return TEST_FAILURE;
  //  }

  return 0;
}
